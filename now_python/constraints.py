import numpy as np
from optimizationProblem import NOW_config
from scipy.optimize import LinearConstraint, NonlinearConstraint
from scipy.sparse import csr_matrix, csc_matrix, lil_array, tril, diags, eye, vstack, hstack, kron
import torch


class constraints:
    def __init__(self, problem = NOW_config(), method=None):
        self.problem = problem
        self.method = method

        self.dt = self.problem._dt 

        self.eps = self.problem._tolIsotropy
        self.N = self.problem.N

        self.indices = self.problem.zeroGradientAtIndex
        self.signs = self.problem.signs

        # Targets
        self.Bhat = self.problem.targetTensor
        self.gMax = self.problem._gMaxConstraint
        self.sMax = self.problem._sMaxConstraint
        self.integralMax = self.problem._integralConstraint
        self.tolMaxwell = self.problem._tolMaxwell

        #Integration and derivative matrices
        self.Theta = build_integration_matrix(self.problem.N, 1)
        self.A1 = build_first_derivative_matrix(self.problem.N, 1)
        self.A2 = build_second_derivative_matrix(self.problem.N, 1)

    def unpack(self, x):
        q = x[:-1].reshape(-1, 3, order='F')
        b = x[-1] 
        return q, b

    def all_linear(self):
        #Zeroing on interval and boundary
        A = np.zeros((2 + len(self.indices), self.N))
        A[0,0] = 1 ; A[1,-1] = 1
        A[2:, :] = self.A1.toarray()[self.indices, :]
        A = np.kron(np.eye(3), A)
        A = np.hstack([A, np.zeros((A.shape[0], 1))])
        A_1 = csc_matrix(A) ; b_1 = np.zeros((A_1.shape[0]))

        #Max gradient strength
        if self.problem.useMaxNorm:
            A1 = self.A1.toarray() ; b = np.ones(int(self.N-1))*self.gMax
            A = np.kron(np.eye(3), A1)
            A = np.hstack([A, np.zeros((A.shape[0], 1))])
            A_2 = csc_matrix(A) ; b_2 = np.repeat(b,3)
        else:
            A_2 = csc_matrix(np.empty((0, self.N*3 + 1))) ; b_2 = []

        #Max slew rate
        A1 = self.A2.toarray() ; b = np.ones(int(self.N))*self.sMax
        A = np.kron(np.eye(3), A1)
        A = np.hstack([A, np.zeros((A.shape[0], 1))])
        A_3 = csc_matrix(A) ; b_3 = np.repeat(b,3)

        #Linear Motion compensation
        if np.any(self.problem.motionCompensation['linear']):
            t = (np.arange(1, self.N + 1) - 0.5)*self.dt
            linear_ind = np.where(self.problem.motionCompensation['linear'])[0]
            A = np.zeros((len(linear_ind), self.N))
            for i in range(len(linear_ind)):
                order = self.problem.motionCompensation['order'][linear_ind[i]]
                A[i,:] = - order * self.dt * t**(order - 1)
            A = np.kron(np.eye(3), A)
            A = np.hstack([A, np.zeros((A.shape[0], 1))])
            A_4 = csc_matrix(A) ; b_4 = np.zeros((A_4.shape[0]))
        else:
            A_4 = csc_matrix(np.empty((0, self.N*3 + 1))) ; b_4 = []

        #Stacking 
        A = vstack([A_1, A_2, A_3, A_4])
        #A = A.toarray() # used if using trust-constr as jac need to match (dense or sparse)
        b = np.concatenate([b_1, b_2, b_3, b_4])
        if self.method == "SLSQP": return [dict( type='ineq', fun=lambda x: b - A @ x ), dict( type='ineq', fun=lambda x: b + A @ x )]
        else: return [LinearConstraint(A, -b, b)]

    def all_nonlinear(self):
        def fun(x):
            q, b = self.unpack(x) ; N = self.N
            g = self.A1 @ q
            B = q.T @ (self.Theta/(N-1)) @ q

            #Tensor encoding   
            c1 = np.array([ (self.eps * b)**2 - np.trace((B-b*self.Bhat).T @ (B-b*self.Bhat)) ])  # np.array([(self.eps * b)**2 - np.linalg.norm(B - b * self.Bhat, 'fro')**2 ]) 
            
            #Max gradient strength
            if not self.problem.useMaxNorm:
                c2 = (self.gMax**2 - np.sum(g**2, axis = 1)).T
            else:
                c2 = np.array([], dtype=float)
            
            #Power constraint 
            c3 = np.array([self.integralMax - g[:,0].T @ g[:,0]])
            c4 = np.array([self.integralMax - g[:,1].T @ g[:,1]])
            c5 = np.array([self.integralMax - g[:,2].T @ g[:,2]])
            
            #Maxwell
            if self.problem.doMaxwellComp:
                M = g.T@(g*self.signs)
                m = np.sqrt(np.trace(M.T @ M)) 
                c6 = np.array([self.tolMaxwell - m])
            else:
                c6 = np.array([], dtype=float)

            #Non linear Motion compensation 
            if not np.all(self.problem.motionCompensation['linear']):
                nonlinear_ind = np.where(~self.problem.motionCompensation["linear"])[0] 
                t = (np.arange(1, self.N + 1) - 0.5)*self.dt
                gamma = 2.6751e+08 # radians / T / s for hydrogen.

                A = np.zeros((1,len(nonlinear_ind)))
                for i in range(len(nonlinear_ind)):
                    order = self.problem.motionCompensation['order'][nonlinear_ind[i]]
                    moment_weighting = - order*self.dt*t**(order - 1)
                    #print(moment_weighting.shape, q.shape)
                    moment_vector = moment_weighting @ q
                    A[0,i] = (self.problem.motionCompensation['maxMagnitude'][nonlinear_ind[i]]*1000**order / (gamma * 1e-6))**2 - np.sum(moment_vector**2) 
                c7 = A.ravel()
            else:
                c7 = np.array([], dtype=float)
            
            return np.hstack([c1,c2,c3,c4,c5,c6,c7]) 
        def jac(x): #(m,n)
            q, b = self.unpack(x) ; N = self.N
            g = self.A1 @ q
            wq = (self.Theta/(N-1)) @ q
            B = q.T @ wq

            #Tensor encoding   
            dc1_db = 2*B - 2*b*self.Bhat
            firstTerm = np.kron(np.eye(3),wq.T)
            secondTerm = np.reshape(firstTerm, (3,3,3*N), order='F')
            secondTerm = np.transpose(secondTerm, (1,0,2))
            secondTerm = np.reshape(secondTerm, (9,3*N))
            dB_dq = firstTerm + secondTerm
            dc1_dq = np.reshape(dc1_db, (1,9))@dB_dq
            dc1_ds = np.array([ 2*b*np.trace(self.Bhat.T@self.Bhat) - 2*np.trace(B.T@self.Bhat) - 2*b*self.eps**2 ]).reshape(1,1)
            dc1_dx = - np.hstack([dc1_dq, dc1_ds])
            
            #Max gradient strength
            if not self.problem.useMaxNorm:
                #dc2_dq = 2*np.vstack([self.A1 @ g[:,0], self.A1 @ g[:,1], self.A1 @ g[:,2]])
                dc2_dq = 2 * np.hstack([
                    self.A1.toarray() * g[:N-1, 0][:, None],   # g[0:N-1, 0] -> (N-1,1)
                    self.A1.toarray() * g[:N-1, 1][:, None],
                    self.A1.toarray() * g[:N-1, 2][:, None],
                ])

                """
                dc2_dq = 2*np.vstack([diags(g[:,0], shape = (N-1, N-1))@self.A1, diags(g[:,1], shape = (N-1, N-1))@self.A1 , diags(g[:,2], shape = (N-1, N-1))@self.A1 ])
                test= diags([g[:,0]], shape = (N-1, N-1))@self.A1
                print(dc2_dq.shape, test.shape, "wtf")
                dc2_dx = - np.hstack([dc2_dq.T, np.zeros((N-1, 1))])
                print(dc2_dq.shape, dc2_dx.shape, test.shape, "wtf")
                """
                dc2_dx = - np.hstack([dc2_dq, np.zeros((N-1, 1))])
            else:
                dc2_dx = np.empty((0,3*N + 1))
            
            #Power constraint 
            dc3_dx = np.zeros((1,3*N+1))
            dc3_dx[0,0:N] = -2*g[:,0].T @ self.A1
            dc4_dx = np.zeros((1,3*N+1))
            dc4_dx[0,N:2*N] = -2*g[:,1].T @ self.A1
            dc5_dx = np.zeros((1,3*N+1))
            dc5_dx[0,2*N:-1] = -2*g[:,2].T @ self.A1
            
            #Maxwell
            if self.problem.doMaxwellComp:
                signedg = g*self.signs
                M = g.T@signedg
                m = np.sqrt(np.trace(M.T@M))
                dc6_dM = 1/m * M
                firstTerm = np.kron(np.eye(3), (self.A1.T @ signedg ).T )
                secondTerm = np.reshape(firstTerm, (3,3,3*N), order='F')
                secondTerm = np.transpose(secondTerm, (1,0,2))
                secondTerm = np.reshape(secondTerm, (9,3*N))
                dM_dq = firstTerm + secondTerm
                dc6_dq = np.reshape(dc6_dM, (1,9)) @ dM_dq
                #print(dc6_dq.shape)
                dc6_dx = - np.hstack([dc6_dq, np.zeros((1, 1))]) 

            else:
                dc6_dx = np.empty((0,3*N + 1))

            #Non linear Motion compensation
            if not np.all(self.problem.motionCompensation['linear']):
                nonlinear_ind = np.where(~self.problem.motionCompensation["linear"])[0] 
                t = (np.arange(1, self.N + 1) - 0.5)*self.dt
                gamma = 2.6751e+08 # radians / T / s for hydrogen.

                A = np.zeros((len(nonlinear_ind),3*N+1))
                for i in range(len(nonlinear_ind)):
                    order = self.problem.motionCompensation['order'][nonlinear_ind[i]]
                    moment_weighting = - order*self.dt*t**(order - 1)
                    moment_vector = moment_weighting @ q
                    A[i,:-1] = -2 * np.kron(moment_vector, moment_weighting)
                dc7_dx = A
            else:
                dc7_dx = np.empty((0,3*N + 1))
    
            return csc_matrix(np.vstack([dc1_dx,dc2_dx,dc3_dx,dc4_dx,dc5_dx,dc6_dx,dc7_dx]))


        if self.method == "SLSQP": return [dict( type='ineq', fun=lambda x: fun(x))]
        else: return [NonlinearConstraint(fun, 0, np.inf, jac=jac)]
    
    def build_contraints(self):
        constraints = []
       
        constraints += self.all_linear()
        constraints += self.all_nonlinear()
        
        return constraints


#Finite difference matrices

def build_integration_matrix(N: int, dt: float) -> csc_matrix:
    """Sparse integration matrix Θ (trapezoidal cumulative integration)."""
    diagonal = np.ones(N)
    diagonal[0] =  0.5   ; diagonal[-1] = 0.5
    Theta = diags(diagonal, shape=(N, N), format='csc') * dt
    return Theta

def build_first_derivative_matrix(N: int, dt: float) -> csc_matrix:
    """Sparse first derivative matrix A₁ using backward difference."""
    diagonals = [-np.ones(N-1), np.ones(N-1)]
    offsets = [0, 1]
    A1 = diags(diagonals, offsets, shape=(N-1, N), format='csc') / dt
    return A1

def build_second_derivative_matrix(N: int, dt: float) -> csc_matrix:
    """Sparse second derivative matrix A₂ using central differences."""
    diagonals = [np.ones(N-1), -2*np.ones(N), np.ones(N-1)]
    offsets = [-1, 0, 1]
    A2 = diags(diagonals, offsets, shape=(N, N), format='csc') / (dt ** 2)
    return A2.tocsc()


if __name__ == "__main__":
    pass



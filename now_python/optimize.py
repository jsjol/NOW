from optimizationProblem import NOW_config, obj_fun
from scipy.optimize import minimize
from constraints import constraints
import numpy as np
import matplotlib.pyplot as plt
import time
import os, json, sys, platform
from datetime import datetime
from types import SimpleNamespace
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"


class optimize:
    def __init__(self, 
                 problem = NOW_config(),
                 ):
        
        self.res = None #holds the result after run()
        self.q, self.g, self.b, self.B = None, None, None, None #results
        
        self.problem = problem 
        self.objective_function = obj_fun

        constrains = constraints(problem=self.problem) #Creates the class which specifies the contraints given the problem, specifies if method to compute jacobian is standard or autograd
        self.constraints = constrains.build_contraints() #Builds the constraints in a format accepted by scipy.optimize.minimize (list of dicts or constraint objects)
        self.Theta, self.A1, self.A2 = constrains.Theta, constrains.A1, constrains.A2 #Integration and derivative matrices
        self.dt = self.problem._dt

        self.intermediate_results = [] #store intermediate results for plotting

        self.elapsed_time = None #Timing of the optimization process
        self.method_used = None #The optimization method used

    def callback(self, xk):
        self.intermediate_results.append(xk.copy())
    
    def callback_trust(self, xk, optimize_result): #callback function for trust-constr method
        self.intermediate_results.append(xk.copy())

    def run(self, method='SLSQP', early_stopping = True):
        t0 = time.perf_counter()
        self.method_used = method
        if method == 'SLSQP':
            if early_stopping:
                chosen_callback = self.callback_with_earlystop
            else:
                chosen_callback = self.callback
            try: 
                self.res = minimize(
                    fun = self.objective_function, 
                    x0 = self.problem.x0, 
                    jac = True, #objfun provides gradient as tuple !! 
                    method = 'SLSQP', 
                    constraints = self.constraints,
                    callback=chosen_callback,
                    options={'maxiter': 2000,
                            'disp': True},
                            #'iprint': 3}, # according to ChatGPT: 0	No output at all (default) 1	Final result summary only (same as disp=True) 2	Iteration-by-iteration summary of the optimization progress 3 or higher	Even more internal details per iteration (gradient, step size, constraint info, etc.)
                    tol=1e-6,
                    #callback=self.callback
                    )
            except StopIteration:
                print("Optimization stopped early.")
                # SLSQP does not return a result when StopIteration is thrown
                # so we must construct the result manually
                self.res = SimpleNamespace(x=self.intermediate_results[-1])
        elif method == 'trust-constr':
            self.res = minimize(
                fun = self.objective_function, 
                x0 = self.problem.x0, 
                jac = True, #objfun provides gradient as tuple !! 
                method = 'trust-constr', 
                constraints = self.constraints,
                options={'maxiter': 2000,
                         'disp': True,
                         'verbose':3},
                callback=self.callback_trust,
                         
            )
        self.elapsed_time = time.perf_counter() - t0
        
    def get_results(self):
        now_gamma = 42.576e6  # Hz/T, samma som now_gamma i MATLAB
        gamma = now_gamma*2*np.pi

        q_raw = self.res.x[:-1].reshape(self.problem.N, 3, order="F") #What is used in the optimization, same as Matlab
        q_phys = q_raw * gamma * 1e-6 #SI-units
        g = (self.A1 / self.dt) @ q_raw #OBS OBS OBS ask why this is q_raw and not q_phys??? In matlab code, they use q_raw for some reason
        slew = (self.A2 @ q_raw) / (self.dt**2) #OBS OBS OBS here aswell, q raw??????
        dt_s = self.dt * 1e-3
        B = dt_s * (q_phys.T @ q_phys)
        b_raw = self.res.x[-1] #OBS only corresponds to the intern variable, but we actually want to calculate b ourselves
        b = np.trace(B) * 1e-6 # s/mm^2, This is the exact b-value per definition, and also what is done in the matlab code

        self.q = q_phys
        self.g = g
        self.B = B
        self.b = b
        self.slew = slew
        
        return self.q, self.g, self.b, self.B
    
    def results_from_x(self, x):
        now_gamma = 42.576e6  # Hz/T, samma som now_gamma i MATLAB
        gamma = now_gamma*2*np.pi

        q_raw = x[:-1].reshape(self.problem.N, 3, order="F") #What is used in the optimization, same as Matlab
        q_phys = q_raw * gamma * 1e-6 #SI-units
        g = (self.A1 / self.dt) @ q_raw #OBS OBS OBS ask why this is q_raw and not q_phys??? In matlab code, they use q_raw for some reason
        slew = (self.A2 @ q_raw) / (self.dt**2) #OBS OBS OBS here aswell, q raw??????
        dt_s = self.dt * 1e-3
        B = dt_s * (q_phys.T @ q_phys)
        #b_raw = self.res.x[-1] #OBS only corresponds to the intern variable, but we actually want to calculate b ourselves
        b = np.trace(B) * 1e-6 # s/mm^2, This is the exact b-value per definition, and also what is done in the matlab code

        return q_phys, g, b, B
    
    def plot_results(self, animate = False, dir=None, x =None):
        if x is not None:
            q, g, b, B = self.results_from_x(x)
            self.res = SimpleNamespace(x=self.intermediate_results[-1],nit = len(self.intermediate_results))
        else:
            q, g, b, B = self.get_results()

        t_g = np.linspace(0.5 * self.problem._dt, self.problem.totalTimeActual - 0.5 * self.problem._dt, endpoint=True, num=self.problem.N - 1)
        t_q = np.linspace(0, self.problem.totalTimeActual, self.problem.N, endpoint = True)

        fig, axs = plt.subplots(2, 1, sharex = True, figsize = (8, 6))
        q_mikrometer = q/10**6
        #subplot for q
        qx = axs[0].plot(t_q, q_mikrometer[:, 0], label='q_x')
        qy = axs[0].plot(t_q, q_mikrometer[:, 1], label='q_y')
        qz = axs[0].plot(t_q, q_mikrometer[:, 2], label='q_z')
        axs[0].set_ylabel("q [µs/m]")
        axs[0].set_title("Integrated gradient q(t)")
        axs[0].legend()
        axs[0].grid(True)

        #subplot for g
        gx = axs[1].plot(t_g, g[:, 0], label='g_x')
        gy = axs[1].plot(t_g, g[:, 1], label='g_y')
        gz = axs[1].plot(t_g, g[:, 2], label='g_z')
        axs[1].set_xlabel("Time [ms]")
        axs[1].set_ylabel("g [mT/m]")
        axs[1].set_title("Magnetic field gradient g(t)")
        axs[1].legend()
        axs[1].grid(True)

        #Information text
        iter_count = float(len(self.intermediate_results)) #getattr(self.res, "nit", None)
        elapsed_time = getattr(self, "elapsed_time", None)
        target = self.problem.targetTensor          # 3x3
        trace_B = np.trace(B)
        trace_target = np.trace(target)
        encoding_tensor = (B / trace_B) * trace_target
        fro = np.linalg.norm(encoding_tensor - target)
        
        info_text = (
            f"Iterations: {iter_count}\n"
            f"Time: {elapsed_time:.3f} s\n"
            f"b: {b:.6f} s/mm^2\n"
            f"Fro-norm = {fro:.6e}\n"
            "Encoding tensor:\n"
            + "\n".join(["  ".join([f"{val: .6f}" for val in row]) for row in encoding_tensor])
        )


        if animate == False:
            plt.tight_layout(rect=[0, 0.22, 1, 1])  #Leave room for the box 
            info_ax = fig.add_axes([0.08, 0.03, 0.84, 0.16])  #[left, bottom, width, height] (distances)
            info_ax.axis('off')
            info_ax.text(0.5, 0.5, info_text,
                        ha='center', va='center', fontsize=8,
                        bbox=dict(facecolor='white', alpha=0.85, boxstyle='round'))
            if dir == None:
                plt.show()
            else:
                plt.savefig(dir)

        elif animate == True:
            from matplotlib.widgets import Slider, Button
            
            plt.tight_layout(rect=[0, 0.22, 1, 1])  #Leave room for the box 
            info_ax = fig.add_axes([0.08, 0.07, 0.84, 0.20])  #[left, bottom, width, height] (distances)
            info_ax.axis('off')
            info_ax.text(0.5, 0.5, info_text,
                        ha='center', va='center', fontsize=8,
                        bbox=dict(facecolor='white', alpha=0.85, boxstyle='round'))

            q_frames = [] ; g_frames = [] ; b_frames = [] ; B_frames = []
            intermediate_results = self.intermediate_results
            for x in intermediate_results:
                _q, _g, _b, _B = self.results_from_x(x)
                q_frames.append(_q)
                g_frames.append(_g)
                b_frames.append(_b)
                B_frames.append(_B)

            # Get line objects for easy updating
            lines_q = axs[0].lines  # [q_x, q_y, q_z]
            lines_g = axs[1].lines  # [g_x, g_y, g_z]

            n_frames = len(q_frames)

            # Make space at bottom for the slider
            plt.subplots_adjust(bottom=0.35)

            # Slider axis: [left, bottom, width, height] in figure coordinates
            ax_slider = plt.axes([0.15, 0.05, 0.7, 0.03])
            frame_slider = Slider(
                ax=ax_slider,
                label="Iteration",
                valmin=0,
                valmax=n_frames - 1,
                valinit=0,
                valstep=1,           # integer steps = frame indices
            )

            def update_slider(val):
                i = int(val)

                q = q_frames[i]
                g = g_frames[i]
                b = b_frames[i]
                B = B_frames[i]

                # update q lines
                for j in range(3):
                    lines_q[j].set_ydata(q[:, j]/10**6) #mikrometer conversion

                # update g lines
                for j in range(3):
                    lines_g[j].set_ydata(g[:, j])

                #Information text
                elapsed_time = getattr(self, "elapsed_time", None)
                target = self.problem.targetTensor          # 3x3
                trace_B = np.trace(B)
                trace_target = np.trace(target)
                encoding_tensor = (B / trace_B) * trace_target
                fro = np.linalg.norm(encoding_tensor - target)

                info_text = (
                    f"Iteration: {i}\n"
                    f"Time: {elapsed_time:.3f} s\n"
                    f"b: {b:.6f} s/mm^2\n"
                    f"Fro-norm = {fro:.6e}\n"
                    "Encoding tensor:\n"
                    + "\n".join(["  ".join([f"{val: .6f}" for val in row]) for row in encoding_tensor])
                )
                info_ax.text(0.5, 0.5, info_text,
                        ha='center', va='center', fontsize=8,
                        bbox=dict(facecolor='white', alpha=0.85, boxstyle='round'))
                

                fig.canvas.draw_idle()  # redraw efficiently

            # Connect the slider to the update function
            frame_slider.on_changed(update_slider)

            plt.show()

    def callback_with_earlystop(self, xk):

        # Save xk for animation
        self.callback(xk)

        # Step size
        if len(self.intermediate_results) > 1:
            step_norm = np.linalg.norm((xk - self.intermediate_results[-2])[:-1]) #OBS OBS OBS IGNORERA b-värdet!!!
        else:
            step_norm = np.inf

        # Objective
        f, grad = self.objective_function(xk)

        # Constraint violation
        violations = []
        for c in self.constraints:
            if isinstance(c, dict):
                val = c["fun"](xk)
                violations.append(np.min(val))
        constraint_violation = max(0, -min(violations)) if violations else 0

        # External stagnation
        if not hasattr(self, "prev_f_ext"):
            self.prev_f_ext = f
            self.no_improve_ext = 0
        else:
            if abs(f - self.prev_f_ext) < 1e-4:   # <-- BIGGER tolerance!
                self.no_improve_ext += 1
            else:
                self.no_improve_ext = 0

            self.prev_f_ext = f

        # Early stop for EXTERNAL iterations
        #print(
        #    f"[ITER {len(self.intermediate_results)-1:4d}]  "
        #    f"Δf = {abs(f - self.prev_f_ext):.3e}   "
        #    f"step = {step_norm:.3e}   "
        #    f"viol = {constraint_violation:.3e}   "
        #    f"no_improve = {self.no_improve_ext}"
        #)
        if (
            constraint_violation < 1e-6 and   # constraints ok
            step_norm < 1e-2 and              # small enough for external loop
            self.no_improve_ext >= 20         # no improvement for 10 external iterations
        ):
            print("EARLY STOPPING (EXTERNAL): Stagnation detected.")
            raise StopIteration



if __name__ == "__main__":
    pass





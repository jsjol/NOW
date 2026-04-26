import numpy as np
import matplotlib.pyplot as plt


def plot_result(result, config, show=True, save_path=None):
    """Plot gradient waveform and q-trajectory."""
    t_g = np.linspace(0.5 * config._dt,
                      config.totalTimeActual - 0.5 * config._dt,
                      num=config.N - 1)
    t_q = np.linspace(0, config.totalTimeActual, config.N)

    q_um = result.q / 1e6  # convert to µm⁻¹
    g = result.g[1:-1, :]  # strip zero padding

    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    axs[0].plot(t_q, q_um[:, 0], label='$q_x$')
    axs[0].plot(t_q, q_um[:, 1], label='$q_y$')
    axs[0].plot(t_q, q_um[:, 2], label='$q_z$')
    axs[0].set_ylabel(r'q [$\mu$m$^{-1}$]')
    axs[0].set_title('Integrated gradient q(t)')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(t_g, g[:, 0], label='$g_x$')
    axs[1].plot(t_g, g[:, 1], label='$g_y$')
    axs[1].plot(t_g, g[:, 2], label='$g_z$')
    axs[1].set_xlabel('Time [ms]')
    axs[1].set_ylabel('g [mT/m]')
    axs[1].set_title('Magnetic field gradient g(t)')
    axs[1].legend()
    axs[1].grid(True)

    B = result.B
    trace_B = np.trace(B)
    trace_target = np.trace(config.targetTensor)
    if trace_B > 0 and trace_target > 0:
        encoding_tensor = (B / trace_B) * trace_target
        fro = np.linalg.norm(encoding_tensor - config.targetTensor)
    else:
        encoding_tensor = B
        fro = float('nan')

    info_text = (
        f"b = {result.b:.4f} s/mm²\n"
        f"κ = {result.kappa:.4f}\n"
        f"||B̂ - B_target||_F = {fro:.2e}\n"
        f"Time: {result.optimizationTime:.2f}s, {result.iter} attempt(s)"
    )

    plt.tight_layout(rect=[0, 0.12, 1, 1])
    info_ax = fig.add_axes([0.08, 0.01, 0.84, 0.10])
    info_ax.axis('off')
    info_ax.text(0.5, 0.5, info_text, ha='center', va='center', fontsize=9,
                 bbox=dict(facecolor='white', alpha=0.85, boxstyle='round'))

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    if show:
        plt.show()

    return fig

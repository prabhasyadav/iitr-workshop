import numpy as np
import matplotlib.pyplot as plt


def groundwater_flow_with_well(
    nx=50, ny=50, Lx=1000.0, Ly=1000.0,
    K=10.0,              # hydraulic conductivity [m/day], assumed isotropic
    Qw=-500.0,           # well rate [m^3/day] (negative = pumping)
    wx=500.0, wy=500.0,  # well coordinates [m]
    h_left=100.0, h_right=90.0,  # Dirichlet heads [m]
    max_iter=8000, tol=1e-6
):
    """
    2-D steady-state confined groundwater flow with a single well (finite differences).
    Left/right: Dirichlet heads; Top/bottom: no-flow (Neumann).

    Returns:
        h  : (ny, nx) head array [m]
        x,y: 1D coordinates [m]
        (ix,iy): indices of the well cell
        qx,qy: Darcy specific discharge components [m/day]
    """
    # Grid spacing
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)

    # Grid coordinates
    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, Ly, ny)

    # Head field
    h = np.zeros((ny, nx))

    # Dirichlet boundaries (left/right)
    h[:, 0] = h_left
    h[:, -1] = h_right

    # Map well location to nearest grid cell (clip to be safe)
    ix = int(np.clip(round(wx / dx), 0, nx - 1))
    iy = int(np.clip(round(wy / dy), 0, ny - 1))

    # Source/sink per cell area (m/day) / (m^2) = m/day/m^2 = 1/day * m → but in discrete form
    # Add as S term in the FD update at the well cell.
    cell_area = dx * dy
    S_term = Qw / (cell_area)   # [m/day] * (1/m^2) → m/day/m^2 (volumetric per area)

    # Iterative solve (Gauss–Seidel-like Jacobi update for Laplace with source at well)
    for it in range(max_iter):
        h_old = h.copy()

        # Sweep interior points; top/bottom rows included with “mirrored” neighbors via index clamping below
        for j in range(0, ny):
            jm = max(j - 1, 0)         # mirror at bottom boundary
            jp = min(j + 1, ny - 1)    # mirror at top boundary
            for i in range(1, nx - 1):  # exclude left/right Dirichlet columns
                # Source at the well cell only
                S = (S_term / K) if (i == ix and j == iy) else 0.0  # divide by K to match your original formulation

                h[j, i] = 0.5 * (
                    (h[j, i + 1] + h[j, i - 1]) / dx**2 +
                    (h[jp, i]   + h[jm, i])     / dy**2 - S
                ) / (1.0/dx**2 + 1.0/dy**2)

        # Re-enforce Dirichlet after each sweep
        h[:, 0]  = h_left
        h[:, -1] = h_right

        error = np.max(np.abs(h - h_old))
        if error < tol:
            print(f"Converged in {it} iterations (max Δh = {error:.2e})")
            break
    else:
        print("Did not converge within max_iter")

    # ---- Compute Darcy specific discharge q = (-K dh/dx, -K dh/dy) in m/day ----
    qx = np.zeros_like(h)
    qy = np.zeros_like(h)

    # dh/dx (central differences interior; one-sided at edges)
    qx[:, 1:-1] = -K * (h[:, 2:] - h[:, :-2]) / (2 * dx)
    qx[:, 0]    = -K * (h[:, 1] - h[:, 0]) / dx
    qx[:, -1]   = -K * (h[:, -1] - h[:, -2]) / dx

    # dh/dy (central differences interior; one-sided at edges)
    qy[1:-1, :] = -K * (h[2:, :] - h[:-2, :]) / (2 * dy)
    qy[0,    :] = -K * (h[1, :] - h[0,   :]) / dy
    qy[-1,   :] = -K * (h[-1, :] - h[-2,  :]) / dy

    return h, x, y, (ix, iy), qx, qy


# ---- Run model ----
h, x, y, (ix, iy), qx, qy = groundwater_flow_with_well()

# ---- Plot: head contours (head lines), streamlines, and quiver ----
X, Y = np.meshgrid(x, y, indexing='xy')  # (ny, nx) alignment with h

plt.figure(figsize=(10, 7.5))

# Head contours (equipotential lines)
cs = plt.contour(X, Y, h, levels=20, colors='k', linewidths=0.9)
plt.clabel(cs, inline=True, fmt="%.2f", fontsize=8)

# Filled head background (optional for context)
cf = plt.contourf(X, Y, h, levels=30, cmap='viridis', alpha=0.85)
plt.colorbar(cf, label='Hydraulic Head (m)')

# Streamlines (instantaneous flow lines) from Darcy q
# Note: streamplot expects U=dx/dt, V=dy/dt arrays with shape (ny, nx)
plt.streamplot(X, Y, qx, qy, color='white', density=1.2, linewidth=1.1, arrowsize=1.2)

# Quiver: sparse arrows for readability
step = 4  # take every Nth vector in both directions
Xs = X[::step, ::step]
Ys = Y[::step, ::step]
Us = qx[::step, ::step]
Vs = qy[::step, ::step]

# Scale quiver by magnitude percentile for nice arrow lengths
mag = np.hypot(Us, Vs)
p95 = np.percentile(mag, 95) if np.any(np.isfinite(mag)) else 1.0
scale = 0.08 * p95 if p95 > 0 else 1.0  # tweak 0.05–0.15 as needed

plt.quiver(Xs, Ys, Us, Vs, color='black', alpha=0.75,
           scale_units='xy', scale=scale, width=0.003)

# Well marker
plt.scatter(x[ix], y[iy], c=('crimson'), s=120, edgecolor='k', zorder=3, label='Pumping well')

plt.title("2-D Steady Groundwater Flow with Pumping Well\nHead lines, Streamlines, and Darcy Quiver")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.legend(loc='lower right')
plt.axis('equal')
plt.tight_layout()
plt.savefig("fig2D_Model.pdf", dpi=300)
plt.show()

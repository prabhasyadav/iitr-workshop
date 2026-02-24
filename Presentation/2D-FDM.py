import numpy as np
import matplotlib.pyplot as plt


def groundwater_flow_with_well(nx=50, ny=50, Lx=1000, Ly=1000,
                               K=10, Qw=-500, wx=500, wy=500):
    
    """
    2-D steady-state groundwater flow with a pumping well.

    PDE:
        ∂/∂x (K ∂h/∂x) + ∂/∂y (K ∂h/∂y) + Q/K = 0

    Parameters:
        nx, ny  : grid cells
        Lx, Ly  : domain length (m)
        K       : hydraulic conductivity (m/day)
        Qw      : pumping rate (m³/day), negative for pumping
        wx, wy  : well location (m)
    """
    #Grid spacing
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)

    # Grid
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)

    # Head array
    h = np.zeros((ny, nx))

    # Dirichlet BCs
    h[:, 0] = 100         # left boundary
    h[:, -1] = 90         # right boundary

    # Map well coordinates to nearest grid location
    ix = int(wx / dx)
    iy = int(wy / dy)

    # Source/sink term (Q per cell area)
    cell_area = dx * dy
    Q_term = Qw / (K * cell_area)

    max_iter = 8000
    tol = 1e-6

    for it in range(max_iter):
        h_old = h.copy()

        for j in range(1, ny - 1):
            for i in range(1, nx - 1):

                # Source/sink term only at well cell
                S = Q_term if (i == ix and j == iy) else 0.0

                h[j, i] = 0.5 * (
                    (h[j, i+1] + h[j, i-1]) / dx**2 +
                    (h[j+1, i] + h[j-1, i]) / dy**2 - S
                ) / (1/dx**2 + 1/dy**2)

        # Convergence check
        error = np.max(np.abs(h - h_old))
        if error < tol:
            print(f"Converged in {it} iterations with error {error:.2e}")
            break
    else:
        print("Did not converge")

    return h, x, y, (ix, iy)


# Run model
h, x, y, (ix, iy) = groundwater_flow_with_well()

# Plot results
plt.figure(figsize=(8, 6))
plt.contourf(x, y, h, 25, cmap='viridis')
plt.colorbar(label='Hydraulic Head (m)')
plt.scatter(x[ix], y[iy], c='red', s=80, label='Pumping Well')
plt.title("2-D Groundwater Flow With Pumping Well")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.legend()
plt.show()
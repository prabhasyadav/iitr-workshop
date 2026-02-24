import numpy as np
import matplotlib.pyplot as plt
from typing import Literal

# code obtained using MS 365 copilot, mdfied to fit the workshop context. 

def steady_1d_with_well(
    L: float = 1000.0,        # domain length [m]
    nx: int = 401,            # number of nodes (refine to better localize the well)
    K: float = 1e-4,          # hydraulic conductivity [m/s]
    b: float = 10.0,          # aquifer thickness [m]
    h0: float = 100.0,        # left boundary head [m]
    hL: float = 95.0,         # right boundary head [m]
    xw: float = 500.0,        # well location [m] along x
    Qw: float = -600.0,       # well rate (negative = pumping)
    Q_units: Literal["m3/s", "m3/day"] = "m3/day",
    width: float = 1.0        # out-of-plane width [m] (1D model thickness in y)
):
    """
    Steady-state 1-D confined groundwater flow with a single well (finite differences).

    Governing equation (steady, confined, no recharge):
        d/dx( T * dh/dx ) + S(x) = 0,   T = K*b
    where S(x) is a distributed source/sink term [m/s].

    Well representation:
        The point well of volumetric rate Qw [m^3/s] is applied to the nearest cell
        as an equivalent sink per unit plan area:
            S_well_cell = -(Qw) / (width * b * dx)   [m/s]
        (negative Qw → pumping → negative S added to RHS as -S)

    Boundary conditions:
        Dirichlet: h(0) = h0,  h(L) = hL

    Returns:
        x : node coordinates [m]
        h : head [m]
        q : Darcy flux [m/s] at nodes (centered interior, one-sided at boundaries)
        iw: index of the well host cell
    """
    # Grid
    x = np.linspace(0.0, L, nx)
    dx = x[1] - x[0]
    T = K * b

    # Map well to nearest node index
    iw = int(np.clip(np.round(xw / dx), 0, nx - 1))

    # Convert Q units to m3/s
    Qw_m3s = Qw / 86400.0 if Q_units.lower() in ["m3/day", "m^3/day", "m3/d", "m3d"] else Qw

    # Build linear system A h = rhs
    # Interior nodes: T*(h_{i+1} - 2h_i + h_{i-1})/dx^2 + S_i = 0
    # → -T/dx^2*h_{i-1} + 2T/dx^2*h_i - T/dx^2*h_{i+1} = -S_i
    A = np.zeros((nx, nx), dtype=float)
    rhs = np.zeros(nx, dtype=float)

    coef = T / dx**2
    for i in range(1, nx - 1):
        A[i, i - 1] = -coef
        A[i, i]     =  2 * coef
        A[i, i + 1] = -coef
        rhs[i]      =  0.0  # no recharge baseline

    # Add well sink into RHS at its host node:
    # S_i [m/s] multiplied by (-1) to RHS as above.
    # S_well_cell is volumetric per-area rate applied to the cell that spans dx in 1D:
    # Volume balance for that cell: Qw_m3s distributed over plan area (width * dx) and thickness b ⇒
    # Equivalent per-area source term (m/s) = Qw / (width * dx * b).
    # Sign: pumping Qw<0 → S_well_cell<0.
    S_well_cell = Qw_m3s / (width * dx * b)  # [m/s]
    rhs[iw] += -S_well_cell  # move +S to LHS → RHS gets -S

    # Dirichlet boundaries
    A[0, 0]   = 1.0
    rhs[0]    = h0
    A[-1, -1] = 1.0
    rhs[-1]   = hL

    # Solve linear system
    h = np.linalg.solve(A, rhs)

    # Darcy flux q = -K * dh/dx (node-based approximation)
    q = np.zeros_like(h)
    q[1:-1] = -K * (h[2:] - h[:-2]) / (2 * dx)
    q[0]    = -K * (h[1]  - h[0]) / dx
    q[-1]   = -K * (h[-1] - h[-2]) / dx

    return x, h, q, iw


if __name__ == "__main__":
    # Run an example
    x, h, q, iw = steady_1d_with_well(
        L=1000.0, nx=401,
        K=1e-4, b=10.0,
        h0=100.0, hL=95.0,
        xw=500.0, Qw=-600.0, Q_units="m3/day",
        width=1.0
    )

    # Plots
    plt.figure(figsize=(9, 4.4))
    plt.plot(x, h, "b-", lw=2, label="Head")
    plt.axvline(x[iw], color="crimson", ls="--", lw=1.5, label="Well location")
    plt.xlabel("x (m)"); plt.ylabel("Head h (m)")
    plt.title("1‑D Steady Confined Groundwater Flow with Pumping Well")
    plt.grid(True, alpha=0.3); plt.legend(); plt.tight_layout()
    plt.show()

    plt.figure(figsize=(9, 4.0))
    plt.plot(x, q * 86400.0, "m-", lw=2)  # m/s → m/day
    plt.axvline(x[iw], color="crimson", ls="--", lw=1.2)
    plt.xlabel("x (m)"); plt.ylabel("Darcy flux q (m/day)")
    plt.title("Darcy Flux Distribution")
    plt.grid(True, alpha=0.3); plt.tight_layout()
    plt.show()
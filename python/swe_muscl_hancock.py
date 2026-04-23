"""
1D Shallow Water Equations solver — MUSCL-Hancock + HLLC Riemann solver.

Port of the MATLAB case_1 algorithm.  All arrays are 0-indexed numpy arrays
of length n+4 (n interior cells + 2 ghost cells on each side).
Interior cells occupy indices 2:n+2; flux interfaces run over 2:n+3.
"""

import numpy as np

GRAV = 9.806
DRY_TOL = 1.0e-6


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def minmod(a, b):
    """Element-wise minmod slope limiter (works on scalars and arrays)."""
    return np.where(a * b <= 0, 0.0,
           np.where(np.abs(a) <= np.abs(b), a, b))


def _hllc(hl, hr, hul, hur):
    """
    HLLC flux at one interface.  hl/hr are well-balanced depths; hul/hur are
    the reconstructed discharges (not rescaled after depth correction).
    Returns (f1, f2): mass flux and x-momentum flux.
    """
    ul = hul / hl if hl > DRY_TOL else 0.0
    ur = hur / hr if hr > DRY_TOL else 0.0

    um = 0.5*(ul + ur) + np.sqrt(GRAV*hl) - np.sqrt(GRAV*hr)
    hm = (0.5*(np.sqrt(GRAV*hl) + np.sqrt(GRAV*hr)) + 0.25*(ul - ur))**2 / GRAV

    Sl = (ur - 2*np.sqrt(GRAV*hr)) if hl == 0.0 else min(ul - np.sqrt(GRAV*hl), um - np.sqrt(GRAV*hm))
    Sr = (ul + 2*np.sqrt(GRAV*hl)) if hr == 0.0 else max(ur + np.sqrt(GRAV*hr), um + np.sqrt(GRAV*hm))

    f1l = hl * ul;  f2l = hl*ul**2 + 0.5*GRAV*hl**2
    f1r = hr * ur;  f2r = hr*ur**2 + 0.5*GRAV*hr**2

    if Sr != Sl:
        f1m = (Sr*f1l - Sl*f1r + Sl*Sr*(hr - hl))     / (Sr - Sl)
        f2m = (Sr*f2l - Sl*f2r + Sl*Sr*(hr*ur - hl*ul)) / (Sr - Sl)
    else:
        f1m = f2m = 0.0

    denom = hr*(ur - Sr) - hl*(ul - Sl)
    Sm = (Sl*hr*(ur - Sr) - Sr*hl*(ul - Sl)) / denom if denom != 0.0 else 0.0

    if Sl == 0.0 and Sr == 0.0:
        return 0.0, 0.0
    elif Sl >= 0.0:
        return f1l, f2l
    elif Sl <= 0.0 <= Sm:
        return f1m, f2m
    elif Sm <= 0.0 <= Sr:
        return f1m, f2m
    else:
        return f1r, f2r


# ---------------------------------------------------------------------------
# Core stage: compute fluxes then update solution by half a time step
# ---------------------------------------------------------------------------

def _compute_fluxes(h, eta, hu, dx, n):
    """
    MUSCL reconstruction + well-balanced hydrostatic fix + HLLC fluxes.
    Returns (f1, f2, hlx, hrx, zbx).
    """
    N = n + 4
    f1  = np.zeros(N)
    f2  = np.zeros(N)
    hlx = np.zeros(N)
    hrx = np.zeros(N)
    zbx = np.zeros(N)

    for i in range(2, n + 3):          # interfaces i=3..n+3 in MATLAB
        # --- slope limiters ---
        limhl  = float(minmod((h[i-1]-h[i-2])/dx,  (h[i]-h[i-1])/dx))
        limeta = float(minmod((eta[i-1]-eta[i-2])/dx, (eta[i]-eta[i-1])/dx))
        limhul = float(minmod((hu[i-1]-hu[i-2])/dx, (hu[i]-hu[i-1])/dx))

        limhr   = float(minmod((h[i]-h[i-1])/dx,   (h[i+1]-h[i])/dx))
        limetar = float(minmod((eta[i]-eta[i-1])/dx, (eta[i+1]-eta[i])/dx))
        limhur  = float(minmod((hu[i]-hu[i-1])/dx,  (hu[i+1]-hu[i])/dx))

        # --- piecewise-linear reconstruction ---
        hl   = h[i-1]  + 0.5*dx*limhl
        etal = eta[i-1] + 0.5*dx*limeta
        hul_ = hu[i-1] + 0.5*dx*limhul

        hr   = h[i]   - 0.5*dx*limhr
        etar = eta[i]  - 0.5*dx*limetar
        hur_ = hu[i]  - 0.5*dx*limhur

        # --- well-balanced interface bed elevation ---
        zl = etal - hl
        zr = etar - hr
        zbx[i] = max(zl, zr)

        hl = max(0.0, etal - zbx[i])
        hr = max(0.0, etar - zbx[i])
        hlx[i] = hl
        hrx[i] = hr

        f1[i], f2[i] = _hllc(hl, hr, hul_, hur_)

    return f1, f2, hlx, hrx, zbx


def _update(h, hu, f1, f2, hlx, hrx, zbx, dt, dx, n):
    """Apply flux divergence + bed-slope source term, advance by 0.5*dt."""
    for i in range(2, n + 2):          # interior cells i=3..n+2 in MATLAB
        sox = 0.5*(hlx[i+1] + hrx[i]) * (zbx[i+1] - zbx[i]) / dx
        h[i]  -= 0.5*(dt/dx)*(f1[i+1] - f1[i])
        hu[i] -= 0.5*(dt/dx)*(f2[i+1] - f2[i]) + 0.5*dt*GRAV*sox
        if h[i] < 0.0:
            h[i]  = 0.0
            hu[i] = 0.0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def make_grid(xmin, xmax, n):
    """Return (x, dx) for n interior cells plus 2 ghost layers on each side."""
    dx = (xmax - xmin) / n
    x  = np.linspace(xmin - 2*dx, xmax + 2*dx, n + 5)   # n+5 points → n+4 cells
    # cell-centre coordinates
    x = xmin - 1.5*dx + np.arange(n + 4)*dx
    # Note: interior cells at indices 2..n+1 (0-based)
    return x, dx


def compute_dt(h, hu, dx):
    """CFL time step: dt = 0.5*dx / max(|u|+sqrt(g*h)) over all cells."""
    safe_h = np.where(h > DRY_TOL, h, 1.0)
    u = np.where(h > DRY_TOL, hu / safe_h, 0.0)
    c = np.sqrt(GRAV * np.maximum(h, 0.0))
    speed = np.max(np.abs(u) + c)
    if speed == 0.0:
        return 1.0
    return 0.5 * dx / speed


def apply_bc(h, eta, hu, bc_left, bc_right, n):
    """
    Fill ghost cells (indices 0,1 on the left; n+2,n+3 on the right).

    bc_left / bc_right options:
      'transmissive' — copy interior value (open boundary, allows outflow)
      'reflective'   — copy depth/eta, negate momentum (closed wall)
      dict({'h':..., 'eta':..., 'hu':...}) — fixed inflow values
    """
    # left ghosts
    if isinstance(bc_left, dict):
        for arr, key in ((h, 'h'), (eta, 'eta'), (hu, 'hu')):
            arr[0] = arr[1] = bc_left[key]
    elif bc_left == 'reflective':
        h[0] = h[2];   h[1] = h[2]
        eta[0] = eta[2]; eta[1] = eta[2]
        hu[0] = -hu[2]; hu[1] = -hu[2]
    else:  # transmissive
        h[0] = h[2];   h[1] = h[2]
        eta[0] = eta[2]; eta[1] = eta[2]
        hu[0] = hu[2];  hu[1] = hu[2]

    # right ghosts
    if isinstance(bc_right, dict):
        for arr, key in ((h, 'h'), (eta, 'eta'), (hu, 'hu')):
            arr[n+3] = arr[n+2] = bc_right[key]
    elif bc_right == 'reflective':
        h[n+3] = h[n+1];   h[n+2] = h[n+1]
        eta[n+3] = eta[n+1]; eta[n+2] = eta[n+1]
        hu[n+3] = -hu[n+1]; hu[n+2] = -hu[n+1]
    else:  # transmissive
        h[n+3] = h[n+1];   h[n+2] = h[n+1]
        eta[n+3] = eta[n+1]; eta[n+2] = eta[n+1]
        hu[n+3] = hu[n+1];  hu[n+2] = hu[n+1]


def step(h, eta, hu, z, dx, dt, n):
    """
    One full MUSCL-Hancock time step (two half-step stages).
    Mutates h, eta, hu in place.
    """
    for _ in range(2):
        f1, f2, hlx, hrx, zbx = _compute_fluxes(h, eta, hu, dx, n)
        _update(h, hu, f1, f2, hlx, hrx, zbx, dt, dx, n)
        eta[:] = h + z


def run(z, h0, hu0, n, dx, nt, bc_left='transmissive', bc_right='transmissive',
        bc_left_override=None):
    """
    Run the 1D solver for nt time steps.

    Parameters
    ----------
    z   : bed elevation, shape (n+4,)
    h0  : initial water depth, shape (n+4,)
    hu0 : initial x-discharge, shape (n+4,)
    n   : number of interior cells
    dx  : cell width
    nt  : number of time steps
    bc_left, bc_right : 'transmissive' or dict {'h':..., 'eta':..., 'hu':...}

    Returns
    -------
    h, eta, hu : final state arrays
    times      : list of cumulative times
    """
    h   = h0.copy()
    hu  = hu0.copy()
    eta = h + z
    times = [0.0]

    for _ in range(nt):
        apply_bc(h, eta, hu, bc_left, bc_right, n)
        dt = compute_dt(h, hu, dx)
        step(h, eta, hu, z, dx, dt, n)
        times.append(times[-1] + dt)

    return h, eta, hu, times

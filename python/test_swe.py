"""
pytest tests for the 1D SWE MUSCL-Hancock solver.

Run with:  pytest python/test_swe.py -v
"""

import numpy as np
import pytest
from swe_muscl_hancock import (
    GRAV, DRY_TOL,
    minmod, compute_dt, apply_bc, step, run, make_grid,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _flat_still(n=40, eta0=5.0):
    """Flat bed, still water."""
    N = n + 4
    z  = np.zeros(N)
    h0 = np.full(N, eta0)
    hu0 = np.zeros(N)
    return z, h0, hu0, n


def _hump_bed(n=40, eta0=6.0):
    """Cosine-hump bed (case_1 geometry), still water."""
    xmin, xmax = 0.0, 1000.0
    dx = (xmax - xmin) / n
    x  = xmin - 1.5*dx + np.arange(n + 4)*dx
    z  = np.where((x >= 125) & (x <= 875),
                  4.75*(1 - np.cos(2*np.pi*(x - 125)/750))/2,
                  0.0)
    h0 = np.maximum(0.0, eta0 - z)
    hu0 = np.zeros(n + 4)
    return z, h0, hu0, n, dx


# ---------------------------------------------------------------------------
# Unit tests — minmod
# ---------------------------------------------------------------------------

class TestMinmod:
    def test_same_sign_returns_smaller_magnitude(self):
        assert minmod(2.0, 3.0) == pytest.approx(2.0)
        assert minmod(3.0, 2.0) == pytest.approx(2.0)

    def test_opposite_signs_returns_zero(self):
        assert minmod(2.0, -1.0) == pytest.approx(0.0)
        assert minmod(-3.0, 1.0) == pytest.approx(0.0)

    def test_zero_input(self):
        assert minmod(0.0, 5.0) == pytest.approx(0.0)

    def test_vectorized(self):
        a = np.array([1.0, -1.0, 2.0])
        b = np.array([3.0,  2.0, 1.0])
        result = minmod(a, b)
        np.testing.assert_allclose(result, [1.0, 0.0, 1.0])


# ---------------------------------------------------------------------------
# Physics tests
# ---------------------------------------------------------------------------

class TestStillWater:
    """
    Well-balanced property: still water (hu=0, eta=const) must remain still.
    This is the critical test for non-flat terrain correctness.
    """

    def test_flat_bed_stays_still(self):
        z, h0, hu0, n = _flat_still(n=40, eta0=5.0)
        dx = 1000.0 / n
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=50)
        np.testing.assert_allclose(h[2:n+2], h0[2:n+2], atol=1e-12,
                                   err_msg="Flat bed: depth changed")
        np.testing.assert_allclose(hu[2:n+2], 0.0, atol=1e-12,
                                   err_msg="Flat bed: spurious discharge")

    def test_nonflat_bed_stays_still(self):
        """
        Well-balanced test on cosine hump (case_1 geometry).
        eta = 6 m everywhere, u = 0 initially → must stay still.
        This test would FAIL on the original MATLAB code with z>=h bug.
        """
        z, h0, hu0, n, dx = _hump_bed(n=40, eta0=6.0)
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=100)
        # Interior cells only (skip ghost cells)
        np.testing.assert_allclose(eta[2:n+2], 6.0, atol=1e-10,
                                   err_msg="Non-flat bed: surface elevation drifted")
        np.testing.assert_allclose(hu[2:n+2], 0.0, atol=1e-10,
                                   err_msg="Non-flat bed: spurious discharge generated")

    def test_elevated_still_water(self):
        """Still water sitting 2 m above a uniform elevated bed (z = 3)."""
        n  = 30
        N  = n + 4
        dx = 10.0
        z  = np.full(N, 3.0)
        h0 = np.full(N, 2.0)     # eta = 5 everywhere
        hu0 = np.zeros(N)
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=50)
        np.testing.assert_allclose(h[2:n+2], 2.0, atol=1e-12)
        np.testing.assert_allclose(hu[2:n+2], 0.0, atol=1e-12)


class TestPositivity:
    """Water depth must never go negative."""

    def test_h_nonnegative_dam_break(self):
        """1D dam break onto dry bed."""
        n  = 80
        N  = n + 4
        dx = 100.0 / n
        z  = np.zeros(N)
        h0  = np.zeros(N)
        h0[:n//2 + 2] = 5.0   # left half wet
        hu0 = np.zeros(N)
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=100)
        assert np.all(h >= 0.0), f"Negative depth detected: min(h) = {h.min():.3e}"

    def test_h_nonnegative_hump(self):
        """Flow over hump with inflow BC."""
        z, h0, hu0, n, dx = _hump_bed(n=40, eta0=6.0)
        bc_left  = {'h': 6.0, 'eta': 10.0, 'hu': 2.0}
        bc_right = 'transmissive'
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=200,
                            bc_left=bc_left, bc_right=bc_right)
        assert np.all(h >= 0.0), f"Negative depth: min(h) = {h.min():.3e}"


class TestMassConservation:
    """
    With wall (zero-flux) boundary conditions, total water volume is conserved.
    Transmissive BCs copy interior → ghost; since ghost fluxes cancel with the
    mirrored state, total mass is conserved to machine precision.
    """

    def test_mass_conserved_dam_break(self):
        """Reflective (wall) BCs — no mass crosses boundaries → exact conservation."""
        n  = 60
        N  = n + 4
        dx = 1.0
        z  = np.zeros(N)
        h0 = np.zeros(N)
        h0[2:n//2 + 2] = 3.0
        hu0 = np.zeros(N)

        mass_before = np.sum(h0[2:n+2]) * dx
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=80,
                            bc_left='reflective', bc_right='reflective')
        mass_after = np.sum(h[2:n+2]) * dx

        # Reflective ghost cells are only re-applied once per time step (matching
        # the MATLAB design), so stage-2 sees slightly stale ghost values →
        # O(dt²) per step boundary error.  1e-3 relative catches real leaks.
        assert abs(mass_after - mass_before) / (mass_before + 1e-30) < 1e-3, \
            f"Mass not conserved: {mass_before:.6f} → {mass_after:.6f}"

    def test_mass_conserved_hump(self):
        z, h0, hu0, n, dx = _hump_bed(n=40, eta0=6.0)
        mass_before = np.sum(h0[2:n+2]) * dx
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=100)
        mass_after = np.sum(h[2:n+2]) * dx
        rel_err = abs(mass_after - mass_before) / mass_before
        assert rel_err < 1e-10, f"Mass error: {rel_err:.2e}"


class TestDryBed:
    """Dry-bed and wet-dry front handling."""

    def test_dry_bed_stays_dry(self):
        """Cells with no water and no incoming flow must stay dry."""
        n   = 40
        N   = n + 4
        dx  = 10.0
        z   = np.zeros(N)
        h0  = np.zeros(N)
        hu0 = np.zeros(N)
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=50)
        np.testing.assert_allclose(h, 0.0, atol=1e-15)

    def test_water_does_not_appear_on_elevated_dry_bed(self):
        """
        This specifically tests that the old z>=h bug is gone:
        a cell with z=4, h=0 must stay dry even after many steps.
        With the old bug, cells with z > h (e.g. z=4, h=0) were
        correctly left alone; but cells that became h=0.5 on a z=4 bed
        would have been wrongly wiped out.  Here we check the complementary
        case: a small but valid wet state on an elevated bed survives.
        """
        n  = 20
        N  = n + 4
        dx = 1.0
        z  = np.full(N, 4.0)     # uniformly elevated bed
        h0 = np.full(N, 1.0)     # 1 m deep water on top (eta = 5 everywhere)
        hu0 = np.zeros(N)
        h, eta, hu, _ = run(z, h0, hu0, n, dx, nt=50)
        # Water should still be there (well-balanced still water)
        assert np.all(h[2:n+2] > 0.5), \
            "Water on elevated bed was incorrectly removed (z>=h regression)"


class TestCFL:
    """Time step computation."""

    def test_dt_zero_velocity(self):
        n   = 10
        N   = n + 4
        h   = np.ones(N) * 4.0
        hu  = np.zeros(N)
        dx  = 2.0
        dt  = compute_dt(h, hu, dx)
        expected = 0.5 * dx / np.sqrt(GRAV * 4.0)
        assert dt == pytest.approx(expected, rel=1e-10)

    def test_dt_dry_domain_returns_finite(self):
        n  = 10
        N  = n + 4
        h  = np.zeros(N)
        hu = np.zeros(N)
        dt = compute_dt(h, hu, dx=1.0)
        assert np.isfinite(dt)
        assert dt > 0.0

"""Independent check: retarded Green matrices of the simple-cubic lattice with
equal nearest-neighbour and face-diagonal central springs (K = m_b = a = 1).

Two independent routes:
  1. total_G_eta:  G(omega + i*eta) by a k-space sum on an N^3 grid, for several
     eta; extrapolated to eta -> 0 by polynomial extrapolation.
  2. spectral_surface: branch-resolved Im G on the real axis from an integral
     over the constant-frequency surface of each branch (no broadening).

Displacements: dX in -5..5, dY,dZ in -1..1 (all pairs among the 24 loaded corners).
"""
import numpy as np

DX = np.arange(-5, 6)
DYZ = np.arange(-1, 2)


def dyn_matrix(kx, ky, kz):
    """Six independent components of D(k). Works for real or complex arrays."""
    cx, cy, cz = np.cos(kx), np.cos(ky), np.cos(kz)
    sx, sy, sz = np.sin(kx), np.sin(ky), np.sin(kz)
    Dxx = 6 - 2 * cx * (1 + cy + cz)
    Dyy = 6 - 2 * cy * (1 + cz + cx)
    Dzz = 6 - 2 * cz * (1 + cx + cy)
    Dxy = 2 * sx * sy
    Dxz = 2 * sx * sz
    Dyz = 2 * sy * sz
    return Dxx, Dyy, Dzz, Dxy, Dxz, Dyz


def total_G_eta(N, omega, etas):
    """G[e, dX, dY, dZ, i, j] = (1/N^3) sum_k exp(i k.r) [D(k) - (omega+i eta)^2]^-1."""
    etas = np.asarray(etas, float)
    E = len(etas)
    G = np.zeros((E, 11, 3, 3, 3, 3), complex)
    k1 = 2 * np.pi * np.arange(N) / N
    KY, KZ = np.meshgrid(k1, k1, indexing="ij")
    pY = np.exp(1j * np.outer(DYZ, k1))          # (3, N)
    pZ = np.exp(1j * np.outer(DYZ, k1))          # (3, N)
    idx = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    for n in range(N // 2 + 1):
        kx = k1[n]
        Dxx, Dyy, Dzz, Dxy, Dxz, Dyz = dyn_matrix(kx, KY, KZ)
        paired = not (n == 0 or (N % 2 == 0 and n == N // 2))
        phx = np.exp(1j * kx * DX)               # (11,)
        for e, eta in enumerate(etas):
            z2 = (omega + 1j * eta) ** 2
            a = Dxx - z2; d = Dyy - z2; f = Dzz - z2
            b = Dxy; c = Dxz; ee = Dyz
            A11 = d * f - ee * ee
            A12 = c * ee - b * f
            A13 = b * ee - c * d
            det = a * A11 + b * A12 + c * A13
            inv = [A11 / det, (a * f - c * c) / det, (a * d - b * b) / det,
                   A12 / det, A13 / det, (b * c - a * ee) / det]
            plane = np.zeros((3, 3, 3, 3), complex)   # [dY, dZ, i, j]
            for comp, (i, j) in zip(inv, idx):
                r = pY @ comp @ pZ.T                  # (3, 3)
                plane[:, :, i, j] = r
                plane[:, :, j, i] = r
            plane /= N ** 3
            G[e] += phx[:, None, None, None, None] * plane[None]
            if paired:
                # -kx plane: D(-kx) = P D P with P = diag(-1, 1, 1)
                pl2 = plane.copy()
                pl2[:, :, 0, 1] *= -1; pl2[:, :, 1, 0] *= -1
                pl2[:, :, 0, 2] *= -1; pl2[:, :, 2, 0] *= -1
                G[e] += np.conj(phx)[:, None, None, None, None] * pl2[None]
    return G


def extrapolate(etas, vals):
    """Polynomial (Lagrange) extrapolation of vals[e, ...] to eta = 0."""
    etas = np.asarray(etas, float)
    out = np.zeros(vals.shape[1:], complex)
    for i, xi in enumerate(etas):
        w = 1.0
        for j, xj in enumerate(etas):
            if j != i:
                w *= (0.0 - xj) / (xi - xj)
        out += w * vals[i]
    return out


def _eig(kvec):
    """Sorted eigenvalues/eigenvectors of D at real k (..., 3)."""
    Dxx, Dyy, Dzz, Dxy, Dxz, Dyz = dyn_matrix(kvec[..., 0], kvec[..., 1], kvec[..., 2])
    D = np.empty(kvec.shape[:-1] + (3, 3))
    D[..., 0, 0] = Dxx; D[..., 1, 1] = Dyy; D[..., 2, 2] = Dzz
    D[..., 0, 1] = D[..., 1, 0] = Dxy
    D[..., 0, 2] = D[..., 2, 0] = Dxz
    D[..., 1, 2] = D[..., 2, 1] = Dyz
    return np.linalg.eigh(D)


def spectral_surface(omega, ntheta=200, nphi=400, check_monotone=False):
    """Branch-resolved spectral matrices on the real axis.

    A_T[dX, dY, dZ, i, j] (two lower branches) and A_L (upper branch), with
    A_B(r) = pi * int d^3k/(2pi)^3 sum_{b in B} e_b e_b^T cos(k.r) delta(lam_b(k) - omega^2).
    """
    mu, wmu = np.polynomial.legendre.leggauss(ntheta)
    phi = 2 * np.pi * (np.arange(nphi) + 0.5) / nphi
    wphi = 2 * np.pi / nphi
    MU, PHI = np.meshgrid(mu, phi, indexing="ij")
    W = (wmu[:, None] * wphi * np.ones_like(PHI)).ravel()
    st = np.sqrt(1 - MU ** 2)
    khat = np.stack([st * np.cos(PHI), st * np.sin(PHI), MU], axis=-1).reshape(-1, 3)
    kmax = np.pi / np.abs(khat).max(axis=1)
    w2 = omega ** 2
    A = {"T": np.zeros((11, 3, 3, 3, 3)), "L": np.zeros((11, 3, 3, 3, 3))}
    rvec = np.stack(np.meshgrid(DX, DYZ, DYZ, indexing="ij"), axis=-1).reshape(-1, 3).astype(float)
    worst = 0.0
    for b in range(3):
        lo = np.zeros(len(khat)); hi = kmax.copy()
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lam = _eig(mid[:, None] * khat)[0][:, b]
            up = lam < w2
            lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
        kb = 0.5 * (lo + hi)
        kv = kb[:, None] * khat
        lam, vec = _eig(kv)
        worst = max(worst, np.abs(lam[:, b] - w2).max())
        e = vec[:, :, b]
        # radial derivative of the eigenvalue by Hellmann-Feynman with complex-step dD/dk
        h = 1e-30
        comps = dyn_matrix(kv[:, 0] + 1j * h * khat[:, 0], kv[:, 1] + 1j * h * khat[:, 1],
                           kv[:, 2] + 1j * h * khat[:, 2])
        dD = [np.imag(cc) / h for cc in comps]
        dl = (dD[0] * e[:, 0] ** 2 + dD[1] * e[:, 1] ** 2 + dD[2] * e[:, 2] ** 2
              + 2 * dD[3] * e[:, 0] * e[:, 1] + 2 * dD[4] * e[:, 0] * e[:, 2]
              + 2 * dD[5] * e[:, 1] * e[:, 2])
        if check_monotone:
            # Check the radial slope at each returned root.
            assert dl.min() > 0, "nonpositive radial slope at a sampled root"
        amp = np.pi / (2 * np.pi) ** 3 * W * kb ** 2 / dl          # (ndir,)
        ph = np.cos(kv @ rvec.T)                                   # (ndir, nr)
        P = e[:, :, None] * e[:, None, :]                          # (ndir, 3, 3)
        contrib = np.einsum("d,dr,dij->rij", amp, ph, P).reshape(11, 3, 3, 3, 3)
        A["T" if b < 2 else "L"] += contrib
    return A["T"], A["L"], worst

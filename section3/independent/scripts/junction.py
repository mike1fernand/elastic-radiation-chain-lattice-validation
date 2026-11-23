"""Scattering of a chain compression wave at the spring contact with the cubic lattice.

Units K = m_b = a = 1.  Chain: m_q, K_q.  Segments i = 1..n (stiffness w[i-1] per spring)
join chain sites i-1 and i; each endpoint is tied to the four corners of the shared face.
Incident wave exp(i(k j - omega t)) from j -> -infinity.
"""
import numpy as np

S3 = np.sqrt(3.0)


def assemble(w):
    n = len(w)
    sites = []          # lattice integer coords of loaded corners
    U = np.zeros((12 * n, 12 * n))
    V = np.zeros((12 * n, n + 1))
    Kself = np.zeros(n + 1)
    c = 0
    for i in range(1, n + 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                sites.append((i - 1, (sy - 1) // 2, (sz - 1) // 2))
                nL = np.array([1.0, sy, sz]) / S3
                nR = np.array([-1.0, sy, sz]) / S3
                sl = slice(3 * c, 3 * c + 3)
                U[sl, sl] = w[i - 1] * (np.outer(nL, nL) + np.outer(nR, nR))
                V[sl, i - 1] += w[i - 1] * nL * nL[0]
                V[sl, i] += w[i - 1] * nR * nR[0]
                Kself[i - 1] += w[i - 1] * nL[0] ** 2
                Kself[i] += w[i - 1] * nR[0] ** 2
                c += 1
    return np.array(sites), U, V, Kself


def site_matrix(sites, G):
    """72x72 matrix from G[dX+5, dY+1, dZ+1, i, j]."""
    ns = len(sites)
    M = np.zeros((3 * ns, 3 * ns), dtype=G.dtype)
    for a in range(ns):
        for b in range(ns):
            d = sites[a] - sites[b]
            M[3 * a:3 * a + 3, 3 * b:3 * b + 3] = G[d[0] + 5, d[1] + 1, d[2] + 1]
    return M


def solve(w, omega, ReG, AT, AL, Kq=1.69, mq=1.0):
    """Return dict with R, T_line, T_qT, T_qL (fractions of incident power)."""
    w = np.asarray(w, float)
    n = len(w)
    sites, U, V, Kself = assemble(w)
    aT = site_matrix(sites, AT)
    aL = site_matrix(sites, AL)
    g0 = site_matrix(sites, ReG) + 1j * (aT + aL)
    k = np.arccos(1 - mq * omega ** 2 / (2 * Kq))
    S = np.linalg.solve(np.eye(12 * n) + g0 @ U, g0 @ V)      # u = S q
    Sig = V.T @ S
    A = np.zeros((n + 1, n + 1), complex)
    for j in range(n + 1):
        A[j, j] = 2 * Kq - mq * omega ** 2 + Kself[j]
        if j > 0:
            A[j, j - 1] = -Kq
        if j < n:
            A[j, j + 1] = -Kq
    A -= Sig
    A[0, 0] -= Kq * np.exp(1j * k)
    A[n, n] -= Kq * np.exp(1j * k)
    b = np.zeros(n + 1, complex)
    b[0] = -2j * Kq * np.sin(k)            # unit incident amplitude
    q = np.linalg.solve(A, b)
    u = S @ q
    f = V @ q - U @ u                       # force on the bulk corners
    Pinc = 0.5 * omega * Kq * np.sin(k)
    R = abs(q[0] - 1.0) ** 2
    Tl = abs(q[n]) ** 2
    PT = 0.5 * omega * np.real(np.conj(f) @ aT @ f)
    PL = 0.5 * omega * np.real(np.conj(f) @ aL @ f)
    Fx = f.reshape(-1, 3)[:, 0].sum()
    return dict(R=R, Tline=Tl, TqT=PT / Pinc, TqL=PL / Pinc,
                balance=R + Tl + (PT + PL) / Pinc, k=k, Fx=Fx, q=q)

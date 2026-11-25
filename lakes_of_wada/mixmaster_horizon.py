"""Particle-horizon scales along a Mixmaster (Bianchi IX) orbit: a Theta-like quantity per axis.

Setting (QINERTIA.md section 7).  Same equations and slice as ``mixmaster``: q_i = a_i^2, x_i = ln q_i, BKL time tau
(dt = abc dtau), integrated from the maximum-expansion slice towards the singularity (S = sum x' < 0).  With c = 1 the
comoving particle horizon along axis i at cosmic time t is  H_i = integral_0^t dt'/a_i , and the proper one
D_i = a_i H_i.  In BKL time  dt/a_i = sqrt(q_j q_k) dtau  and  t = integral abc dtau, so both are integrals along the orbit
taken from the current point to the singularity (tau -> infinity).  We take the Theta-like scale of McCulloch and Gine as
Theta_i = 2 D_i and  a0_i = 2 c^2 / Theta_i = 1 / D_i  (the paper's Theta is the co-moving diameter of the observable universe;
here the horizon is axis dependent and measured along the orbit).  This is a definition made here, not the paper's.

Computation.  The integrals from the current point to the end of the run are accumulated in reverse (no cancellation), and the
rest of the way to the singularity is added with the Kasner tail of the final epoch (x_i' = S p_i):
  integral_tau^inf sqrt(q_j q_k) = sqrt(q_j q_k) / (|S| (1 - p_i)/2),     t_tail = 2 abc / |S|.
In a Kasner epoch that continues to the singularity, D_i / t = 1/(1 - p_i) exactly (a_i = A t^p_i); along a Mixmaster orbit the
later bounces change it, which is the quantity of interest.
"""
import numpy as np

from .mixmaster import force, potential, turnaround_state


def horizon_orbits(bp, bm, h0=0.05, tau_budget=60.0, tau_cap=120.0, kappa=1e-3, thin=10, nmax=30000, x0=None, v0=None):
    """Integrate the maximum-expansion states at (bp, bm) (or explicit ``x0``, ``v0`` of shape (m, 3)) and return a list of
    dicts, one per orbit, with arrays over the thinned samples: ``tau``, ``x`` (m, 3), ``p`` (= v / S), ``kap`` (= |V| / sum v^2),
    ``t`` (cosmic time to the singularity), ``H`` (comoving horizons, 3 columns), ``D`` (proper horizons a_i H_i), ``R`` (= D / t),
    plus ``kasner_end`` (True if the run ended in a Kasner epoch, which the tail estimate needs).  Orbits stop once tau >= tau_budget and
    they are in a Kasner epoch (kap < kappa, S < 0), or at tau_cap."""
    if x0 is None:
        x0, v0, V0 = turnaround_state(np.ravel(bp), np.ravel(bm))
        ok = V0 < 0
        x0, v0 = x0[ok], v0[ok]
    m = len(x0)
    x, v = x0.copy(), v0.copy()
    tau = np.zeros(m)
    maxrec = nmax // thin + 2
    rec_tau = np.full((m, maxrec), np.nan)
    rec_x = np.full((m, maxrec, 3), np.nan)
    rec_p = np.full((m, maxrec, 3), np.nan)
    rec_kap = np.full((m, maxrec), np.nan)
    rec_dG = np.zeros((m, maxrec, 3))
    rec_dT = np.zeros((m, maxrec))
    nrec = np.zeros(m, dtype=int)
    accG = np.zeros((m, 3))
    accT = np.zeros(m)
    end_S = np.zeros(m)
    end_kas = np.zeros(m, dtype=bool)
    alive = np.arange(m)

    def integrands(xa):
        e = np.exp(0.5 * xa)                       # a_i
        g = np.stack([e[:, 1] * e[:, 2], e[:, 2] * e[:, 0], e[:, 0] * e[:, 1]], axis=1)   # sqrt(q_j q_k) = b c, c a, a b
        return g, e[:, 0] * e[:, 1] * e[:, 2]

    for it in range(1, nmax + 1):
        if alive.size == 0:
            break
        xa, va = x[alive], v[alive]
        h = (h0 / (1.0 + np.max(np.abs(va), axis=1)))[:, None]
        g0, a0_ = integrands(xa)
        k1x, k1v = va, force(xa)
        k2x, k2v = va + 0.5 * h * k1v, force(xa + 0.5 * h * k1x)
        k3x, k3v = va + 0.5 * h * k2v, force(xa + 0.5 * h * k2x)
        k4x, k4v = va + h * k3v, force(xa + h * k3x)
        xa = xa + h / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
        va = va + h / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        g1, a1_ = integrands(xa)
        accG[alive] += 0.5 * h * (g0 + g1)
        accT[alive] += 0.5 * h[:, 0] * (a0_ + a1_)
        x[alive], v[alive] = xa, va
        tau[alive] += h[:, 0]
        S = va.sum(axis=1)
        kap = np.abs(potential(xa)) / np.maximum(np.sum(va ** 2, axis=1), 1e-300)
        kas = (kap < kappa) & (S < 0)
        done = ((tau[alive] >= tau_budget) & kas) | (tau[alive] >= tau_cap)
        sel = np.arange(len(alive)) if it % thin == 0 else np.where(done)[0]
        if sel.size:
            idx = alive[sel]
            k = nrec[idx]
            rec_tau[idx, k] = tau[idx]
            rec_x[idx, k] = xa[sel]
            rec_p[idx, k] = va[sel] / np.where(S[sel] == 0, 1.0, S[sel])[:, None]
            rec_kap[idx, k] = kap[sel]
            rec_dG[idx, k] = accG[idx]
            rec_dT[idx, k] = accT[idx]
            accG[idx] = 0.0
            accT[idx] = 0.0
            nrec[idx] += 1
        end_S[alive[done]] = S[done]
        end_kas[alive[done]] = kas[done]
        alive = alive[~done]
    out = []
    for i in range(m):
        n = nrec[i]
        if n < 2 or end_S[i] == 0:
            continue
        x_ = rec_x[i, :n]
        p_ = rec_p[i, :n]
        e = np.exp(0.5 * x_)
        g_end = np.array([e[-1, 1] * e[-1, 2], e[-1, 2] * e[-1, 0], e[-1, 0] * e[-1, 1]])
        abc_end = e[-1].prod()
        Sabs = abs(end_S[i])
        tail_H = g_end / (0.5 * Sabs * np.maximum(1.0 - p_[-1], 1e-6))
        tail_t = abc_end / (0.5 * Sabs)
        # reverse cumulative sums of the increments recorded at each sample, plus the tail
        H = np.cumsum(rec_dG[i, :n][::-1], axis=0)[::-1] + tail_H
        t = np.cumsum(rec_dT[i, :n][::-1])[::-1] + tail_t
        # the increment stored at sample k covers the steps *up to* k, so the horizon at sample k excludes it
        H = H - rec_dG[i, :n]
        t = t - rec_dT[i, :n]
        D = e * H
        out.append(dict(tau=rec_tau[i, :n], x=x_, p=p_, kap=rec_kap[i, :n], t=t, H=H, D=D, R=D / t[:, None],
                        kasner_end=bool(end_kas[i]), start=(x0[i], v0[i])))
    return out


def epochs(orbit, kappa=1e-3):
    """Kasner epochs of one orbit from the samples (enter at kap < kappa, leave above 10 kappa, as in ``mixmaster``).  Returns a list of
    dicts with the sample of minimum kap in each epoch: ``i`` (sample index), ``tau``, ``p``, ``R`` (D_i / t), ``kasner`` (= 1/(1 - p))."""
    kap, S_neg = orbit["kap"], True
    inK, cur, res = False, [], []
    for i, k in enumerate(kap):
        now = (k < 10 * kappa) if inK else (k < kappa)
        if now:
            cur.append(i)
        elif inK and cur:
            res.append(cur)
            cur = []
        inK = now
    if inK and cur:
        res.append(cur)
    out = []
    for c in res:
        j = c[int(np.argmin(kap[c]))]
        p = orbit["p"][j]
        out.append(dict(i=j, tau=float(orbit["tau"][j]), p=p, R=orbit["R"][j], kasner=1.0 / (1.0 - p)))
    return out

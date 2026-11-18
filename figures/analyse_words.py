#!/usr/bin/env python3
"""Reduced-word coding of the Mixmaster exit basins (persistent exit rule).

    python3 figures/analyse_words.py [n]        n x n grid on [-1,1]^2 of (beta_+, beta_-), default 160 (about 3 minutes)

The wall word of a trajectory is the sequence of bounces (wall i = x_i' going from + to -).  It is a reduced word in the
(inf,inf,inf) group of HOPF_BKL.md (consecutive letters differ).  Measured:
  1  how often the exit axis equals the last letter, and the exceptions;
  2  H(label | last k letters) and H(label | first k letters), plug-in entropy in bits, groups of at least 20 states;
  3  mixing of the cylinders: the fraction of cylinders of the first k letters (at least 30 states) in which all three labels
     occur with at least 10 % each.
Writes figures/mixmaster_words.png and figures/mixmaster_words.csv.
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from lakes_of_wada import mixmaster as M            # noqa: E402

n = int(sys.argv[1]) if len(sys.argv) > 1 else 160
cache = HERE / ("mixmaster_words_%d.npz" % n)
a = (np.arange(n) + 0.5) / n * 2 - 1
BP, BM = np.meshgrid(a, a, indexing="ij")
if cache.exists():
    z = np.load(cache); lab, W = z["lab"], z["W"]
else:
    lab, W = M.exit_words(BP, BM)
    np.savez_compressed(cache, lab=lab, W=W)
Lw = (W > 0).sum(axis=-1)
dec = lab > 0
print("grid %d x %d: valid %d, decided %d (%.1f %%), median word length of decided %d" %
      (n, n, (lab >= 0).sum(), dec.sum(), 100 * dec.sum() / (lab >= 0).sum(), np.median(Lw[dec])))
words = [tuple(int(c) for c in W[i, j, :Lw[i, j]]) for i, j in zip(*np.where(dec))]
labs = [int(lab[i, j]) for i, j in zip(*np.where(dec))]
rows = []

# 1. last letter
last = Counter((w[-1] if w else 0, y) for w, y in zip(words, labs))
tot = sum(v for (l, y), v in last.items() if l > 0)
hit = sum(v for (l, y), v in last.items() if l == y)
print("1. exit axis == last letter: %d of %d (%.1f %%)" % (hit, tot, 100 * hit / tot))
rows.append(["last_letter_equals_exit", "", hit / tot, tot])
exc = Counter()
for w, y in zip(words, labs):
    if len(w) >= 2 and w[-1] != y:
        exc[(w[-2], w[-1], y)] += 1
print("   exceptions by (second to last, last, exit):", dict(exc.most_common(6)))


def cond_entropy(keys, ys, minsize=20):
    g = defaultdict(list)
    for k, y in zip(keys, ys):
        g[k].append(y)
    H = N = 0.0
    for v in g.values():
        if len(v) < minsize:
            continue
        c = np.array(list(Counter(v).values()), float) / len(v)
        H += -len(v) * np.sum(c * np.log2(c)); N += len(v)
    return (H / N if N else np.nan), int(N), len(g)


c = np.array(list(Counter(labs).values()), float) / len(labs)
H0 = float(-np.sum(c * np.log2(c)))
print("2. H(label) = %.3f bits" % H0)
print("   k   H(label | last k)  covered   H(label | first k)  covered")
for k in range(1, 9):
    tl = [(w[-k:], y) for w, y in zip(words, labs) if len(w) >= k]
    hd = [(w[:k], y) for w, y in zip(words, labs) if len(w) >= k]
    Ht, Nt, _ = cond_entropy(*zip(*tl)) if tl else (np.nan, 0, 0)
    Hh, Nh, _ = cond_entropy(*zip(*hd)) if hd else (np.nan, 0, 0)
    print("   %d   %.3f             %5d     %.3f              %5d" % (k, Ht, Nt, Hh, Nh))
    rows += [["H_given_last", k, Ht, Nt], ["H_given_first", k, Hh, Nh]]

# 3. mixing of cylinders
print("3. cylinders of the first k letters with >= 30 states: how many have all three labels with >= 10 %% each")
for k in range(1, 8):
    g = defaultdict(list)
    for w, y in zip(words, labs):
        if len(w) >= k:
            g[w[:k]].append(y)
    big = [v for v in g.values() if len(v) >= 30]
    mixed = sum(1 for v in big if all(np.mean(np.array(v) == i) >= 0.10 for i in (1, 2, 3)))
    cov = sum(len(v) for v in big)
    print("   k=%d: %d cylinders in total, %d with >= 30 states, %d mixed (%.0f %%), covering %d states" %
          (k, len(g), len(big), mixed, 100 * mixed / max(len(big), 1), cov))
    rows.append(["mixed_cylinders", k, mixed / max(len(big), 1), len(big)])

# figure
first3 = np.full((n, n), -1)
codes = {}
for i, j in zip(*np.where(dec)):
    w = tuple(int(c) for c in W[i, j, :3])
    if len(w) == 3 and Lw[i, j] >= 3:
        first3[i, j] = codes.setdefault(w, len(codes))
fig, axs = plt.subplots(1, 3, figsize=(16, 5.2))
ext = (-1, 1, -1, 1)
axs[0].imshow(np.where(lab > 0, lab, 0).T, origin="lower", extent=ext, cmap=ListedColormap(["#ddd", "#d1495b", "#2e86ab", "#edae49"]), vmin=0, vmax=3)
axs[0].set_title("exit axis (grey: no persistent exit)")
im = axs[1].imshow(np.where(first3 >= 0, first3, np.nan).T, origin="lower", extent=ext, cmap="tab20")
axs[1].set_title("cylinder of the first 3 bounces (%d cylinders)" % len(codes))
lw = np.where(dec, Lw, np.nan)
im = axs[2].imshow(lw.T, origin="lower", extent=ext, cmap="viridis", vmax=np.nanpercentile(lw, 98))
plt.colorbar(im, ax=axs[2], label="number of bounces before the exit"); axs[2].set_title("word length")
for ax in axs:
    ax.set_xlabel("beta+"); ax.set_ylabel("beta-")
fig.tight_layout(); fig.savefig(HERE / "mixmaster_words.png", dpi=110)
with open(HERE / "mixmaster_words.csv", "w", newline="") as fh:
    w_ = csv.writer(fh); w_.writerow(["kind", "k", "value", "n"]); w_.writerows(rows)

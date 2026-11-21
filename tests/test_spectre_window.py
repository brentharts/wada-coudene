import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from lakes_of_wada import spectre_window as SW


def test_piece_counts_follow_the_shape_matrix():
    c = [SW.count("O", 0, d) for d in range(10)]
    assert c[:4] == [1, 13, 103, 811]
    assert SW.count("G", 1, 3) == 118 and SW.count("O", 1, 5) == SW.count("O", 0, 5)
    assert abs(c[9] / c[8] - (4 + np.sqrt(15))) < 1e-3          # Perron root of the shape matrix, lam^2 = 4 + sqrt 15
    assert abs(SW.count("O", 0, 8) / SW.count("G", 0, 8) - 6.873) < 0.01


def test_points_stay_inside_the_radius_bound_and_depths_agree():
    d = np.load(Path(__file__).resolve().parent.parent / "figures" / "spectre_window_maps.npz")
    rad, lam = float(d["rad"]), float(d["lam"])
    p5, _ = SW.points("O", 0, 5); p6, _ = SW.points("O", 0, 6)
    assert np.abs(p5).max() < rad and np.abs(p6).max() < rad
    from scipy.spatial import cKDTree
    dist, _ = cKDTree(p6).query(p5)                              # depth-5 points lie within rad * lam^-5 of depth-6 points
    assert dist.max() < rad * lam ** -5


def test_label_image_labels_children_and_marks_outside():
    pts, kid = SW.points("O", 0, 5)
    im, box, px = SW.label_image(pts, kid, 120)
    assert set(np.unique(im)) <= set(range(-1, 13)) and (im == -1).any() and len(set(np.unique(im)) - {-1}) == 13
    assert im[0, 0] == -1 and im[-1, -1] == -1


def test_neighbour_window_partition_tables_are_not_wada():
    """Committed output of figures/spectre_window_covering.py: three or more pieces meet at a minority of boundary pixels and the fraction grows with r."""
    import csv
    from pathlib import Path
    fig = Path(__file__).resolve().parent.parent / "figures"
    for c in (0, 2, 5):
        rows = [r for r in csv.reader(open(fig / f"spectre_window_partition_N3_c{c}.csv")) if r and r[0] not in ("r", "probe_pieces")]
        p3 = [float(r[6]) for r in rows]
        assert len(p3) == 4 and all(a < b for a, b in zip(p3, p3[1:])) and max(p3) < 0.25


def test_phason_flip_exponent_close_to_two_minus_seam_dimension():
    """Gaussian displacement of the internal point: the exit fraction of Z_0(O) scales like delta^(2 - D), D = 1.1110 (to about 0.05 at this resolution)."""
    from scipy import ndimage as ndi
    pts, kid = SW.points("O", 0, 6)
    lo, hi = pts.min(0), pts.max(0); c = 0.5 * (lo + hi); h = 0.5 * float((hi - lo).max()) + 0.6
    lab, box, px = SW.label_image(pts, kid, 900, box=(c[0] - h, c[0] + h, c[1] - h, c[1] + h), reach=1.5)
    inw = lab >= 0
    ds = np.array([0.02, 0.03, 0.045, 0.07])
    f = [1 - ndi.gaussian_filter(inw.astype(float), d / px, mode="constant")[inw].mean() for d in ds]
    alpha = np.polyfit(np.log(ds), np.log(f), 1)[0]
    D = np.log(5 + 2 * np.sqrt(6)) / np.log(4 + np.sqrt(15))
    assert abs(alpha - (2 - D)) < 0.06

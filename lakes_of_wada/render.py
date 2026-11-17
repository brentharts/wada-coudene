"""Image output."""
import numpy as np
from matplotlib import pyplot as plt


def save_image(pic, path, cmap="magma", background="white", blank_labels=None):
    """Save ``pic.color`` as an image with a white (default) background.

    Pixels whose basin label is in ``blank_labels`` are drawn in ``background``.
    Default: label 0 (undecided) always, plus label 1 (the origin basin) for
    plane pictures, where it is the unbounded outer region.
    """
    if blank_labels is None:
        blank_labels = (0, 1) if pic.params.get("kind") == "plane" else (0,)
    img = np.ma.masked_where(np.isin(pic.label, blank_labels) | (pic.label < 0), pic.color)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(background)
    plt.imsave(path, img, cmap=cm)


def save_scalar(field, path, cmap="viridis", vmin=None, vmax=None, background="white"):
    """Save a float field (NaN -> background), e.g. a maverick score or Wada map."""
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(background)
    plt.imsave(path, np.ma.masked_invalid(field), cmap=cm, vmin=vmin, vmax=vmax)

"""lakes_of_wada: dynamical Lakes of Wada after Coudene, in pure numpy/scipy."""
from . import maps, classify, views, elliptic, entire, wadatest, maverick, construction, scattering, sweep, mixmaster, yangmills, quantised_inertia, mixmaster_horizon
from .maps import Anosov, Perturbation, torus_map, plane_map
from .classify import orbit_to_targets, Result
from .render import save_image, save_scalar
from .pipeline import Picture, torus_picture, plane_picture, entire_picture

__all__ = ["maps", "classify", "views", "elliptic", "Anosov", "Perturbation",
           "torus_map", "plane_map", "orbit_to_targets", "Result", "Picture",
           "torus_picture", "plane_picture", "save_image", "save_scalar", "entire", "entire_picture", "wadatest", "maverick", "construction", "scattering", "sweep", "mixmaster", "yangmills", "quantised_inertia", "mixmaster_horizon"]

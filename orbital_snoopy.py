import orbital as _orb
from scipy.constants import kilo as _kilo

sol = _orb.bodies.Body(
		mass=_orb.constants.sun_mass,
		mu=_orb.constants.solar_mass_parameter,
		mean_radius=_orb.constants.sun_radius_equatorial,
		equatorial_radius=_orb.constants.sun_radius_equatorial,
		polar_radius=_orb.constants.sun_radius_equatorial,
		apoapsis_names='aphelion',
		periapsis_names='perihelion',
		plot_color='#d0d010')


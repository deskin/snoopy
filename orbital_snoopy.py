import orbital as _orb
from scipy.constants import kilo as _kilo
import astropy.time as _astrotime

sol = _orb.bodies.Body(
		mass=_orb.constants.sun_mass,
		mu=_orb.constants.solar_mass_parameter,
		mean_radius=_orb.constants.sun_radius_equatorial,
		equatorial_radius=_orb.constants.sun_radius_equatorial,
		polar_radius=_orb.constants.sun_radius_equatorial,
		apoapsis_names='aphelion',
		periapsis_names='perihelion',
		plot_color='#d0d010')

def earth_orbit_2017():
		r = _orb.Position(x=-5.028580283544673e7, y=-1.422272856162977e8, z=-1.550277271285653e4)
		r = r * _kilo
		v = _orb.Velocity(x=2.758282089460933e1, y=-1.005988238736911e1, z=4.542609417019783e-4)
		v = v * _kilo
		t = _astrotime.Time(2457905.5, format='jd', scale='tdb')

		o = _orb.elements.KeplerianElements.from_state_vector(r=r, v=v, body=sol, ref_epoch=t)

		return o


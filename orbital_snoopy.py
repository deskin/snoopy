def OrbitalPy():
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

import pykep as _pk

def orbit_2018_av2_pk():
    e = _pk.epoch(2458137.5, 'jd')
    o = _pk.planet.keplerian(
            e,
            (1.029147132 * _pk.AU, 0.02933445, 0.122007 * _pk.DEG2RAD, 347.610893 * _pk.DEG2RAD, 110.815416 * _pk.DEG2RAD, 19.215145973 * _pk.DEG2RAD),
            _pk.MU_SUN,
            0.1,
            10,
            20,
            '2018 AV2')

    return o

def norm(v1, v2):
    from math import sqrt
    return sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)

def porkchop(o1, o2, t_begin_range, t_end_range):
    out = []
    i = 0.0
    while i < t_begin_range[2]:
        t1 = _pk.epoch(t_begin_range[0].mjd2000 + i)
        outrow = []
        out.append((t1, outrow))
        j = 0.0
        while j < t_end_range[2]:
            t2 = _pk.epoch(t_end_range[0].mjd2000 + j)
            dt = (t2.jd - t1.jd) * _pk.DAY2SEC
            r1, v1 = o1.eph(t1)
            r2, v2 = o2.eph(t2)
            l = _pk.lambert_problem(r1=r1, r2=r2, tof=dt, mu=_pk.MU_SUN, max_revs = 1)
            vb = l.get_v1()[0]
            ve = l.get_v2()[0]
            dv1 = norm(v1, vb)
            dv2 = norm(v2, ve)
            outrow.append((t2, dt, dv1 + dv2, dv1, dv2))

            j += t_end_range[1]

        i += t_begin_range[1]

    return out

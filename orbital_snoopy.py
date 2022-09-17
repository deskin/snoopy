# Copyright 2022 Deskin Miller
# Licensed under the MIT License. See COPYING for details

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
        # t1 transfer departure date
        t1 = _pk.epoch(t_begin_range[0].mjd2000 + i)
        outrow = []
        out.append((t1, outrow))
        j = 0.0
        while j < t_end_range[2]:
            # t2 transfer arrival date
            t2 = _pk.epoch(t_end_range[0].mjd2000 + j)
            dt = (t2.jd - t1.jd) * _pk.DAY2SEC
            # Calculate ephemeris for originating body at departure and
            # destination body at arrival
            r1, v1 = o1.eph(t1)
            r2, v2 = o2.eph(t2)

            l = _pk.lambert_problem(r1=r1, r2=r2, tof=dt, mu=_pk.MU_SUN, max_revs = 1)
            # vb beginning velocity on transfer orbit
            # ve end velocity on transfer orbit
            vb = l.get_v1()[0]
            ve = l.get_v2()[0]

            # dv1 delta-V maneuver to transfer ignoring originating body gravity
            # dv2 delta-V maneuver to match destination orbit
            dv1 = norm(v1, vb)
            dv2 = norm(v2, ve)
            outrow.append((t2, dt, dv1 + dv2, dv1, dv2))

            j += t_end_range[1]

        i += t_begin_range[1]

    return out

# (2041-Mar-28 00:00:00, (2041-Nov-30 00:00:00, 21340800.0, 480.0902720093982, 192.937420833875, 287.15285117552315))
# (2041-Apr-09 00:00:00, (2041-Dec-07 00:00:00, 20908800.0, 465.410659561673, 174.42675438345546, 290.9839051782176))

def porkchop_2041_transfer():
    return (
        _pk.planet.jpl_lp('earth'),
        orbit_2018_av2_pk(),
        (_pk.epoch_from_string('2041-03-01 00:00:00'), 1.0, 62.1),
        (_pk.epoch_from_string('2041-11-01 00:00:00'), 1.0, 62.1))

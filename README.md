# Approximate Analysis of 2018 AV2 Orbit, Next Close Approach, Hohmann Transfer Window for TOMTOM/ Project Snoopy
Deskin Miller (deskinm at umich dot edu)  
September 2022

---
## Summary
By analyzing the orbital elements of 2018 AV2 we approximate the next closest Earth approach to be 2041 September 1 (+/- 1 day). We find an approximate Hohmann transfer window open between "2041 March 1 -- 2041 May 1" and arrival between "2041 November 1 -- 2042 January 1" with Earth-departure "$C_3 < 1 km^2/s^2$" and rendezvous "delta-V $< 500 m/s$". This can be used for evaluating project feasibility, further refinement of orbital dynamics, and initial constraints on spacecraft selection.

## Structure of this Document
Background and Orbital Parameters; 2041 Close Approach; Hohmann Transfer Window; Limitations; Potential future work; References

## Background and 2018 AV2 Orbital Elements
The object designated 2018 AV2 was detected in late 2017. The orbital elements used are obtained from [[1]](#1-pseudo-mpec-for-za9872d--2018-av2):

Epoch 2018 Jan 19 00:00:00/ 2458137.5 JD  
Semimajor axis $a = 1.029147132 +/- 2.1e-6 AU$  
Eccentricity $e = 0.02933445 +/- 8.09e-7$ (scalar)  
Inclination $i = 0.122007 +/- 0.000010 deg$  
Argument of periapsis $\omega = 110.815416 +/- 0.0013 deg$  
Longitude of the ascending node $\Omega = 347.610893 +/- 0.0009 deg$  
Mean anomaly $M = 19.215145973 +/- 0.00050 deg$  

Note: we do not use the uncertainty estimates in any calculations here but they are listed for reference.

## 2041 Close Approach
Using the Keplerian orbits for Earth and 2018 AV2 in [pykep [2]](#2-pykep) we find the closest approach:
``` py
import pykep as pk
import orbital_snoopy as sn

o = sn.orbit_2018_av2_pk()
earth = pk.planet.jpl_lp('earth')

e_start = pk.epoch_from_string('2022-09-01 00:00:00')
e_end = pk.epoch_from_string('2050-01-01 00:00:00')

def earth_2018_av2_close_approach(t1, t2, step):
    e = t1
    mn = (None, 1.0e20)
    while e.jd < t2.jd:
        r1, v1 = o.eph(e)
        r2, v2 = earth.eph(e)
        d = sn.norm(r1, r2)
        if d < mn[1]:
            mn = (e, d)
        e = pk.epoch(e.mjd2000 + step)
    return mn

mn = earth_2018_av2_close_approach(e_start, e_end, 1.0)
print ((mn[0], mn[1] / pk.AU))
```
Prints `(2041-Sep-01 00:00:00, 0.035282046333263914)` i.e. according to the approximated orbits the minimum separation will be around 0.035 AU on 2041 September 1.

## Hohmann Transfer Window
Using pykep's Lambert problem solver we can evaluate candidate Hohmann transfer dates for feasibility:
``` py
# orbital_snoopy.py
def porkchop(o1, o2, t_begin_range, t_end_range):
    out = []
    i = 0.0
    while i < t_begin_range[2]:
        t1 = pk.epoch(t_begin_range[0].mjd2000 + i)
        outrow = []
        out.append((t1, outrow))
        j = 0.0
        while j < t_end_range[2]:
            t2 = pk.epoch(t_end_range[0].mjd2000 + j)
            dt = (t2.jd - t1.jd) * pk.DAY2SEC
            r1, v1 = o1.eph(t1)
            r2, v2 = o2.eph(t2)
            l = pk.lambert_problem(r1=r1, r2=r2, tof=dt, mu=_pk.MU_SUN, max_revs = 1)
            vb = l.get_v1()[0]
            ve = l.get_v2()[0]
            dv1 = norm(v1, vb)
            dv2 = norm(v2, ve)
            outrow.append((t2, dt, dv1 + dv2, dv1, dv2))

            j += t_end_range[1]

        i += t_begin_range[1]

    return out

def porkchop_2041_transfer():
    return (
        pk.planet.jpl_lp('earth'),
        orbit_2018_av2_pk(),
        (pk.epoch_from_string('2041-03-01 00:00:00'), 1.0, 90.1),
        (pk.epoch_from_string('2041-11-01 00:00:00'), 1.0, 90.1))

# user code
args = sn.porkchop_2041_transfer()
pork = sn.porkchop(*args)

print((pork[27][0],) + pork[27][1][29])
```
Prints `(2041-Mar-28 00:00:00, 2041-Nov-30 00:00:00, 21340800.0, 480.0902720093982, 192.937420833875, 287.15285117552315)`

The first two elements are the transfer departure and arrival dates.

## Limitations

## Potential Future Work

## References

### 1 "Pseudo-MPEC" for ZA9872D = 2018 AV2
[https://www.projectpluto.com/pluto/mpecs/2018av2.htm](https://www.projectpluto.com/pluto/mpecs/2018av2.htm) created 2020 Feb 11, retrieved 2022 Sep 16

### 2 pykep
pykep Python module, pykep Development Team. [http://esa.github.io/pykep/index.html](http://esa.github.io/pykep/index.html)
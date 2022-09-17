# Approximate Analysis of 2018 AV2 Orbit, Next Close Approach, Hohmann Transfer Window for TOMTOM/ Project Snoopy
Deskin Miller (deskinm at umich dot edu)  
September 2022

---
## Summary
By analyzing the orbital elements of 2018 AV2 we approximate the next closest Earth approach to be 2041 September 1 (+/- 1 day). We find an approximate Hohmann transfer window open between 2041 March 1 -- 2041 May 1 and arrival between 2041 November 1 -- 2042 January 1 with Earth-departure "$C_3 < 1 km^2/s^2$" and rendezvous "delta-V $< 500 m/s$". This can be used for evaluating project feasibility, further refinement of orbital dynamics, and initial constraints on spacecraft design.

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
        # Obtain ephemeris at time e
        r1, v1 = o.eph(e)
        r2, v2 = earth.eph(e)
        # Calculate the distance and update the minimum found
        d = sn.norm(r1, r2)
        if d < mn[1]:
            mn = (e, d)
        e = pk.epoch(e.mjd2000 + step)
    return mn

mn = earth_2018_av2_close_approach(e_start, e_end, 1.0)
print ((mn[0], mn[1] / pk.AU))
```
Prints `(2041-Sep-01 00:00:00, 0.035282046333263914)` i.e. according to the approximated orbits the minimum separation will be around 0.035 AU/ 5278119 km on 2041 September 1.

## Hohmann Transfer Window
Using pykep's Lambert problem solver we can evaluate candidate Hohmann transfer dates for feasibility:
``` py
# orbital_snoopy.py
def porkchop(o1, o2, t_begin_range, t_end_range):
    out = []
    i = 0.0
    while i < t_begin_range[2]:
        # t1 transfer departure date
        t1 = pk.epoch(t_begin_range[0].mjd2000 + i)
        outrow = []
        out.append((t1, outrow))
        j = 0.0
        while j < t_end_range[2]:
            # t2 transfer arrival date
            t2 = pk.epoch(t_end_range[0].mjd2000 + j)
            dt = (t2.jd - t1.jd) * pk.DAY2SEC
            # Calculate ephemeris for originating body at departure and
            # destination body at arrival
            r1, v1 = o1.eph(t1)
            r2, v2 = o2.eph(t2)
            l = pk.lambert_problem(r1=r1, r2=r2, tof=dt, mu=_pk.MU_SUN, max_revs = 1)
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

def porkchop_2041_transfer():
    return (
        pk.planet.jpl_lp('earth'),
        orbit_2018_av2_pk(),
        (pk.epoch_from_string('2041-03-01 00:00:00'), 1.0, 62.1),
        (pk.epoch_from_string('2041-11-01 00:00:00'), 1.0, 62.1))

# user code
args = sn.porkchop_2041_transfer()
pork = sn.porkchop(*args)

print((pork[27][0],) + pork[27][1][29])
```
Prints `(2041-Mar-28 00:00:00, 2041-Nov-30 00:00:00, 21340800.0, 480.0902720093982, 192.937420833875, 287.15285117552315)`

The first two elements are the transfer departure and arrival dates. The third element is the transfer time in seconds. The fourth element is the sum of the fifth and sixth; the fifth element is the Earth-escape (departure) velocity. The sixth element is the delta-V at 2018 AV2 arrival needed to match orbits.

## Limitations
This is a preliminary calculation with many known unknowns (and probably orders of magnitude more unknown unknowns). Some of the former:

 Caclulating as a purely Keplerian orbit does not arrive at a 2037 closest approach, which differs from what is stated in [[1]](#1-pseudo-mpec-for-za9872d--2018-av2). The source does not state its calculation methods, although it does provide a link to orbitalsimulator.com which arrives at a closest approach in 2037. This author makes no claim that the actual close approach is in 2041; rather the point is that simply knowing or being told a close approach date is not enough, and we need access to high-quality orbit prediction as the basis for mission planning.

Next, due to the relatively close orbit of 2018 AV2, the delta-V requirements are modest but the synodic period of 23.7 years is long. This long period leads to accumulating error either due to uncertainty in initial conditions or approximations in orbit prediction. Care must be taken for any mission design to account for the error bars.

One might hope for additional astronomical observations of 2018 AV2 which could reduce the error. However, this will likely prove quite challenging. The 2017-2018 observations were only within +/- two months of close approach, and likely a closer distance than will be achieved in 2041. The object is quite faint. A launch likely would need to precede re-acquisition, which means it won't have the benefit of updated observations until on-orbit.

Finally, the Hohmann transfer calculations disregard Earth gravity completely. A real launch and Earth-escape maneuver likely would occur later than the calculated departure date, since while under Earth's gravitational influence any craft would first have high velocity relative to Earth and then slow down to the asymptotic escape speed as it climbed out of Earth's gravity well.

## Potential Future Work
A higher accuracy orbital simulation would help immensely.

Given the inverse relationship between synodic period and transfer delta-V, a more thorough general study could be useful.

By looking at the existing observations of 2018 AV2 along with its best-fit orbit, it should be possible to place some likely constraints on the time and geometry of when and where to attempt astronomical observations of 2018 AV2 again.

Finally, from [[1]](#1-pseudo-mpec-for-za9872d--2018-av2) it's confidently stated that 2018 AV2 is not Snoopy, but is likely space junk. It would be worthwhile to investigate what their reasoning was for this conclusion. Some characteristics like area/ mass ratio may likely be calculated from public data about the lunar ascent module, for comparison with observations. Similarly a comparison with the observations made of J002E3, likely the Apollo 12 S-IVB stage, could be illuminating.

## References

### 1 "Pseudo-MPEC" for ZA9872D = 2018 AV2
[https://www.projectpluto.com/pluto/mpecs/2018av2.htm](https://www.projectpluto.com/pluto/mpecs/2018av2.htm) created 2020 Feb 11, retrieved 2022 Sep 16

### 2 pykep
pykep Python module, pykep Development Team. [http://esa.github.io/pykep/index.html](http://esa.github.io/pykep/index.html)
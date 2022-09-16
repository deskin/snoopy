# Copyright 2022 Deskin Miller
# Licensed to you under the MIT License. See COPYING for details

import pytest
from math import isclose

import orbital as orb
import orbital_snoopy as sn

class TestSnoopy:
		def test_earth_orbit_2017(self):
				o = sn.earth_orbit_2017()
				assert isclose(150855105915.9792, orb.utilities.norm(o.r))

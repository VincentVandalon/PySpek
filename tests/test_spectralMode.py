# Copyright 2012 Vincent Vandalon
#
# This file is part of <++>. <++> is free software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# <++> is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
# License for more details. You should have received a copy of
# the GNU General Public License along with <++>. If not, see
# <http://www.gnu.org/licenses/>. "

from PySpek import PySpek
from PySpek import ErroneousSpectrumType
import pytest


def test_size():
    """ Test correct opening and determination of size """
    spec = PySpek('1310-0304299.SPE')
    assert spec.get_size() == (1340, 100)
    spec = PySpek('1412-0401898.SPE')
    assert spec.get_size() == (335, 1)


def test_readSpec():
    """ Test the correct loading of the image """
    spec = PySpek('1412-0401898.SPE')
    x, y = spec.readSpec()
    assert len(x) == len(y)
    assert len(x) == 335  # Check correct readout
    xUncorrected, yUncorrected = spec.readSpec(perSecond=False)
    assert y[-1] * spec.aqTime == yUncorrected[-1]  # Check time scaling
    assert yUncorrected[-1] == 135  # Check correct readout
    assert abs(x[-1] - 664.47490223) < 0.00001  # idem
    with pytest.raises(ErroneousSpectrumType):
        specImage = PySpek('1310-0304299.SPE')
        specImage.readSpec()

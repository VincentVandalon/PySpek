# Copyright 2012 Vincent Vandalon
#
# This file is part of the NonlinearModel. NonlinearModel is
# free software: you can
# redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your
# option) any later version.
# NonlinearModel is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public
# License along with NonlinearModel.  If not, see
# <http://www.gnu.org/licenses/>.

# The image reading algorithm was inspired by
# http://scipy-cookbook.readthedocs.org/items/Reading_SPE_files.html
# This project aims to create a complete toolkit to read all types of SPE files

import numpy as N
import time

# Class to read out SPE files in spectral mode (assuming full binning
# in the vertical direction


class PySpek(object):

    def __init__(self, fname):
        self._fid = open(fname, 'rb')
        self._load_size()

    def _load_size(self):
        self._xdim = N.int64(self.read_at(42, 1, N.int16)[0])
        self._ydim = N.int64(self.read_at(656, 1, N.int16)[0])

    def _load_date_time(self):
        rawdate = self.read_at(20, 9, N.int8)
        rawtime = self.read_at(172, 6, N.int8)
        strdate = ''
        for ch in rawdate:
            strdate += chr(ch)
        for ch in rawtime:
            strdate += chr(ch)
        self._date_time = time.strptime(strdate, "%d%b%Y%H%M%S")

    def get_size(self):
        return (self._xdim, self._ydim)

    def read_ycal(self):
        return self.read_cal(3489)

    def read_xcal(self):
        return self.read_cal(3000)

    def read_cal(self, offset):
        self.xDimDet = self.read_at(6, 1, N.int16)
        self.xdim = self.read_at(42, 1, N.int16)
        self.polynom_coeff = self.read_at(offset+263, 6, N.double)
        self.datatype = self.read_at(108, 1, N.int16)
        self.yoffset = self.read_at(3489, 1, N.int16)
        self.noscans = self.read_at(34, 1, N.int16)

        # if self.read_at(4,1,N.int16) == -1:
        #   self.aqTime=self.read_at(660,1,N.int32)[0]/1000.
        # else:
        #   self.aqTime=self.read_at(4,1,N.int16)[0]/1000.

        self.aqTime = self.read_at(10, 1, N.float32)

        if self.aqTime == 0:
            self.aqTime = 1
            print('Time was not read',
                  'correctly from SPE file,',
                  'counts/s is not set correctly')

        if self.noscans == -1:
            self.noscans = self.read_at(664, 1, N.int32)

        self.nrFrames = self.read_at(1446, 1, N.int32)
        self.gain = self.read_at(198, 1, N.uint16)

        # Collect the ROI patterns
        self.NumROI = self.read_at(1510, 1, N.int16)
        self.RIOObjects = []
        for i in range(0, self.NumROI):
            roi = ROIObject()
            self.RIOObjects.append(roi)
            roi.startx = self.read_at(1512+12*i, 1, N.int16)
            roi.endx = self.read_at(1514+12*i, 1, N.int16)
            roi.groupx = self.read_at(1516+12*i, 1, N.int16)
            roi.effendx = roi.endx/roi.groupx

            roi.starty = self.read_at(1518+12*i, 1, N.int16)
            roi.endy = self.read_at(1520+12*i, 1, N.int16)
            roi.groupy = self.read_at(1522+12*i, 1, N.int16)
            roi.effendy = roi.endy/roi.groupy

    def convertPixels(self, pixels):
        wls = []
        currentROI = self.RIOObjects[0]
        for p in pixels:
            # Figure out in which ROI this pixel is
            if p > currentROI.effendx:
                currentROI = self.RIOObjects.next()
            p = currentROI.startx+currentROI.groupx*p
            wls.append(self.polynom_coeff[0]
                       + self.polynom_coeff[1]
                       * (p+1) + self.polynom_coeff[2]
                       * (p+1)**2)
        return N.array(wls)

    def read_at(self, pos, size, ntype):
        self._fid.seek(pos)
        return N.fromfile(self._fid, ntype, size)

    def _load_img(self, startFrame, endFrame):
        if startFrame == endFrame:
            imgSize = self._xdim * self._ydim
            if self.datatype == 0:
                img = self.read_at(4100, imgSize, N.float32)
            elif self.datatype == 1:
                img = self.read_at(4100, imgSize, N.int32)
            elif self.datatype == 2:
                img = self.read_at(4100, imgSize, N.int16)
            elif self.datatype == 3:
                img = self.read_at(4100, imgSize, N.uint16)
            else:
                raise Exception("Unknown datatype in header")

            x = self.convertPixels(N.arange(0, self._xdim))
            return x, img.reshape((self._ydim, self._xdim))
        else:
            frames = []
            for i in range(startFrame, endFrame):
                imgSize = self._xdim * self._ydim
                if self.datatype == 0:
                    # TODO: get size of float32 from python, NOT HARDCODED
                    img = self.read_at(4100+i*imgSize*4, imgSize, N.float32)
                elif self.datatype == 1:
                    # TODO: get size of int32 from python, NOT HARDCODED
                    img = self.read_at(4100+i*imgSize*4, imgSize, N.int32)
                elif self.datatype == 2:
                    img = self.read_at(4100+i*imgSize, imgSize, N.int16)
                elif self.datatype == 3:
                    img = self.read_at(4100+i*imgSize, imgSize, N.uint16)
                else:
                    raise Exception("Unknown datatype in header")

                frames.append(img.reshape((self._ydim, self._xdim))[0])
            x = self.convertPixels(N.arange(0, self._xdim))
            return x, frames

    def close(self):
        self._fid.close()

    def readSpec(self, startFrame=-1, endFrame=-1, perSecond=True):
        """ Read first spectrum in file or a sequence of spectra """
        if self.get_size()[1] != 1:
            raise ErroneousSpectrumType()

        self.read_xcal()
        if perSecond:
            x, y = self._load_img(startFrame, endFrame)
            y = y/self.aqTime
        else:
            x, y = self._load_img(startFrame, endFrame)
        # If single spectrum, strip the extra dimension
        if N.shape(y)[0] == 1:
            y = y[0]
        return x, y

    def readSpecs(self, perSecond=True):
        """ Read all spectra in file """
        if self.get_size()[1] != 1:
            raise ErroneousSpectrumType()
        self.read_xcal()  # needed to ensure access to nrFrames...
        return self.readSpec(0, self.nrFrames, perSecond)

    def readImage(self):
        """ Read file as image """
        if self.get_size()[1] < 2:
            raise ErroneousSpectrumType()
        self.read_xcal()
        x, img = self._load_img(-1, -1)
        return x, img


class ErroneousSpectrumType(Exception):
    pass


class ROIObject:
    i = 0

    def __init__(self):
        self.i += 1

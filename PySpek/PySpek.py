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

            x = self.convertPixels(N.arange(0, self._xdim*self._ydim))
            return x, img.reshape((self._ydim, self._xdim))[0]
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
            x = self.convertPixels(N.arange(0, self._xdim*self._ydim))
            return x, frames

    def close(self):
        self._fid.close()

    def readSpec(self, startFrame=-1, endFrame=-1, perSecond=True):
        if self.get_size()[1] !=1:
            raise ErroneousSpectrumType()

        """ Read some or 1 spectra in file """
        self.read_xcal()
        if perSecond:
            x, y = self._load_img(startFrame, endFrame)[::]
            y = y/self.aqTime
        else:
            x, y = self._load_img(startFrame, endFrame)[::]
        return x, y

    def readSpecs(self, perSecond=True):
        """ Read all spectra in file """
        if self.get_size()[1] !=1:
            raise ErroneousSpectrumType()
        self.read_xcal()
        return self.readSpec(0, self.nrFrames, perSecond)

    def readImage(self):
        """ Read file as image """
        if self.get_size()[1] < 2:
            raise ErroneousSpectrumType()
        # TODO: select right data type somehow
        # scale = 1
        self.read_xcal()
        # NOTE: axis are swapped
        x = self.convertPixels(N.arange(0, self._xdim))
        img = self.read_at(4100, self._xdim * self._ydim, N.uint16)
        print(self._xdim, self._ydim)
        return x, img.reshape((self._ydim, self._xdim))


class ErroneousSpectrumType(Exception):
    pass


class ROIObject:
    i = 0

    def __init__(self):
        self.i += 1

# ////////////////////////////////////////////////////////////////////////
# /Complete Header info
# /
# /                                    Decimal Byte
# /                                       Offset
# /                                    -----------
# /  SHORT   ControllerVersion              0  Hardware Version
# /  SHORT   LogicOutput                     2  Definition of Output BNC
# /  WORD    AmpHiCapLowNoise         4  Amp Switching Mode
# /  WORD    xDimDet                          6  Detector x dimension of chip.
# /  SHORT   mode                               8  timing mode
# /  float   exp_sec                               10  alternitive exposure, in sec.
# /  SHORT   VChipXdim                      14  Virtual Chip X dim
# /  SHORT   VChipYdim                      16  Virtual Chip Y dim
# /  WORD    yDimDet                         18  y dimension of CCD or detector.
# /  char    date[DATEMAX]                  20  date
# /  SHORT   VirtualChipFlag                30  On/Off
# /  char    Spare_1[2]                          32
# /  SHORT   noscan                            34  Old number of scans - should always be -1
# /  float   DetTemperature                    36  Detector Temperature Set
# /  SHORT   DetType                           40  CCD/DiodeArray type
# /  WORD    xdim                                42  actual # of pixels on x axis
# /  SHORT   stdiode                            44  trigger diode
# /  float   DelayTime                            46  Used with Async Mode
# /  WORD    ShutterControl                  50  Normal, Disabled Open, Disabled Closed
# /  SHORT   AbsorbLive                       52  On/Off
# /  WORD    AbsorbMode                    54  Reference Strip or File
# /  SHORT   CanDoVirtualChipFlag       56  T/F Cont/Chip able to do Virtual Chip
# /  SHORT   ThresholdMinLive              58  On/Off
# /  float   ThresholdMinVal                    60  Threshold Minimum Value
# /  SHORT   ThresholdMaxLive             64  On/Off
# /  float   ThresholdMaxVal                   66  Threshold Maximum Value
# /  SHORT   SpecAutoSpectroMode     70  T/F Spectrograph Used
# /  float   SpecCenterWlNm                  72  Center Wavelength in Nm
# /  SHORT   SpecGlueFlag                  76  T/F File is Glued
# /  float   SpecGlueStartWlNm             78  Starting Wavelength in Nm
# /  float   SpecGlueEndWlNm               82  Starting Wavelength in Nm
# /  float   SpecGlueMinOvrlpNm            86  Minimum Overlap in Nm
# /  float   SpecGlueFinalResNm            90  Final Resolution in Nm
# /  SHORT   PulserType                      94  0=None, PG200=1, PTG=2, DG535=3
# /  SHORT   CustomChipFlag               96  T/F Custom Chip Used
# /  SHORT   XPrePixels                       98  Pre Pixels in X direction
# /  SHORT   XPostPixels                    100  Post Pixels in X direction
# /  SHORT   YPrePixels                     102  Pre Pixels in Y direction
# /  SHORT   YPostPixels                    104  Post Pixels in Y direction
# /  SHORT   asynen                           106  asynchron enable flag  0 = off
# /  SHORT   datatype                         108  experiment datatype
# /                                             0 =   FLOATING POINT
# /                                             1 =   LONG INTEGER
# /                                             2 =   INTEGER
# /                                             3 =   UNSIGNED INTEGER //---- SY 03-22-2004 ???? should be UINT16 or WORD or SHORT ????
# /  SHORT   PulserMode                   110  Repetitive/Sequential
# /  USHORT  PulserOnChipAccums           112  Num PTG On-Chip Accums
# /  DWORD   PulserRepeatExp              114  Num Exp Repeats (Pulser SW Accum)
# /  float   PulseRepWidth                118  Width Value for Repetitive pulse (usec)
# /  float   PulseRepDelay                122  Width Value for Repetitive pulse (usec)
# /  float   PulseSeqStartWidth           126  Start Width for Sequential pulse (usec)
# /  float   PulseSeqEndWidth             130  End Width for Sequential pulse (usec)
# /  float   PulseSeqStartDelay           134  Start Delay for Sequential pulse (usec)
# /  float   PulseSeqEndDelay             138  End Delay for Sequential pulse (usec)
# /  SHORT   PulseSeqIncMode              142  Increments: 1=Fixed, 2=Exponential
# /  SHORT   PImaxUsed                    144  PI-Max type controller flag
# /  SHORT   PImaxMode                    146  PI-Max mode
# /  SHORT   PImaxGain                    148  PI-Max Gain
# /  SHORT   BackGrndApplied              150  1 if background subtraction done
# /  SHORT   PImax2nsBrdUsed              152  T/F PI-Max 2ns Board Used
# /  WORD    minblk                       154  min. # of strips per skips
# /  WORD    numminblk                    156  # of min-blocks before geo skps
# /  SHORT   SpecMirrorLocation[2]        158  Spectro Mirror Location, 0=Not Present
# /  SHORT   SpecSlitLocation[4]          162  Spectro Slit Location, 0=Not Present
# /  SHORT   CustomTimingFlag             170  T/F Custom Timing Used
# /  char    ExperimentTimeLocal[TIMEMAX] 172  Experiment Local Time as hhmmss\0
# /  char    ExperimentTimeUTC[TIMEMAX]   179  Experiment UTC Time as hhmmss\0
# /  SHORT   ExposUnits                   186  User Units for Exposure
# /  WORD    ADCoffset                    188  ADC offset
# /  WORD    ADCrate                      190  ADC rate
# /  WORD    ADCtype                      192  ADC type
# /  WORD    ADCresolution                194  ADC resolution
# /  WORD    ADCbitAdjust                 196  ADC bit adjust
# /  WORD    gain                         198  gain
# /  char    Comments[5][COMMENTMAX]      200  File Comments
# /  WORD    geometric                    600  geometric ops: rotate 0x01,
# /                                             reverse 0x02, flip 0x04
# /  char    xlabel[LABELMAX]             602  intensity display string
# /  WORD    cleans                       618  cleans
# /  WORD    NumSkpPerCln                 620  number of skips per clean.
# /  SHORT   SpecMirrorPos[2]             622  Spectrograph Mirror Positions
# /  float   SpecSlitPos[4]               626  Spectrograph Slit Positions
# /  SHORT   AutoCleansActive             642  T/F
# /  SHORT   UseContCleansInst            644  T/F
# /  SHORT   AbsorbStripNum               646  Absorbance Strip Number
# /  SHORT   SpecSlitPosUnits             648  Spectrograph Slit Position Units
# /  float   SpecGrooves                  650  Spectrograph Grating Grooves
# /  SHORT   srccmp                       654  number of source comp. diodes
# /  WORD    ydim                         656  y dimension of raw data.
# /  SHORT   scramble                     658  0=scrambled,1=unscrambled
# /  SHORT   ContinuousCleansFlag         660  T/F Continuous Cleans Timing Option
# /  SHORT   ExternalTriggerFlag          662  T/F External Trigger Timing Option
# /  long    lnoscan                      664  Number of scans (Early WinX)
# /  long    lavgexp                      668  Number of Accumulations
# /  float   ReadoutTime                  672  Experiment readout time
# /  SHORT   TriggeredModeFlag            676  T/F Triggered Timing Option
# /  char    Spare_2[10]                  678
# /  char    sw_version[FILEVERMAX]       688  Version of SW creating this file
# /  SHORT   type                         704   1 = new120 (Type II)
# /                                             2 = old120 (Type I )
# /                                             3 = ST130
# /                                             4 = ST121
# /                                             5 = ST138
# /                                             6 = DC131 (PentaMax)
# /                                             7 = ST133 (MicroMax/SpectroMax)
# /                                             8 = ST135 (GPIB)
# /                                             9 = VICCD
# /                                            10 = ST116 (GPIB)
# /                                            11 = OMA3 (GPIB)
# /                                            12 = OMA4
# /  SHORT   flatFieldApplied             706  1 if flat field was applied.
# /  char    Spare_3[16]                  708
# /  SHORT   kin_trig_mode                724  Kinetics Trigger Mode
# /  char    dlabel[LABELMAX]             726  Data label.
# /  char    Spare_4[436]                 742
# /  char    PulseFileName[HDRNAMEMAX]   1178  Name of Pulser File with
# /                                             Pulse Widths/Delays (for Z-Slice)
# /  char    AbsorbFileName[HDRNAMEMAX]  1298 Name of Absorbance File (if File Mode)
# /  DWORD   NumExpRepeats               1418  Number of Times experiment repeated
# /  DWORD   NumExpAccums                1422  Number of Time experiment accumulated
# /  SHORT   YT_Flag                     1426  Set to 1 if this file contains YT data
# /  float   clkspd_us                   1428  Vert Clock Speed in micro-sec
# /  SHORT   HWaccumFlag                 1432  set to 1 if accum done by Hardware.
# /  SHORT   StoreSync                   1434  set to 1 if store sync used
# /  SHORT   BlemishApplied              1436  set to 1 if blemish removal applied
# /  SHORT   CosmicApplied               1438  set to 1 if cosmic ray removal applied
# /  SHORT   CosmicType                  1440  if cosmic ray applied, this is type
# /  float   CosmicThreshold             1442  Threshold of cosmic ray removal.
# /  long    NumFrames                   1446  number of frames in file.
# /  float   MaxIntensity                1450  max intensity of data (future)
# /  float   MinIntensity                1454  min intensity of data (future)
# /  char    ylabel[LABELMAX]            1458  y axis label.
# /  WORD    ShutterType                 1474  shutter type.
# /  float   shutterComp                 1476  shutter compensation time.
# /  WORD    readoutMode                 1480  readout mode, full,kinetics, etc
# /  WORD    WindowSize                  1482  window size for kinetics only.
# /  WORD    clkspd                      1484  clock speed for kinetics & frame transfer
# /  WORD    interface_type              1486  computer interface
# /                                             (isa-taxi, pci, eisa, etc.)
# /  SHORT   NumROIsInExperiment         1488  May be more than the 10 allowed in
# /                                             this header (if 0, assume 1)
# /  char    Spare_5[16]                 1490
# /  WORD    controllerNum               1506  if multiple controller system will
# /                                             have controller number data came from.
# /                                             this is a future item.
# /  WORD    SWmade                      1508  Which software package created this file
# /  SHORT   NumROI                      1510  number of ROIs used. if 0 assume 1.
# /                                      1512 - 1630  ROI information
# /  struct ROIinfo { //---- SY 03-22-2004 ???? should be WORD or SHORT ????
# /   unsigned int startx                left x start value.
# /   unsigned int endx                  right x value.
# /   unsigned int groupx                amount x is binned/grouped in hw.
# /   unsigned int starty                top y start value.
# /   unsigned int endy                  bottom y value.
# /   unsigned int groupy                amount y is binned/grouped in hw.
# /  } ROIinfoblk[ROIMAX]                   ROI Starting Offsets:
# /                                                  ROI  1 = 1512
# /                                                  ROI  2 = 1524
# /                                                  ROI  3 = 1536
# /                                                  ROI  4 = 1548
# /                                                  ROI  5 = 1560
# /                                                  ROI  6 = 1572
# /                                                  ROI  7 = 1584
# /                                                  ROI  8 = 1596
# /                                                  ROI  9 = 1608
# /                                                  ROI 10 = 1620
# /  char    FlatField[HDRNAMEMAX]       1632  Flat field file name.
# /  char    background[HDRNAMEMAX]      1752  background sub. file name.
# /  char    blemish[HDRNAMEMAX]         1872  blemish file name.
# /  float   file_header_ver             1992  version of this file header
# /  char    YT_Info[1000]               1996-2996  Reserved for YT information
# /  LONG    WinView_id                  2996  == 0x01234567L if file created by WinX
# /
# /------------------------------------------------------------------------
# /
# /                        START OF X CALIBRATION STRUCTURE
# /
# /  double        offset                3000  offset for absolute data scaling
# /  double        factor                3008  factor for absolute data scaling
# /  char          current_unit          3016  selected scaling unit
# /  char          reserved1             3017  reserved
# /  char          string[40]            3018  special string for scaling
# /  char          reserved2[40]         3058  reserved
# /  char          calib_valid           3098  flag if calibration is valid
# /  char          input_unit            3099  current input units for
# /                                            "calib_value"
# /  char          polynom_unit          3100  linear UNIT and used
# /                                            in the "polynom_coeff"
# /  char          polynom_order         3101  ORDER of calibration POLYNOM
# /  char          calib_count           3102  valid calibration data pairs
# /  double        pixel_position[10]    3103  pixel pos. of calibration data
# /  double        calib_value[10]       3183  calibration VALUE at above pos
# /  double        polynom_coeff[6]      3263  polynom COEFFICIENTS
# /  double        laser_position        3311  laser wavenumber for relativ WN
# /  char          reserved3             3319  reserved
# /  unsigned char new_calib_flag        3320  If set to 200, valid label below
# /  char          calib_label[81]       3321  Calibration label (NULL term'd)
# /  char          expansion[87]         3402  Calibration Expansion area
# /
# /------------------------------------------------------------------------
# /
# /                        START OF Y CALIBRATION STRUCTURE
# /
# /  double        offset                3489  offset for absolute data scaling
# /  double        factor                3497  factor for absolute data scaling
# /  char          current_unit          3505  selected scaling unit
# /  char          reserved1             3506  reserved
# /  char          string[40]            3507  special string for scaling
# /  char          reserved2[40]         3547  reserved
# /  char          calib_valid           3587  flag if calibration is valid
# /  char          input_unit            3588  current input units for
# /                                            "calib_value"
# /  char          polynom_unit          3589  linear UNIT and used
# /                                            in the "polynom_coeff"
# /  char          polynom_order         3590  ORDER of calibration POLYNOM
# /  char          calib_count           3591  valid calibration data pairs
# /  double        pixel_position[10]    3592  pixel pos. of calibration data
# /  double        calib_value[10]       3672  calibration VALUE at above pos
# /  double        polynom_coeff[6]      3752  polynom COEFFICIENTS
# /  double        laser_position        3800  laser wavenumber for relativ WN
# /  char          reserved3             3808  reserved
# /  unsigned char new_calib_flag        3809  If set to 200, valid label below
# /  char          calib_label[81]       3810  Calibration label (NULL term'd)
# /  char          expansion[87]         3891  Calibration Expansion area
# /
# /                         END OF CALIBRATION STRUCTURES
# /
# /    ---------------------------------------------------------------------
# /
# /  char    Istring[40]                 3978  special intensity scaling string
# /  char    Spare_6[76]                 4018
# /  SHORT   AvGainUsed                  4094  avalanche gain was used
# /  SHORT   AvGain                      4096  avalanche gain value
# /  SHORT   lastvalue                   4098  Always the LAST value in the header
# / ******************************************************************************/
# /                        /* 4100 Start of Data */

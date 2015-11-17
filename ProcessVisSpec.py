#Copyright 2012 Vincent Vandalon
#
#This file is part of the NonlinearModel. NonlinearModel is free software: you can
#redistribute it and/or modify it under the terms of the GNU
#General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your
# option) any later version.
#NonlinearModel is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
#the GNU General Public License for more details.
#You should have received a copy of the GNU General Public
#License along with NonlinearModel.  If not, see
#<http://www.gnu.org/licenses/>.  
import numpy as N
import scipy as sp
import matplotlib.pyplot as plt
from scipy import stats
import os
from PySpek import *
import tkFileDialog #for selecting input file
import Tkinter #for root directory

plt.rcParams['lines.markersize']=7
plt.rcParams['axes.formatter.limits']=(-4, 4)

inputFiles=['1412-0401898.SPE']

for filName in inputFiles:
   ax=plt.subplot(111)
   basename,ext= os.path.splitext(filName)
   spec1=PySpek(filName)
   x,y=spec1.readSpec(perSecond=True)
   spec1.close()

   plt.plot(x,y,marker='',linestyle='-',color='red')
   plt.xlabel('Wavelength (nm)')
   plt.ylabel('Intensity ( $s^{-1}$pixel$^{-1}$)')
   plt.savefig('%s.pdf'%basename)

plt.clf()

filName='image.spe'
ax=plt.subplot(111)
basename,ext= os.path.splitext(filName)
image=PySpek(filName)
xdata,imdat=image.readImage()
image.close()
plt.imshow(imdat)
xgrid=N.arange(0,len(xdata))
plt.xticks(xgrid[::200],map(lambda s:'%.0f'%s,xdata[::200]))
plt.xlabel('Wavelength (nm)')
plt.savefig('%s.pdf'%"image")

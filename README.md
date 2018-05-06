# PySpek
Python class to interpret Princeton Instruments / Acton SPE files. It basically allows 
to parse the binary spectrum files into your python data processing scripts.

PySpek can handle SPE files containing spectra or images, even if multiple spectra are contained in a single SPE file.
The functionality for spectral, i.e. fully binning the chip perpendicular to the disperion direciton, 
has been tested extensively. The funcitonality for imaging is being developed and input is very welcome.

Typical usage will look something like this:
~~~
spec1=PySpek(filName)
x,y=spec1.readSpec(perSecond=True)
spec1.close()
~~~

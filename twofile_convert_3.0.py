#!/Python27/
import comtypes, comtypes.client
from ctypes import *
from comtypes.automation import *
import sys
import math
import numpy
from scipy.integrate import *
import csv
from numpy import ma
from matplotlib.widgets import Cursor, Button
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#__author__     = "Dr Berwyck Poad"
#__version__	= "3.0"
#__date__		= "20 January, 2013" 
#This script is designed to take two RAW files, extract the data and 
#subtract them from one another. This is intended for use with the manual  
#CD experiments (not employing polarisation switching). Version 3.0 pulls
#this program into code alignment with the BGSUB routine.
#
#Usage:
#------
#python twofile_convert_3.0.py <RAWFILE1.raw> <RAWFILE2.raw>
#
#Output:
#-------
#
#

 
#Default values for MSFileReader.GetMassListFromScanNumber call:
scanNum = 1
scanFilter = u''
scanIntensityCutoffType = 0 # 0 = none, 1=Abs, 2=Rel. to basepk
scanIntensityCutoffValue = 0
scanMaxNumberOfPeaks = 0
scanCentroidResult = 0
pl = VARIANT() #Unused variable
pf = VARIANT()
hm = c_double() #High Mass
lm = c_double() #Low Mass
# ml set up later
arsize = c_long()
points = c_long()
mf = VARIANT()

if (len(sys.argv) < 3):
	sys.exit(-1)
else:
	fInName = sys.argv[1]
	fInNam = sys.argv[2]

#~ if (len(sys.argv) < 3):
	#~ fOutName = fInName + ".csv"
#~ else:
	#~ fOutName = sys.argv[2]
 

 
# Set up the COM object interface
xr = comtypes.client.CreateObject('MSFileReader.XRawfile')
xr.open(fInName)
res = xr.SetCurrentController(0,1)
#print "res: " + str(res)
ns = c_long()
xr.GetNumSpectra(ns)
xr.GetLowMass(lm)
xr.GetHighMass(hm)
xr.GetMassListFromScanNum(  #long winded way to get number of data points
    c_long(1),scanFilter,
    c_long(scanIntensityCutoffType),
    c_long(scanIntensityCutoffValue),
    c_long(scanMaxNumberOfPeaks),
    c_long(scanCentroidResult),
    c_double(0),
    mf,pl,points            #points is all I want from this!
)
print "---------CD File Ppointsr--------------"
print "----------v 0.1 by BLJP--------------"
print "Input File 1:          "+fInName
print "Input File 2:          "+fInNam
print "Number of Spectra:     "+str(ns.value)
print "High Mass:             "+str(hm.value)
print "Low Mass:              "+str(lm.value)
print "Number of Data Points: "+str(points.value)
print "------------------------------------"
minmz=lm.value      # by default we want the whole scan change
maxmz=hm.value      # these values if you want less data
spectra=ns.value
 
#initialise our data array
dataodd = numpy.empty(shape=((ns.value/2), points.value), dtype=float)
dataoutosum = numpy.empty(shape=((ns.value/2), points.value-1), dtype=float)
my_array = numpy.empty(shape=(ns.value/2, points.value-1), dtype=float)
chrom_o = numpy.empty(shape=(2, ns.value/2) , dtype=float)		#chromatogram
#outStrings = []
print dataoutosum.shape
print dataout.shape	

#first we will extract the scan from the first file 
for i in range(1,ns.value+1):
	ml = VARIANT()
	xr.GetMassListFromScanNum(
		c_long(i),scanFilter,
		c_long(scanIntensityCutoffType),
		c_long(scanIntensityCutoffValue),
		c_long(scanMaxNumberOfPeaks),
		c_long(scanCentroidResult),
		c_double(0),
		ml,pl,arsize
		)
	dataout = numpy.array(ml.value)
	#print arsize.value
	#print dataout.shape
	dataoutosum[1]=dataoutosum[1]+dataout[1]
	chrom_o[0,i-1] = i
	chrom_o[1,i-1] = sum(dataout[1])
del(ml,pl)
xr.close()	

#For some silly reason, we have to recall the MSFileReader comtype again.
#At least some copy-pasta works for this!
xr = comtypes.client.CreateObject('MSFileReader.XRawfile')
xr.open(fInNam)
res = xr.SetCurrentController(0,1)
#print "res: " + str(res)
ns = c_long()
xr.GetNumSpectra(ns)
xr.GetLowMass(lm)
xr.GetHighMass(hm)
xr.GetMassListFromScanNum(  #long winded way to get number of data points
    c_long(1),scanFilter,
    c_long(scanIntensityCutoffType),
    c_long(scanIntensityCutoffValue),
    c_long(scanMaxNumberOfPeaks),
    c_long(scanCentroidResult),
    c_double(0),
    mf,pf,points            #points is all I want from this!
)

#now for the even scan numbers
dataoutesum = numpy.empty(shape=((ns.value), points.value), dtype=float)
chrom_e = numpy.empty(shape=(2, (ns.value + 1)), dtype=float)		#chromatogram
print dataoutesum.shape
 

for i in range(1,ns.value+1):
	me = VARIANT()
	xr.GetMassListFromScanNum(
		c_long(i),scanFilter,
		c_long(scanIntensityCutoffType),
		c_long(scanIntensityCutoffValue),
		c_long(scanMaxNumberOfPeaks),
		c_long(scanCentroidResult),
		c_double(0),
		me,pf,arsize
		)
	dataoute = numpy.array(me.value)
	#print dataoute.shape
	#print arsize.value
	dataoutesum[1]=dataoutesum[1]+dataoute[1]
	chrom_e[0,i-1] = i
	chrom_e[1,i-1] = sum(dataoute[1])



xr.close()      #close off the RAW file after use

#some more diagnostic bits...
#print dataoute.shape
#print dataoutesum.shape
#print max(dataoutosum[1])
#print max(dataoutesum[1])

#---------DATA PROCESSING------------------

dataoutosum[1] = dataoutosum[1]/(ns.value)			#average the odd spectra
#dataoutosum[1] = (dataoutosum[1]/sum(dataoutosum[1])) 	#normalise the odd column to total counts
#dataout[1] = ma.masked_where(dataout[1]>0.05*max(dataout[1]), dataout[1])
dataoutesum[1] = dataoutesum[1]/(ns.value)			#average even spectra
#dataoutesum[1] = (dataoutesum[1]/sum(dataoutesum[1])) 	#normalise the even column to total counts

data = numpy.empty(shape=(2, points.value), dtype=float)
#calculate the residual
data[1]=dataoutosum[1]-dataoutesum[1]

#------------PLOTTING SECTION----------------

#plot the CD mass spectrum
plt.figure(1, figsize=(11.7,8.3))	#size is lanscape A4 in inches
plt.suptitle('Subtraction: ' + fInName + ' - ' + fInNam)
plt.subplot(211)
plt.plot(dataout[0], data[1], color='k')
plt.title('CDMS Spectrum')
plt.axis([minmz, maxmz, min(data[1])+0.1*min(data[1]), max(data[1])+0.1*max(data[1])]) #plot the data, with the intensity max as the largest value
plt.xlabel('m/z')
#plt.xaxis.set_major_locator(MaxNLocator(10))
plt.ylabel('Residual Intensity')
plt.grid(True)

#plot the chromatogram
plt.subplot(212)
plt.title('TIC')
plt.plot(chrom_o[0], chrom_o[1], label=fInName, color='r')
plt.plot(chrom_e[0], chrom_e[1], label=fInNam, color='b')
plt.axis([0 , ns.value, 0, max(chrom_o[1])+0.1*max(chrom_o[1])]) 
plt.legend(loc='best')
plt.xlabel('Scan Number')
#plt.xaxis.set_major_locator(MaxNLocator(10))
plt.ylabel('Total Counts')
plt.grid(True)
#plt.show()
plt.savefig(fInName+"-"+fInNam+"processed.pdf", dpi=600)

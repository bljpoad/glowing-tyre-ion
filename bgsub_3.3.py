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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#__author__     = "Dr Berwyck Poad"
#__version__	= "3.3"
#__date__		= "25 January, 2013" 
#This script is designed to separate sequential scans within one RAW file
#and allow you to subtract them from one another. This is particularly 
#useful for doing laser on - laser off background subtraction.
#
#Usage:
#------
#python bgsub_3.0.py RAWFILE.raw 
#
#Output:
#-------
#A pdf file will be output with the subtracted spectra and their TICs
#The file will be called RAWFILE.rawprocessed.pdf
#

 
#Default values for MSFileReader.GetMassListFromScanNumber call:
scanFilter = u''				# unicode null string for filter type
scanIntensityCutoffType = 0 	# 0 = none, 1=Abs, 2=Rel. to basepk
scanIntensityCutoffValue = 0 	# 0 = none
scanMaxNumberOfPeaks = 0 		# 0 = All peaks
scanCentroidResult = 0			# 0 = no centroiding
pl = VARIANT() 					#<- Unused variables
pf = VARIANT()					#<- but it doesn't work without them
hm = c_double() 				# High Mass
lm = c_double() 				# Low Mass
arsize = c_long() 				# array sizes
points = c_long() 				# array size for setting thigs up
mf = VARIANT()					# Mass 'flags'

# This will stop the program if it gets no RAW file name
if (len(sys.argv) < 2):
	sys.exit(-1)
else:
	fInName = sys.argv[1]


# Set up the COM object interface
# Refer to the MSFileReader API for more details and 
# fun tools
xr = comtypes.client.CreateObject('MSFileReader.XRawfile')
xr.open(fInName)
res = xr.SetCurrentController(0,1) 
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
print
print
print "----------------------------------------"
print "|      Laser Background Subtract       |"
print "|         Version 3.3 by BLJP          |"
print "|        Date: 25 January 2013         |" 
print "----------------------------------------"
print
print "Input File:             "+fInName
print "Number of Spectra:      "+str(ns.value)
print "High Mass:              "+str(hm.value)
print "Low Mass:               "+str(lm.value)
print "Number of Data Points:  "+str(points.value)
print 
print
minmz=lm.value      # by default we want the whole scan. Change
maxmz=hm.value      # these values if you want less data
 
#initialise our data arrays
dataodd = numpy.empty(shape=((ns.value/2), points.value), dtype=float)
dataoutosum = numpy.empty(shape=((ns.value/2), points.value-1), dtype=float)
my_array = numpy.empty(shape=(ns.value/2, points.value-1), dtype=float)
chrom_o = numpy.empty(shape=(2, ns.value/2) , dtype=float)		#chromatogram
outStrings = []

# first we will extract the odd scan numbers 
# perhaps there is a more elegant way to do this...
for i in range(1,ns.value,2):
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
	dataodd = numpy.array(ml.value)
	rng = [(dataodd[0]>=minmz)&(dataodd[0]<=maxmz)]
	dataout = numpy.array((dataodd[0][rng],dataodd[1][rng]))
	dataoutosum[1]=dataoutosum[1]+dataout[1]
	#dud_array = dataout.transpose
	my_array[i/2] = dataout[1]
	chrom_o[0,i/2] = i
	chrom_o[1,i/2] = sum(dataout[1])

#now for the even scan numbers
dataeven = numpy.empty(shape=((ns.value/2), points.value), dtype=float)
dataoutesum = numpy.empty(shape=((ns.value/2), points.value-1), dtype=float)
chrom_e = numpy.empty(shape=(2, ns.value/2), dtype=float)		#chromatogram
my_even_array = numpy.empty(shape=(ns.value/2, points.value-1), dtype=float)

for i in range(2,ns.value+1,2):
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
	dataeven = numpy.array(me.value)
	rng = [(dataeven[0]>=minmz)&(dataeven[0]<=maxmz)]
	dataoute = numpy.array((dataeven[0][rng],dataeven[1][rng]))
	dataoutesum[1]=dataoutesum[1]+dataoute[1]
	my_even_array[i/2-1] = dataoute[1]
	#print dataoute[1]
	#print my_even_array[i/2-1]
	chrom_e[0,i/2-1] = i
	chrom_e[1,i/2-1] = sum(dataoute[1])

xr.close()      #close off the RAW file after use


#---------DATA PROCESSING------------------

#Figure out which trace is laser on and laser off
#We are assuming here that the laser off trace should
#have a higher integrated TIC
if sum(chrom_o[1]) > sum(chrom_e[1]):
	laser_on_c = chrom_e
	laser_off_c = chrom_o
	laser_on_ms = dataoutesum
	laser_off_ms = dataoutosum
	laser_on_array = numpy.transpose(my_even_array)
	laser_off_array = numpy.transpose(my_array)
else:
	laser_on_c = chrom_o
	laser_off_c = chrom_e
	laser_on_ms = dataoutosum
	laser_off_ms = dataoutesum
	laser_on_array = numpy.transpose(my_array)
	laser_off_array = numpy.transpose(my_even_array)
	
#Some statistics for later:
rsd_off = numpy.std(laser_off_c[1]) / numpy.average(laser_off_c[1]) 
rsd_on = numpy.std(laser_on_c[1]) / numpy.average(laser_on_c[1]) 
rsd_total = math.sqrt(pow(rsd_off,2) + pow(rsd_on,2)) * 100



recovery =  numpy.empty(shape=(1, ns.value/2), dtype=float)
recovery[0] = laser_on_c[1] / laser_off_c[1]

print 'Recovery mean:'+str(numpy.average(recovery[0])*100)
print 'Recovery median:'+str(numpy.median(recovery[0])*100)
print 'Recovery Std Dev:'+str(numpy.std(recovery[0])*100)
	
#Averaging 
laser_on_ms[1] = laser_on_ms[1]/(ns.value/2)				#average the laser_on_ms spectra
laser_off_ms[1] = laser_off_ms[1]/(ns.value/2)				#average laser off spectra

#Uncomment the next two lines if you want to normalise
#laser_on_ms[1] = (laser_on_ms[1]/max(laser_on_ms[1])) 		#normalise laser_on_ms to base peak
#laser_off_ms[1] = (laser_off_ms[1]/max(laser_off_ms[1])) 	#normalise laser off to base peak

std_dev_on = numpy.empty(shape=(points.value-1, 1), dtype=float)
for j in range (1, points.value-1):
	std_dev_on[j]=numpy.std(laser_on_array[j])

std_dev_on = numpy.transpose(std_dev_on)
	
#calculate the residual (laser on - laser off)
residual = numpy.empty(shape=(2, points.value-1), dtype=float)
residual_err = numpy.empty(shape=(1, points.value-1), dtype=float)
err_tic_on = numpy.empty(shape=(1, ns.value/2), dtype=float)
err_tic_off = numpy.empty(shape=(1, ns.value/2), dtype=float)

residual[0] = dataout[0]
residual[1] = laser_on_ms[1] - laser_off_ms[1]
residual_err[0] = math.sqrt(pow(numpy.std(laser_on_ms[1]),2) + pow(numpy.std(laser_off_ms[1]),2))
err_tic_on[0] = numpy.std(laser_on_ms[1])
err_tic_off[0] = numpy.std(laser_off_ms[1]) 



#calculate the recovery amount for both the parent and total ion count
base = [(laser_off_ms[1]>=max(laser_off_ms[1]))]
base_peak = numpy.array((residual[0][base],laser_off_ms[1][base]))
parent = [(residual[0]>(sum(base_peak[0]-1)))&(residual[0]<sum(base_peak[0])+1)]
parent_on = numpy.array((residual[0][parent], laser_on_ms[1][parent]))
parent_off = numpy.array((residual[0][parent], laser_off_ms[1][parent]))

#parent_recovery = numpy.empty(shape=(1, ns.value/2), dtype=float)
#parent_recovery[0] = parent_on[1] / parent_off[1]

print parent_on
print parent_off

#more feedback to the terminal window
print 'Parent detected at m/z: {:.2f}'.format(sum(base_peak[0]))
print 'Integration window:     {:.2f} +/- 1Da'.format(sum(base_peak[0]))
print
print
print '--------------STATISTICS---------------'
print 'Parent ion recovery:    {:.2%}'.format(sum(parent_on[1]) / sum(parent_off[1]))
print 'Total ion recovery :    {:.2%}' .format(sum(laser_on_c[1])/sum(laser_off_c[1]))
print 
print 'Laser Off Mean:         {:.2f}'  .format(numpy.average(laser_off_c[1]))
print '          RSD:          {:.2%}' .format(numpy.std([laser_off_c[1]]) / numpy.average([laser_off_c[1]]))
print 
print 'Laser On  Mean:         {:.2f}' .format(numpy.average(laser_on_c[1]))
print '          RSD:          {:.2%}' .format(numpy.std([laser_on_c[1]]) / numpy.average(laser_on_c[1]))
print 
print 'RSD Laser On - Off:     {:.2%}'  .format(rsd_total / 100)
print '----------------------------------------'

#------------PLOTTING SECTION----------------
# Use matplotlib.pyplot to plot out our results. 
# http://matplotlib.org/api/pyplot_api.html is an 
# awesome reference for what this can do.

#Figure out the largest value to plot in the MS
if math.fabs(min(residual[1])) > math.fabs(max(residual[1])):
	y_axis_limits = math.fabs(min(residual[1]))
else:
	y_axis_limits = math.fabs(max(residual[1]))

	
#PLOT MASS SPECTUM
plt.figure('Laser Background Subtraction', figsize=(11.7,8.3))
plt.suptitle('File: ' + fInName)
plt.subplot(211)
plt.plot(residual[0], residual[1])
#plt.errorbar(residual[0], residual[1], yerr=residual_err[0], fmt='--')
plt.title('Laser On - Laser Off Mass Spectrum')
plt.axis([minmz, maxmz, -(y_axis_limits+0.1*y_axis_limits), (y_axis_limits+0.1*y_axis_limits)])
plt.xlabel('m/z',style='italic')
#plt.xaxis.set_major_locator(MaxNLocator(10)) #why doesn't this work?
plt.ylabel('Intensity Difference (Counts)')
plt.grid(True)

#PLOT CHROMATOGRAM
plt.subplot(212)
plt.title('TIC')
plt.plot(laser_off_c[0], laser_off_c[1], label=r'laser off', color='b')
plt.plot(laser_on_c[0], laser_on_c[1], label=r'laser on', color='g')
plt.axis([0 , ns.value, 0, max(laser_off_c[1])+0.1*max(laser_off_c[1])]) 
plt.legend(loc='best')
plt.xlabel('Scan Number')
#plt.xaxis.set_major_locator(MaxNLocator(10)) #why doesn't this work?
plt.grid(True)
plt.ylabel('Total Counts')
plt.savefig(fInName+"plots.pdf", dpi=600)

#------------CSV OUTPUT SECTION----------------
#This will provide a CSV file with the stats as
#a headder and then the two mass spectra with the 
#residual. Useful if you want to plot the data in 
#your favorite plotting program

#build out output array 
spectra = numpy.empty(shape=(4, points.value-1), dtype=float)
spectra[0] = dataout[0] 			#m/z
spectra[1] = laser_on_ms[1]			#laser on spectrum
spectra[2] = laser_off_ms[1]		#laser off spectrum
spectra[3] = residual[1]			#laser on - laser off
spectra = numpy.transpose(spectra)	#we need to transpose it otherwise we get 4 long rows

outfile = open(fInName+'data.csv','wb',arsize.value)
#headder information for CSV file. 
outfile.write("Input File:            "+fInName +'\n')
outfile.write("Number of Spectra:     "+str(ns.value) +'\n')
outfile.write("High Mass:             "+str(hm.value)+'\n')
outfile.write("Low Mass:              "+str(lm.value)+'\n')
outfile.write("Parent detected at:    {:.2f}" .format(sum(base_peak[0]))+'\n')
outfile.write("Number of Data Points: "+str(points.value)+'\n')
outfile.write("Total Ion recovery:    {:.2%}" .format(sum(laser_on_c[1])/sum(laser_off_c[1]))+'\n')
outfile.write("Parent Ion recovery:   {:.2%}" .format(sum(parent_on[1])/sum(parent_off[1]))+'\n')
outfile.write("Laser On RSD:          {:.2%}" .format(numpy.std([laser_on_c[1]]) / numpy.average(laser_on_c[1]))+'\n')
outfile.write("Laser Off RSD:         {:.2%}" .format(numpy.std([laser_off_c[1]]) / numpy.average(laser_off_c[1]))+'\n\n')
outfile.write('m/z, Laser On, Laser Off, On - Off\n')
writer = csv.writer(outfile)

for row in spectra:
    writer.writerow(row)
outfile.close()

#Uncomment the following line for semi-interactive mode.
#This window will pop up after all output files have been saved. 
#plt.show()
	
print
print
print
print '   Laser Background Subtract finished'
print '----------------------------------------'
print 'Output files: '+ fInName+'plots.pdf'
print '              '+fInName+'data.csv'
print '----------------------------------------'
print 
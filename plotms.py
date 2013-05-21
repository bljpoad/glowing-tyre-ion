# Plotting routine for Mass Spectra acquired using duty-cycle acquisition mode.
# Date: 14 January, 2013
# Usage: python <this scripts name> <CSV file name>
# -------------------------------

import numpy as np
from numpy import ma
import csv
import sys
import matplotlib.pyplot as plt

if (len(sys.argv) < 2):
	sys.exit(-1)
else:
	fInName = sys.argv[1]

# Reading the file in assuming that there are 13 header lines and:
# data[0] is m/z
# data[1] is laser on MS
# data[2] is laser off MS
# data[3] is laser on - laser off
#
# Then we normalise the MS columns

data = np.genfromtxt(fInName, delimiter=',', skip_header=13, unpack=True)
#data[1]=(data[1]/max(data[1]))*100 #normalise the data column to total counts
#data[2]=(data[2]/max(data[2]))*100


#--------PLOTTING SECTION----------
plt.figure(figsize=(10,9))
plt.subplot(211)
#plt.suptitle('File: ' + fInName)
plt.plot(data[0], data[2]+1500, label = 'DL-TRP')
plt.plot(data[0], data[1], label = 'L-TRP')
plt.plot(data[0], data[3]+3000, label = 'D-TRP')
plt.axis([60, 220, min(data[1])-1, max(data[1])+0.05*max(data[1])]) #plot the data, with the intensity max as the largest value
plt.yticks([])
plt.xticks(np.arange(100, 301, 50))
plt.legend(loc = 'upper left')
plt.xlabel('m/z', style = 'italic')
plt.ylabel('Normalised Intensity (arb.)')
plt.grid(True)

#Residual plot, on a symmetric log axis for 'clarity'
# plt.subplot(212)
# plt.plot(data[0],data[3])
# plt.axis([min(data[0]), max(data[0]), min(data[3])+0.05*min(data[3]), max(data[3])+0.05*max(data[3])]) #plot the data, with the intensity max as the largest value
# plt.yscale('symlog')
# plt.xticks(np.arange(np.int_(min(data[0])), np.int_(max(data[0])), 25))
# plt.grid(True)
# plt.xlabel('m/z', style = 'italic')
# plt.ylabel('Absolute Residual Intensity (counts)')
plt.show()
#plt.savefig(fInName+'plot.pdf', dpi=600)
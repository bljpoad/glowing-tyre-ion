glowing-tyre-ion
================

Laser Chemistry Laboratory at UOW
---------------------------------
This is a repo containing scripts used by the Poad faction of the LCL. To use them, you will need the MSfile reader tools, available from Thermo Scientific at http://sjsupport.thermofinnigan.com/public/detail.asp?id=703. 

Programs:
BGSUB: main program to process interleaved laser on - laser off MS measurements within the same RAW file. Output is a PDF spectrum with the subtracted spectra and both the raw and rolling averaged chromatogram. 

TWOFILE_CONVERT: does the same as the above, but with two independent RAW files.

PLOTMS: This particular program is meant to plot out CSV data files that have been generated by the bgsub program. Useful for making figures for talks, reports, etc.

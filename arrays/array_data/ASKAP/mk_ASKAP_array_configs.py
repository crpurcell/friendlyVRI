#!/usr/bin/env python

telescope = "ASKAP"
latitude_deg = -26.697000722
diameter_m = 12.0

import os
import sys
from util_misc import ascii_dat_read


#-----------------------------------------------------------------------------#
def main():

    # Read the pad lookup table
    col, dummy = ascii_dat_read("ASKAP_pads.txt", delim=",",
                                       doFloatCols=[1, 2])
    padDict = {}
    pad = 1
    for E, N in zip(col[1], col[2]):
        padDict[pad] = (E-463334.04, N-7047071.01)
        pad += 1
    
    # Read the array configuration file
    col, dummy = ascii_dat_read("ASKAP_configs.txt", delim=" ")
    for confName, padStr in zip(col[1], col[2]):
        padLst = [int(x) for x in padStr.split(",")]

        outFileName = "ASKAP_%s.config" % confName
        FH = open(outFileName, "w")
        FH.write("#" + "-"*78 + "#\n")
        FH.write("#\n")
        FH.write("# Array definition file for the %s %s configuration.\n"
                 % (telescope, confName))
        FH.write("#\n")
        FH.write("#" + "-"*78 + "#\n")
        FH.write("\n")
        FH.write("# Name of the telescope\n")
        FH.write("telescope = %s\n" % telescope)
        FH.write("\n")
        FH.write("# Name of the configuration\n")
        FH.write("config = %s\n" % confName)
        FH.write("\n")
        FH.write("# Latitude of the array centre\n")
        FH.write("latitude_deg = %f\n" % latitude_deg)
        FH.write("\n")
        FH.write("# Antenna diameter\n")
        FH.write("diameter_m = %f\n" % diameter_m)
        FH.write("\n")
        FH.write("# Antenna coordinates (offset E, offset N)\n")

        for pad in padLst:
            FH.write("%f, %f\n" % (padDict[pad][0], padDict[pad][1]))      
        FH.close()
        
        
#-----------------------------------------------------------------------------#
main()

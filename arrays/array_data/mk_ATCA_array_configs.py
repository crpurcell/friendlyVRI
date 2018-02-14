#!/usr/bin/env python

telescope = "ATCA"
latitude_deg = -30.312906
diameter_m = 22.0

import os
import sys
from util_misc import ascii_dat_read


#-----------------------------------------------------------------------------#
def main():

    # Read the station lookup table
    col, dummy = ascii_dat_read("ATCA_stations.txt", delim=" ",
                                       doFloatCols=[2, 3])
    statDict = {}
    for station, N, W in zip(col[1], col[2], col[3]):
        statDict[station] = (-W+1622.449, N)

    # Read the array configuration file
    col, dummy = ascii_dat_read("ATCA_configs.txt", delim=" ",
                                doFloatCols=[2, 3, 4, 5, 6, 7])

    for confName, A1, A2, A3, A4, A5, A6 in zip(col[1], col[2], col[3], col[4],
                                                col[5], col[6], col[7]):

        if A1=='':
            continue
        outFileName = "ATCA_%s.config" % confName
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
        FH.write("%f, %f\n" % (statDict[A1][0], statDict[A1][1]))
        FH.write("%f, %f\n" % (statDict[A2][0], statDict[A2][1]))
        FH.write("%f, %f\n" % (statDict[A3][0], statDict[A3][1]))
        FH.write("%f, %f\n" % (statDict[A4][0], statDict[A4][1]))
        FH.write("%f, %f\n" % (statDict[A5][0], statDict[A5][1]))
        FH.write("%f, %f\n" % (statDict[A6][0], statDict[A6][1]))        
        FH.close()
        
    for confName, A1, A2, A3, A4, A5 in zip(col[1], col[2], col[3], col[4],
                                            col[5], col[6]):

        if A1=='':
            continue
        confName += "_No_6"
        outFileName = "ATCA_%s.config" % confName
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

        FH.write("%f, %f\n" % (statDict[A1][0], statDict[A1][1]))
        FH.write("%f, %f\n" % (statDict[A2][0], statDict[A2][1]))
        FH.write("%f, %f\n" % (statDict[A3][0], statDict[A3][1]))
        FH.write("%f, %f\n" % (statDict[A4][0], statDict[A4][1]))
        FH.write("%f, %f\n" % (statDict[A5][0], statDict[A5][1]))
        
        FH.close()

#-----------------------------------------------------------------------------#
main()

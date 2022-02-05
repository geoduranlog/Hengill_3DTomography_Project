#!/usr/bin/env python2

import re
import math
import time
import numpy as np
import cPickle
#from obspy.core import UTCDateTime
#from obspy.core import read
#from obspy.core.util import gps2DistAzimuth
from obspy.geodetics import gps2dist_azimuth   #  ADN
#from obspy.taup.taup import getTravelTimes
#from obspy.arclink.client import Client
#from obspy.fdsn import Client as Client_WS
#import pg  #ADN
import os
import shutil
import glob
import string
import matplotlib.pyplot as plt
import optparse

global VERY_SMALL_DOUBLE, SMALL_DOUBLE
VERY_SMALL_DOUBLE = 1.0e-30
SMALL_DOUBLE = 1.0e-8

#------------------------------------------------------------------------------------------
class stationobj(object):

      def __init__(self,staorg,staali,net,loc,lat,lon,ele,dep,SCP,SCS,noP,noS,ICC):
          self.staorg = str(staorg)
          self.staali = str(staali)
          self.net = str(net)
          self.loc = str(loc)
          self.lat = float(lat)
          self.lon = float(lon)
          self.ele = float(ele) #in km
          self.dep = float(dep)
          self.SCP = float(SCP)
          self.SCS = float(SCS)
          self.noP = int(noP)
          self.noS = int(noS)
          self.ICC = int(ICC)      #Velest specific number, same ICC -> station merged; highest ICC -> reference station.

      def __str__(self):
          return '%-7s %-7s %-7s %-7s %8.4f %9.4f %6.3f %6.3f' % (self.staorg,self.staali,self.net,self.loc,self.lat,self.lon,self.ele,self.dep)
#------------------------------------------------------------------------------------------
#class stationobj(object):
#
#      def __init__(self,staorg,staali,net,loc,lat,lon,ele):
#          self.staorg = str(staorg)
#          self.staali = str(staali)
#          self.net = str(net)
#          self.loc = str(loc)
#          self.lat = float(lat)
#          self.lon = float(lon)
#          self.ele = float(ele)
#
#      def __str__(self):
#          return '%-7s %-7s %-7s %-7s %8.4f %9.4f %6.3f' % (self.staorg,self.staali,self.net,self.loc,self.lat,self.lon,self.ele)
#
#------------------------------------------------------------------------------------------
def StringIsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
#------------------------------------------------------------------------------------------
def LoadSIMULPSStationFile(ifile,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,phase,whitelist):

    statlist = []

    icc = 0

    #Open
    fp = open(ifile,"r")

    for line in fp:

        if(line[9:10]!="N") and (line[9:10]!="S"):
           print "    Ignore line:",line
           continue
        if(line[20:21]!="E") and (line[20:21]!="W"):
           print "    Ignore line:",line
           continue

        #stat,netw,lat,lon,ele,nobsP,nobsS        
        net = ''
        loc = ''
        staorg = line[2:6]
        staali = line[2:6] #Alias = original
        lat  = float(line[7:9])
        latc = line[9:10]
        latm = float(line[11:16])
        lon  = float(line[17:20])
        lonc = line[20:21]
        lonm = float(line[22:27])
        ele = float(line[27:32])  #elev in m..
        if(phase=='P'):
           stcP = float(line[32:37])
           stcS = 0
        if(phase=='S'):
           stcP = 0
           stcS = float(line[32:37])

        #Check if number of observations is provided:
        nobsP=-9
        nobsS=-9

        #Convert coordinates:
        lat = lat + latm/60.0
        lon = lon + lonm/60.0

        if(latc == 'S'):
           lat = lat*-1.0

        if(lonc == 'W'):
           lon = lon*-1.0

        #Check if station was read before:
        for j in range(len(statlist)):
            if(statlist[j].staorg == staorg):
               print "WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME:",statlist[j].staorg,statlist[j].net
               fp1.write("WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME: %-s %-s\n" % (statlist[j].staorg,statlist[j].net))
            if(statlist[j].staali == staali):
                print "WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME:",statlist[j].staali,statlist[j].net
                fp1.write("WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME: %-s %-s\n" % (statlist[j].staali,statlist[j].net))

        #Check if nobs provided AND enough observations:
        #if(nobsP>=0) and (nobsP<minNobs):
        #   print "WARNING: Number of observations not sufficient for station",staorg,nobsP,'--> skip station'
        #   continue

        #Now reset nobsP in case information was not provided:
        if(nobsP<0):
           nobsP = 0
        if(nobsS<0):
           nobsS = 0

        #Check if station whitelisted:
        if(len(whitelist)>0):
           whitelisted = False
           for z in range(len(whitelist)):
               if(whitelist[z] == staorg) or (whitelist[z] == staali):
                  whitelisted = True
                  break
        else:
           whitelisted = True

        #Check if station is in box:
        if(lat>=latmin) and (lat<=latmax) and (lon>=lonmin) and (lon<=lonmax) and (whitelisted):
           #Append station: staorg,staali,net,loc,lat,lon,ele,dep,SCP,SCS
           #print "--> Select Station",staorg,staali
           icc += 1
           statlist.append(stationobj(staorg,staali,net,loc,lat,lon,ele/1000.0,0.0,stcP,stcS,nobsP,nobsS,icc))

    fp.close()

    return statlist
#------------------------------------------------------------------------------------------
def LoadVELESTStationFile(ifile,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,minNobs,phase,whitelist):

    statlist = []

    icc = 0

    #Open
    fp = open(ifile,"r")

    for line in fp:

        if(line[0:1]=="("):
           print "    Ignore line:",line
           continue
        if(len(line)<=46):
           print "    Ignore line:",line
           continue

        #stat,netw,lat,lon,ele,nobsP,nobsS
        net = ''
        loc = ''
        staorg = line[0:4]
        staali = line[0:4] #Alias = original
        lat  = float(line[4:11])
        latc = line[11:12]
        lon  = float(line[13:21])
        lonc = line[21:22]
        ele = float(line[23:28])  #elev in m..

        if(phase=='P'):
           stcP = float(line[35:40])
           stcS = 0.0

        if(phase=='PS'):
           stcP = float(line[35:40])
           stcS = float(line[42:47])
 
        if(phase=='SA'):
           stcS = float(line[35:40])
           stcP = 0.0
           #print "--> Phase-Type:",phase,"StatCor-S:",stcS

        #Check if number of observations is provided:
        if(len(line)>=65) and (StringIsInt(line[58:59])) and (StringIsInt(line[64:65])):
           if(phase=='PS'):
              nobsP = int(line[54:59])
              nobsS = int(line[60:65])
           if(phase=='P'):
              nobsP = int(line[54:59])
              nobsS = 0
           if(phase=='SA'):
              nobsS = int(line[54:59])
              nobsP = int(line[54:59])  #We haveto asigne S-obs to P in order to pass the test of minimum phases (only performed for P)
        else:
           nobsP=-9
           nobsS=-9

        if(latc == 'S'):
           lat = lat*-1.0

        if(lonc == 'W'):
           lon = lon*-1.0

        #Check if station was read before:
        for j in range(len(statlist)):
            if(statlist[j].staorg == staorg):
               print "WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME:",statlist[j].staorg,statlist[j].net
               fp1.write("WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME: %-s %-s\n" % (statlist[j].staorg,statlist[j].net))
            if(statlist[j].staali == staali):
                print "WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME:",statlist[j].staali,statlist[j].net
                fp1.write("WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME: %-s %-s\n" % (statlist[j].staali,statlist[j].net))

        #Check if nobs provided AND enough observations:
        if(nobsP>=0) and (nobsP<minNobs):
           print "WARNING: Number of observations not sufficient for station",staorg,nobsP,phase,'--> skip station'
           continue

        #Now reset nobsP in case information was not provided:
        if(nobsP<0):
           nobsP = 0
        if(nobsS<0):
           nobsS = 0

        #Check if station whitelisted:
        if(len(whitelist)>0):
           whitelisted = False
           for z in range(len(whitelist)):
               if(whitelist[z] == staorg) or (whitelist[z] == staali):
                  whitelisted = True
                  break
        else:
           whitelisted = True

        #Check if station is in box:
        if(lat>=latmin) and (lat<=latmax) and (lon>=lonmin) and (lon<=lonmax) and (whitelisted):
           #Append station: staorg,staali,net,loc,lat,lon,ele,dep,SCP,SCS
           #print "--> Select Station",staorg,staali
           icc += 1
           statlist.append(stationobj(staorg,staali,net,loc,lat,lon,ele/1000.0,0.0,stcP,stcS,nobsP,nobsS,icc))

    fp.close()

    return statlist
#------------------------------------------------------------------------------------------
def loadstations(ifile,dbgmode,latmin,latmax,lonmin,lonmax):

    list = []

    icc=0

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<5):
           print "Line not considered:", line

        else:
           #Test if line is comment:
           t = string.split(line,None,-1)
           if(t[0].startswith("#") == False):
              staorg = string.strip(line[0:5])
              staali = string.strip(line[9:13])
              lat    = float(line[15:23])
              lon    = float(line[23:32])
              ele    = float(line[33:39])
              net    = ""
              loc    = ""

              #Check if station is in box:
              if(lat>=latmin) and (lat<=latmax) and (lon>=lonmin) and (lon<=lonmax):
                 list.append(stationobj(staorg,staali,net,loc,lat,lon,ele,1,1))
                 if(dbgmode == 1):
                    icc += 1
                    print stationobj(staorg,staali,net,loc,lat,lon,ele,1,1,icc)

    #Close
    fp.close()

    return list
#------------------------------------------------------------------------------------------
def LoadSC3exStationFile(ifile,dbgmode,latmin,latmax,lonmin,lonmax,logfile_DUPL,whitelist):
    list = []

    icc=0

    #Open
    fp = open(ifile,"r")

    fp1 = open(logfile_DUPL,"w")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<5):
           print "Line not considered:", line

        else:
           #Test if line is comment:
           t = string.split(line,None,-1)
           if(t[0].startswith("#") == False):
              staorg = string.strip(t[1])
              staali = string.strip(t[1])
              #print "Load station",staorg,staali
              lat    = float(t[2])
              lon    = float(t[3])
              ele    = float(t[4])
              net    = str(t[0])
              loc    = ""
              stcP   = float('0.0')
              stcS   = float('0.0')

              #Check if station was read before:
              for j in range(len(list)):
                  if(list[j].staorg == staorg):
                     print "WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME:",list[j].staorg,list[j].net
                     fp1.write("WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME: %-s %-s\n" % (list[j].staorg,list[j].net))
                  if(list[j].staali == staali):
                     print "WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME:",list[j].staali,list[j].net
                     fp1.write("WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME: %-s %-s\n" % (list[j].staali,list[j].net))

              #Check if station whitelisted:
              if(len(whitelist)>0):
                 whitelisted = False
                 for z in range(len(whitelist)):
                     if(whitelist[z] == staorg) or (whitelist[z] == staali):
                        whitelisted = True
                        break
              else:
                 whitelisted = True

              #Check if station is in box:
              if(lat>=latmin) and (lat<=latmax) and (lon>=lonmin) and (lon<=lonmax) and (whitelisted):
                 print "--> Select Station",staorg,staali
                 icc += 1
                 list.append(stationobj(staorg,staali,net,loc,lat,lon,ele,0.0,stcP,stcS,1,1,icc))
                 if(dbgmode == 1):
                    print stationobj(staorg,staali,net,loc,lat,lon,ele,0.0,stcP,stcS,1,1,icc)

    #Close
    fp.close()
    fp1.close()

    return list
#------------------------------------------------------------------------------------------
def LoadWhiteList(ifile):

    whitelist = []

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<1):
           print "Line not considered:", line

        else:
           #Test if line is comment:
           t = string.split(line,None,-1)
           if(t[0].startswith("#") == False):
              whitelist.append(string.strip(t[0]))
              continue

    fp.close()

    return whitelist
#------------------------------------------------------------------------------------------
def LoadAliasStationFile(ifile,dbgmode,latmin,latmax,lonmin,lonmax,logfile_DUPL,whitelist):

    icc=0

    list = []

    #Open
    fp = open(ifile,"r")

    fp1 = open(logfile_DUPL,"w")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<5):
           print "Line not considered:", line

        else:
           #Test if line is comment:
           t = string.split(line,None,-1)
           if(t[0].startswith("#") == False):
              staorg = string.strip(line[0:5])
              staali = string.strip(line[9:13])
              #print "Load station",staorg,staali
              lat    = float(line[15:23])
              lon    = float(line[23:32])
              ele    = float(line[33:39])
              net    = str(line[40:53])
              loc    = ""
              stcP   = float(line[54:61])
              stcS   = float(line[62:69])

              #Check if station was read before:
              for j in range(len(list)):
                  if(list[j].staorg == staorg):
                     print "WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME:",list[j].staorg,list[j].net
                     fp1.write("WARNING: DUBPLICAT ORIGINAL STATION <-> ORIGINAL NAME: %-s %-s\n" % (list[j].staorg,list[j].net))
                  if(list[j].staali == staali):
                     print "WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME:",list[j].staali,list[j].net
                     fp1.write("WARNING: DUBPLICAT ALIAS    STATION <-> ALIAS    NAME: %-s %-s\n" % (list[j].staali,list[j].net))

              #Check if station whitelisted:
              if(len(whitelist)>0):
                 whitelisted = False
                 for z in range(len(whitelist)):
                     if(whitelist[z] == staorg) or (whitelist[z] == staali):
                        whitelisted = True
                        break
              else:
                 whitelisted = True

              #Check if station is in box:
              if(lat>=latmin) and (lat<=latmax) and (lon>=lonmin) and (lon<=lonmax) and (whitelisted):
                 print "--> Select Station",staorg,staali
                 icc += 1
                 list.append(stationobj(staorg,staali,net,loc,lat,lon,ele,0.0,stcP,stcS,1,1,icc))
                 if(dbgmode == 1):
                    print stationobj(staorg,staali,net,loc,lat,lon,ele,0.0,stcP,stcS,1,1,icc)

    #Close
    fp.close()
    fp1.close()

    return list
#------------------------------------------------------------------------------------------
def merge_StationlistPS(stationlist,stationlistSdelay):

    #staorg,staali,net,loc,lat,lon,ele,dep,SCP,SCS,noP,noS,ICC

    #Loop over possible S-delays:
    for i in range(len(stationlistSdelay)):
          
        #Loop over all P-delays:
        for j in range(len(stationlist)):

            #Check if stations match:
            if(stationlistSdelay[i].staorg == stationlist[j].staorg):
               stationlist[j].SCS = stationlistSdelay[i].SCS
               stationlist[j].noS = stationlistSdelay[i].noS
               break

    return stationlist
#------------------------------------------------------------------------------------------
##################################################################################################################
#Main
#Get command line arguments:
oparser = optparse.OptionParser()

oparser.add_option('-s', '--stationlist',        action='store', dest='statlist',    help='Station list (MPX format; default: /Users/tdiehl/lib/sed_stations.GSE_SED.alias)')
oparser.add_option('-i', '--inputformat',        action='store', dest='inform',      help='Format of input station list; default: aliaslist, SC3exp (SC3DB_ExtractEvents.py), VELEST, SIMULPS, other to be implemented')
oparser.add_option('-o', '--output',             action='store', dest='outfi',       help='Output file with converted coordinates: convert_station.out (default)')
oparser.add_option('-e', '--exportformat',       action='store', dest='outform',     help='Output format: VELEST (default), SIMULPS, NLL_Grid2Time, NLL_StatDelay, NLL_LocAlias, hypoDD')
oparser.add_option('-c', '--complementarySdelays', action='store', dest='SdelayFile',help='Complementary S-delays from a SIMULPS station file (S as P); default: None')
oparser.add_option('-v', '--complementarySdelaysVELEST', action='store', dest='SdelayFileVELEST',help='Complementary S-delays from a VELEST station file (S as P); default: None')
oparser.add_option('-a', '--aliaslist', action='store', dest='aliaslist',help='ALIAS-File, is given, station names of VELEST/SIMULPS lists are replaced with originalnam; default: None')
oparser.add_option('-w', '--whitelist',      action='store', dest='whitelistFile'   ,help='File with whitelist of stations (original or alias) to extract from station list; default: None')

(options, arg) = oparser.parse_args()

#List with station coordinates:
#Check if list is included in list of input arguments:
if(options.statlist != None):
   statfile = options.statlist
else:
   statfile = "/Users/tdiehl/lib/sed_stations.GSE_SED.alias"

#Check if whitelist is defined:
if(options.whitelistFile != None):
   whitelistFile = options.whitelistFile
else:
   whitelistFile = "None"

#Check if complementary SIMULPS file defined:
if(options.SdelayFile != None):
   SdelayFile = options.SdelayFile
else:
   SdelayFile = "None"

#Check if complementary VELEST file defined:
if(options.SdelayFileVELEST != None):
   SdelayFileVELEST = options.SdelayFileVELEST
else:
   SdelayFileVELEST = "None"

#Input format:
#Check if format is included in list of input arguments:
if(options.inform != None):
   inform = options.inform
else:
   inform = "aliaslist"

#Output file:
#Check if output file is included in list of input arguments:
if(options.outfi != None):
   outfi = options.outfi
else:
   outfi='convert_station.out'
print "--> Output    File       :",outfi

#Check output format:
#Check if output file is included in list of input arguments:
if(options.outform != None):
   outform = options.outform
else:
   outform = 'VELEST'
print "--> Output    Format     :",outform

#Check if complementary ALIAS file defined:
if(options.aliaslist != None):
   aliasfile = options.aliaslist
else:
   aliasfile = "None"

#Define Region:
#region = 'CH-Tomo-Update'
#region = 'Alparray'
region = 'COSEISMIQ'

#To merge stations with virtually same location, define maximum horizontal distance (m):
maxdistH = 100.0

#For VELEST sta file with nobs specified, only extract stations with minimum nobs:
minNobsVel = 1

########################################################################
print "Regional settings for:",region

if(region == 'CH-Tomo-Update'):
   #geographic selection (New CH-Tomo model):
   latmin = 44.5
   latmax = 49.5
   lonmin = 4.5
   lonmax = 12.1 #To include RISI

   #Short distance conversion point used for SIMULPS (SIMULPS/VELEST convention: East is negative):
   sd_lat = 47.00000
   sd_lon = -8.5000

if(region == 'Alparray'):
   #geographic selection (now restriction at the moment):
   latmin = 25.0
   latmax = 70.0
   lonmin = -10.0
   lonmax = 40.0

   #Short distance conversion point used for SIMULPS (SIMULPS/VELEST convention: East is negative):
   sd_lat = 47.00000
   sd_lon = -8.5000

if(region=="COSEISMIQ"):
   #geographic selection (now restriction at the moment):
   latmin = 60
   latmax = 70
   lonmin = -25
   lonmax = -15

   #Short distance conversion point used for SIMULPS (SIMULPS/VELEST convention: East is negative):
   sd_lat = 64.02
   sd_lon = +21.35

#Log-files
logfile_DUPLICA = outfi + '.DuplicateStations'
logfile_ALIASDUPLICA = outfi + '.DuplicateStations.aliasfile'

#Read whiletlist:
if(whitelistFile != "None") and (os.path.isfile(whitelistFile)):
   print "--> Load Station Whitelist (stations not included will not be extracted)",whitelistFile
   whitelist = LoadWhiteList(whitelistFile)
else:
   whitelist = []

#Read alias file:
aliaslist = []
if(aliasfile != "None") and (os.path.isfile(aliasfile)):
   aliaslist = LoadAliasStationFile(aliasfile,0,latmin,latmax,lonmin,lonmax,logfile_ALIASDUPLICA,whitelist)

#Read stations:
if(inform == 'aliaslist'):
   stationlist = LoadAliasStationFile(statfile,0,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,whitelist)
if(inform == 'SC3exp'):
   stationlist = LoadSC3exStationFile(statfile,0,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,whitelist)
if(inform == 'VELEST'):
   stationlist = LoadVELESTStationFile(statfile,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,minNobsVel,'PS',whitelist)
if(inform == 'VELEST') and (SdelayFileVELEST != "None") and (os.path.isfile(SdelayFileVELEST)):
   stationlistSdelay = LoadVELESTStationFile(SdelayFileVELEST,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,minNobsVel,'SA',whitelist)
   stationlist = merge_StationlistPS(stationlist,stationlistSdelay)
if(inform == 'SIMULPS'):
   stationlist = LoadSIMULPSStationFile(statfile,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,'P',whitelist)
if(inform == 'SIMULPS') and (SdelayFile != "None") and (os.path.isfile(SdelayFile)):
   print "--> Load S-delays from additional file:",SdelayFile
   stationlistSdelay = LoadSIMULPSStationFile(SdelayFile,latmin,latmax,lonmin,lonmax,logfile_DUPLICA,'S',whitelist)
   stationlist = merge_StationlistPS(stationlist,stationlistSdelay)

#Replace station names with original names in potential alias list:
aliasmis = 0
if(len(aliaslist)>0) and (outform != "SIMULPS") and (outform != "VELEST"):
   #Loop over stations list
   for i in range(len(stationlist)):
       sidx = -1
       for j in range(len(aliaslist)):
           if(stationlist[i].staali == aliaslist[j].staali):
              sidx = j
              break
    
       #Now replace or warn if station is not included in alias file:
       if(sidx>=0):
          print "  --> Replace alias name",stationlist[i].staorg,"by original name",aliaslist[sidx].staorg
          stationlist[i].staorg=aliaslist[sidx].staorg 
       else:
          aliasmis += 1
          print "  --> ERROR:  alias name missing in alias2org list:",stationlist[i].staali


#Write basic event information + converted OT:
fp1 = open(outfi,"w")

#Write possible header:
if(outform == "VELEST"):
   fp1.write("%-s\n" % ("(a4,f7.4,a1,1x,f8.4,a1,1x,i5,1x,i1,1x,i3,1x,f5.2,2x,f5.2)"))

if(outform == "SIMULPS"):
   sd_latd = int(sd_lat)
   sd_latm = abs(sd_lat-sd_latd)*60.0
   sd_lond = int(sd_lon)
   sd_lonm = abs(sd_lon-sd_lond)*60.0
   fp1.write("       %03d  %05.2f %04d %05.2f   0    0     0\n" % (sd_latd,sd_latm,sd_lond,sd_lonm))
   fp1.write("%4d\n" % (len(stationlist)))

#Loop over all stations in list:
for i in range(len(stationlist)):

    #First calculate the epicentral distance to all previous stations:
    #distinf = gps2dist_azimuth(elat,elon,stations[statidx].lat,stations[statidx].lon)
    if(i >= 1):
       for j in range(i):
           dist=gps2dist_azimuth(stationlist[i].lat,stationlist[i].lon,stationlist[j].lat,stationlist[j].lon)
           if(dist[0] <= maxdistH):
              print "  --> Stations located within 100 m:",dist[0],stationlist[i].ele-stationlist[j].ele,stationlist[i].staorg,stationlist[i].ICC,stationlist[j].staorg,stationlist[j].ICC

    latdv = stationlist[i].lat
    latcv = "N"
    if(latdv<0):
       latcv = "S"
       latdv = latdv*(-1.0)

    londv = stationlist[i].lon
    loncv = "E"
    if(londv<0):
       loncv = "W"
       londv = londv*(-1.0)

    #staorg,staali,net,loc,lat,lon,ele
    if(outform == "VELEST"):
       fp1.write("%-4s%7.4f%1s %8.4f%1s %5d 1 %3d %5.2f %5.2f\n" % (stationlist[i].staali,latdv,latcv,londv,loncv,stationlist[i].ele*1000.0,stationlist[i].ICC,0.0,0.0))

    if(outform == "SIMULPS"):
       latd = int(latdv)
       latm = abs(latdv-latd)*60.0
       lond = int(londv)
       lonm = abs(londv-lond)*60.0
       if(stationlist[i].ele*1000.0 < -999):
          print "--> WARNING, elevation < -999 not possible for SIMULPS at the moment, need to be removed and NStat in header needs to be adjusted",stationlist[i].staali
       fp1.write("  %-4s%02d%1s%05.2f %03d%1s%05.2f%5d %4.2f %4.2f %2d\n" % (stationlist[i].staali,latd,latcv,latm,lond,loncv,lonm,stationlist[i].ele*1000.0,0.00,0.00,1))

    if(outform == 'NLL_Grid2Time'):
       #GTSRCE MABI LATLON 46.05490 10.51400 0.0 1.853
       fp1.write("GTSRCE %-5s LATLON %9.5f %10.5f 0.0 %6.3f\n" % (stationlist[i].staorg,stationlist[i].lat,stationlist[i].lon,stationlist[i].ele))

    if(outform == 'NLL_LocAlias'):
       #LOCALIAS code alias yearStart monthStart dayStart yearEnd monthEnd dayEnd
       #LOCALIAS BNAL BNALP 0 0 0 9999 99 99
       fp1.write("LOCALIAS %-8s %-8s 0 0 0 9999 99 99\n" % (stationlist[i].staali,stationlist[i].staorg))

    if(outform == 'NLL_StatDelay'):
       #LOCDELAY name phase n_readings delay
       #LOCDELAY BHE01 S 1  1.16
       #Write for P:
       fp1.write("LOCDELAY %-5s P  %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noP,stationlist[i].SCP))
       fp1.write("LOCDELAY %-5s Pg %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noP,stationlist[i].SCP))
       fp1.write("LOCDELAY %-5s Pn %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noP,stationlist[i].SCP))
       fp1.write("LOCDELAY %-5s P1 %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noP,stationlist[i].SCP))

       #Write for S:
       fp1.write("LOCDELAY %-5s S  %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noS,stationlist[i].SCS))
       fp1.write("LOCDELAY %-5s Sg %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noS,stationlist[i].SCS))
       fp1.write("LOCDELAY %-5s Sn %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noS,stationlist[i].SCS))
       fp1.write("LOCDELAY %-5s S1 %6d %7.3f\n" % (stationlist[i].staorg,stationlist[i].noS,stationlist[i].SCS))

    if(outform == 'hypoDD'):
       #BFO        48.330100     8.329600    589.0
       #BNI        45.052000     6.678000   1395.0
       fp1.write("%-9s %10.6f %12.6f %8.1f\n" % (stationlist[i].staorg,stationlist[i].lat,stationlist[i].lon,stationlist[i].ele*1000.0))

#VELEST need one empty line at the end:
if(outform == "VELEST"):
   fp1.write("\n")

#Close file:
fp1.close()

print ""
print "Extracted stations exported to     :",outfi
print "Duplicate station codes reported in:",logfile_DUPLICA

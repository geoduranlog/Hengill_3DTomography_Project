#!/usr/bin/env python2

import re
import math
import time
import numpy as np
import cPickle
from obspy.core import UTCDateTime
#from obspy.core.util import gps2DistAzimuth
from obspy.geodetics import gps2dist_azimuth
#import pg
import os
import shutil
import glob
import string
import copy
import optparse

global VERY_SMALL_DOUBLE, SMALL_DOUBLE
VERY_SMALL_DOUBLE = 1.0e-30
SMALL_DOUBLE = 1.0e-8


###########################################################################################
# Classes:
#------------------------------------------------------------------------------------------
class hypo(object):

      def __init__(self,yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,matchIDX):
          self.yy = int(yy)
          self.mm = int(mm)
          self.dd = int(dd)
          self.hh = int(hh)
          self.mi = int(mi)
          self.ss = float(ss)
          self.timestamp = float(timestamp)
          self.lat = float(lat)
          self.lon = float(lon)
          self.dep = float(dep)
          self.lonerr = float(lonerr)
          self.laterr = float(laterr)
          self.deperr = float(deperr)
          self.mag = float(mag)
          self.mtype = str(mtype)
          self.mnobs = int(mnobs)
          self.merr  = float(merr)
          self.mmeth = str(mmeth)
          self.rms = float(rms)
          self.gap = int(gap)
          self.mdist = float(mdist)
          self.nobs = int(nobs)
          self.etype = str(etype)
          self.lqual = str(lqual)
          self.chx   = float(chx)                 #CH-coordinates
          self.chy   = float(chy)                 #CH-coordinates
          self.agency = str(agency)
          self.evID1 = str(evID1)                 #Filename
          self.evID2 = long(evID2)                #DD-ID
          self.evPID = str(evPID)                 #publicID of event
          self.orID1 = str(orID1)                 #publicID of origin
          self.methodID = str(methodID)
          self.earthmodelID = str(earthmodelID)
          self.author = str(author)
          self.region = str(region)
          self.fstatus = str(fstatus)
          self.difflon = float(difflon)
          self.difflat = float(difflat)
          self.diffdep = float(diffdep)
          self.diffepi = float(diffepi)
          self.diffoti = float(diffoti)
          self.diffazi = float(diffazi)
          self.latalt = float(latalt)
          self.lonalt = float(lonalt)
          self.depalt = float(depalt)
          self.matchIDX = int(matchIDX)

      def __str__(self):
          return '%04d/%02d/%02d %02d:%02d:%06.3f %9.4f %8.4f %6.2f %4.1f %-5s %6.3f %3d %7.1f %3d %-2s %1s %3.0f %3.0f %-7s %-15s %-35s %-3s %15d %-15s %-s %-s' % (self.yy,self.mm,self.dd,self.hh,self.mi,self.ss,self.lon,self.lat,self.dep,self.mag,self.mtype,self.rms,self.gap,self.mdist,self.nobs,self.etype,self.lqual,self.chx,self.chy,self.agency,self.author,self.region,self.fstatus,self.evID2,self.methodID[0:15],self.evID1,self.orID1)

#------------------------------------------------------------------------------------------
class phase(object):
    
      def __init__(self,stat,ptype,qual,tt,timestamp,evID,dist):
          self.stat = str(stat)
          self.ptype = str(ptype)
          self.qual = int(qual)
          self.tt   = float(tt)
          self.timestamp = float(timestamp)
          self.evID = str(evID)
          self.dist = float(dist)
#------------------------------------------------------------------------------------------
class station(object):

      def __init__(self,stat,netw,lat,lon,ele,nobsP,nobsS,nobsPw,nobsSw):
          self.stat = str(stat)
          self.netw = str(netw)
          self.lat = float(lat)
          self.lon = float(lon)
          self.ele = float(ele)
          self.nobsP = int(nobsP)
          self.nobsS = int(nobsS)
          self.nobsPw = float(nobsPw)
          self.nobsSw = float(nobsSw)
#------------------------------------------------------------------------------------------
###########################################################################################
#Subroutines/functions:

def writeCNVPhases(fp,Phases):
    
    nop = len(Phases)

    in6 = int(nop/6)
    rem = nop - (in6*6)

    #stat,ptype,qual,tt,timestamp,evID

    cnt = 0

    #Print full lines:
    for i in range(in6):
        fp.write("%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt,Phases[cnt+3].stat,Phases[cnt+3].ptype,Phases[cnt+3].qual,Phases[cnt+3].tt,Phases[cnt+4].stat,Phases[cnt+4].ptype,Phases[cnt+4].qual,Phases[cnt+4].tt,Phases[cnt+5].stat,Phases[cnt+5].ptype,Phases[cnt+5].qual,Phases[cnt+5].tt))
        cnt = cnt + 6

    #Print rest:
    if(rem == 1):
       fp.write("%-4s%-1s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt))
       cnt = cnt + 1

    if(rem == 2):
       fp.write("%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt))
       cnt = cnt + 2       

    if(rem == 3):
       fp.write("%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt))
       cnt = cnt + 3


    if(rem == 4):
       fp.write("%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt,Phases[cnt+3].stat,Phases[cnt+3].ptype,Phases[cnt+3].qual,Phases[cnt+3].tt))
       cnt = cnt + 4

    if(rem == 5):
       fp.write("%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f%-4s%-1s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt,Phases[cnt+3].stat,Phases[cnt+3].ptype,Phases[cnt+3].qual,Phases[cnt+3].tt,Phases[cnt+4].stat,Phases[cnt+4].ptype,Phases[cnt+4].qual,Phases[cnt+4].tt))
       cnt = cnt + 5   

    #End of event:
    fp.write("\n")

    #Checksum:
    if(nop != cnt):
       print "CHECKSUM ERROR for event",Phases[0].evID
       sys.exit()
    #else:
    #   print "Checksum correct"

    return
#------------------------------------------------------------------------------------------
def writeCNVHeader(fp,Hypo):

    #Convert year:
    if(Hypo.yy >= 2000):
       vyy = Hypo.yy - 2000
    else:
       vyy = Hypo.yy - 1900

    #Convert coordinates:
    latc = 'N'
    lat  = Hypo.lat
    if(lat < 0):
       lat = lat*-1.0
       latc = 'S'
    lonc = 'E'
    lon  = Hypo.lon
    if(lon < 0):
       lon = lon*-1.0
       lonc = 'W'

    fp.write("%02d%02d%02d %02d%02d %05.2f %7.4f%1s %8.4f%1s %6.2f %6.2f %6d %9.2f  EVID: %-s\n" % (vyy,Hypo.mm,Hypo.dd,Hypo.hh,Hypo.mi,Hypo.ss,lat,latc,lon,lonc,Hypo.dep,Hypo.mag,Hypo.gap,Hypo.rms,Hypo.evID1))

    return
#------------------------------------------------------------------------------------------
def getstatidx(stations,stat2find):

    idx = -1

    for n in range(len(stations)):
        if(stations[n].stat == stat2find):
           idx = n
           break

    return idx
#------------------------------------------------------------------------------------------
def loadVelestStations(ifile):

    statlist = []

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
        netw = ''
        stat = line[0:4]
        lat  = float(line[4:11])
        latc = line[11:12]
        lon  = float(line[13:21])
        lonc = line[21:22]
        elev = float(line[23:28])/1000.0

        if(latc == 'S'):
           lat = lat*-1.0

        if(lonc == 'W'):
           lon = lon*-1.0

        statlist.append(station(stat,netw,lat,lon,elev,0,0,0.0,0.0))

    fp.close()

    return statlist
#------------------------------------------------------------------------------------------
def loadCNV(ifile,ptype,dbgmode):

    evlist = []
    phlist = []

    #Open
    fp = open(ifile,"r")
    ID1 = ''

    for line in fp:

        #Split string to get event ID        
        s = string.split(line,"EVID:",-1)

        #Test if line is empty:
        if(line[69:74]=="EVID:"):
           yys=int(line[0:2])
           if(yys <= 60):
              yy=2000+yys
           else:
              yy=1900+yys

           mm=int(line[2:4])
           dd=int(line[4:6])
           hh=int(line[7:9])
           mi=int(line[9:11])
           ss=float(line[12:17])
           ID1=string.strip(s[1])
           lat=float(line[18:25])
           latc=line[25:26]
           if(latc=="S"):
              lat=lat*-1
           lon=float(line[27:35])
           lonc=line[35:36]
           if(lonc=="W"):
              lon=lon*-1
           dep=float(line[37:43])
           mag=float(line[44:50])
           gap=int(line[54:57])
           rms=float(line[62:67])

           #Get swiss coordinates:
           chx,chy = celleb(lon,lat)
           if(chy < 62.0) or (chy > 302.0) or (chx < 480.0) or (chx > 847.5):
              chy = 999
              chx = 999

           #Unknown:
           mat    = '?'
           loc    = '?'
           mdi=-9
           nob=-9
           ety='?'
           lqa='?'
           ID2 = long(999999999)
           agy    = "SED_VE"
           region = "None"
           author = "VE"
           ID3    = "None"          #Origin publicID (string)
           ID4    = "None"          #Event  publicID (string)
           fstat  = "VE"
           laterr = -9
           lonerr = -9
           deperr = -9
           mnobs  = -9
           merr   = -9
           modlID = '?'
           magmeth = '?'
           difflon = -9
           difflat = -9
           diffdep = -9
           diffepi = -9
           diffoti = -9
           diffazi = -9
           latalt = -9
           lonalt = -9
           depalt = -9

           #Check second:
           if(ss == 60.0):
              print 'Second == 60.0 --> substract 0.00001'
              ss = ss - 0.00001

           wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.999995 -> 60.0 Fix: increase precision

           evlist.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,-9))
           #if(dbgmode == 0):
              #print hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt)
           
           #next line:
           continue

        else:
           #Check if line containes one phase:
           if(len(ID1)>0) and (len(line.rstrip())==12):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),float(line[6:12]),wdate.timestamp+float(line[6:12]),ID1,-9))                       

           #Check if line containes two phase:
           if(len(ID1)>0) and (len(line.rstrip())==24):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),float(line[6:12]),wdate.timestamp+float(line[6:12]),ID1,-9))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),float(line[18:24]),wdate.timestamp+float(line[18:24]),ID1,-9))

           #Check if line containes three phase:
           if(len(ID1)>0) and (len(line.rstrip())==36):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),float(line[6:12]),wdate.timestamp+float(line[6:12]),ID1,-9))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),float(line[18:24]),wdate.timestamp+float(line[18:24]),ID1,-9))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),float(line[30:36]),wdate.timestamp+float(line[30:36]),ID1,-9))

           #Check if line containes four phase:
           if(len(ID1)>0) and (len(line.rstrip())==48):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),float(line[6:12]),wdate.timestamp+float(line[6:12]),ID1,-9))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),float(line[18:24]),wdate.timestamp+float(line[18:24]),ID1,-9))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),float(line[30:36]),wdate.timestamp+float(line[30:36]),ID1,-9))
              phlist.append(phase(line[36:40],line[40:41],int(line[41:42]),float(line[42:48]),wdate.timestamp+float(line[42:48]),ID1,-9))
 
           #Check if line containes five phase:
           if(len(ID1)>0) and (len(line.rstrip())==60):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),float(line[6:12]),wdate.timestamp+float(line[6:12]),ID1,-9))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),float(line[18:24]),wdate.timestamp+float(line[18:24]),ID1,-9))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),float(line[30:36]),wdate.timestamp+float(line[30:36]),ID1,-9))
              phlist.append(phase(line[36:40],line[40:41],int(line[41:42]),float(line[42:48]),wdate.timestamp+float(line[42:48]),ID1,-9))
              phlist.append(phase(line[48:52],line[52:53],int(line[53:54]),float(line[54:60]),wdate.timestamp+float(line[54:60]),ID1,-9))

           #Check if line containes five phase:
           if(len(ID1)>0) and (len(line.rstrip())==72):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),float(line[6:12]),wdate.timestamp+float(line[6:12]),ID1,-9))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),float(line[18:24]),wdate.timestamp+float(line[18:24]),ID1,-9))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),float(line[30:36]),wdate.timestamp+float(line[30:36]),ID1,-9))
              phlist.append(phase(line[36:40],line[40:41],int(line[41:42]),float(line[42:48]),wdate.timestamp+float(line[42:48]),ID1,-9))
              phlist.append(phase(line[48:52],line[52:53],int(line[53:54]),float(line[54:60]),wdate.timestamp+float(line[54:60]),ID1,-9))
              phlist.append(phase(line[60:64],line[64:65],int(line[65:66]),float(line[66:72]),wdate.timestamp+float(line[66:72]),ID1,-9))

    #Close
    fp.close()

    return evlist,phlist

#------------------------------------------------------------------------------------------
def loadrefevents(ifile,dbgmode):

    list = []

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<10):
           print "Line not considered:", line
        else:
           #Test if line is comment:
           #t=line.split()
           t = string.split(line,None,-1)
           s = string.split(line,"|",-1)

           if(t[0].startswith("#") == False):
              lon=float(line[0:8])
              lat=float(line[9:16])
              dep=float(line[17:22])
              mag=float(line[23:27])
              mat=string.strip(line[28:30])
              yy=int(line[31:35])
              mm=int(line[36:38])
              dd=int(line[39:41])
              hh=int(line[42:44])
              mi=int(line[45:47])
              ss=float(line[48:53])
              loc=string.strip(line[54:58])
              rms=float(line[60:65])
              gap=int(line[66:69])
              mdi=float(line[70:75])
              nob=int(line[76:79])
              ety=string.strip(line[80:81])
              lqa=string.strip(line[82:83])
              ID1=string.strip(line[84:104])    #File-ID (string)
              ID2=long(line[105:114])           #DD-evID (long)     

              #Get swiss coordinates:
              chx,chy = celleb(lon,lat)
              if(chy < 62.0) or (chy > 302.0) or (chx < 480.0) or (chx > 847.5):
                chy = 999
                chx = 999

              #Check if event was extracted from KP file structure:
              #if(ID1[0:2]=="KP"):
              if(loc != "SEDS"):
                 agy    = "SED_KP"
                 region = "None"
                 author = "KP"
                 ID3    = "None"          #Origin publicID (string)
                 ID4    = "None"          #Event  publicID (string)
                 fstat  = "REF"
                 laterr = -9
                 lonerr = -9
                 deperr = -9
                 mnobs  = -9
                 merr   = -9
                 modlID = '?'
                 magmeth = '?'
                 difflon = -9
                 difflat = -9
                 diffdep = -9
                 diffepi = -9
                 diffoti = -9
                 diffazi = -9

              if(loc == "SEDS"):
                 ID3    = string.strip(s[5])
                 ID4    = string.strip(s[6])
                 agy    = "SED"
                 region = string.strip(s[4])
                 author = string.strip(s[3])
                 fstat  = "REF"
                 laterr = -9
                 lonerr = -9
                 deperr = -9
                 mnobs  = -9
                 merr   = -9
                 modlID = '?'
                 magmeth = '?'
                 difflon = -9
                 difflat = -9
                 diffdep = -9
                 diffepi = -9
                 diffoti = -9
                 diffazi = -9
                 latalt = -9
                 lonalt = -9
                 depalt = -9


              #Check if event was extracted from SC3 DB structure:
              #STILL MISSING 
              #else:

              #Check if format is correct:
              if(ss < 60.0) and (mi < 60) and (hh < 24) and (dd < 32) and (mm < 13):
                 wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.99995 -> 60.0 Fix: increase precision
                 #wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss))

                 #print ID1
                 list.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,-9))
                 if(dbgmode == 1):
                    print hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,-9)
              else:
                 print "Warning: Corrupt origin time ",datetime2(yy,mm,dd,hh,mi,ss),ID1
           else:
              print "Line not considered:", line

    #Close
    fp.close()

    return list
#------------------------------------------------------------------------------------------
def celleb(lon,lat):

    #**********************
    #  lon
    #  SCHWEIZ. PROJEKTIONSSYSTEM  FORMELN VON H. ODERMATT
    #  TRANSFORMATION ELLIPSOID - EBENE
    #  L,B  laenge und breite in grad
    #  Y,X LANDESKOORDINATEN IN KILO-METER y= e-w; x=n-s
    #  MY  MERIDIANKONVERGENZ ( SEXAG. SEK.)
    #
    #  from nlloc/eth_custom/get_region_name_nr.c
    #*************************

    #To be consistent with orignal code:
    l = lon      #Laenge = longitude
    b = lat      #Breite = latitude

    #Initials:
    y = 999
    x = 999

    #Define TOP 8
    TOP = 8

    bb = 169028.66
    bl = 26782.5

    a = l
    a = a*3600.0 - bl
    c = b
    c = c * 3600.0 - bb

    d = []
    e = []
    f = []
    rw = []
    iw = []

    d.append(float( 0.0))                 #d00
    d.append(float( 0.68382546262761))    #d01
    d.append(float(-3.91798328045E-8))    #d02
    d.append(float( 1.4965410352E-15))    #d03
    d.append(float(-8.039471422E-23))     #d04
    d.append(float( 7.0021390E-30))       #d05
    d.append(float(-5.586904E-37))        #d06
    d.append(float( 4.0402E-44))          #d07
    d.append(float(-3.06E-51))            #d08

    e.append(float( 0.0))                 #e00
    e.append(float( 2.3635916074715E-2))  #e01
    e.append(float( 0.0))                 #e02
    e.append(float( 4.527219881E-17))     #e03
    e.append(float(-3.89081120E-24))      #e04
    e.append(float( 2.3407700E-31))       #e05
    e.append(float(-1.59674E-38))         #e06
    e.append(float( 1.287E-45))           #e07
    e.append(float( 0.0))                 #e08

    f.append(float( 0.0))                 #f00
    f.append(float( 4.515344386039E1))    #f01
    f.append(float( 1.17912305209E-4))    #f02
    f.append(float( 5.8474201864E-10))    #f03
    f.append(float( 2.73386187E-15))      #f04
    f.append(float( 1.4308547E-20))       #f05
    f.append(float( 7.66562E-26))         #f06
    f.append(float( 4.2445E-31))          #f07
    f.append(float( 2.40E-36))            #f08

    rw.append(float( 0.0))
    iw.append(float( 0.0))

    p = 30.91849390613 * a
    q = c * f[8]

    for i in range(TOP-1,0,-1):
        q = c * (q + f[i])

    rw.append(float(q)) #rw[1]
    iw.append(float(p)) #iw[1]

    for i in range(2,TOP+1,+1):
        rw.append(float(q*rw[i-1]-p*iw[i-1]))
        iw.append(float(p*rw[i-1]+q*iw[i-1]))

    dx = d[TOP]*rw[TOP]
    dy = d[TOP]*iw[TOP]
    my = 0.0

    for i in range(TOP-1,0,-1):
        dx=dx+d[i]*rw[i]
        dy=dy+d[i]*iw[i]
        my=my+e[i]*iw[i]

    dx=dx+200000.0
    dy=dy+600000.0
    dy= dy/1000.0
    dx= dx/1000.0
    y =  dy
    x =  dx

    #print "y/x:",y,"/",x 
    return y,x
#------------------------------------------------------------------------------------------
def export2kml(ofile,hyplist):

    lS=len(hyplist)

    #Open file:
    fpkml = open(ofile,"w")

    #Generate header of kml:
    fpkml.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    fpkml.write("<kml xmlns=\"http://earth.google.com/kml/2.0\"> <Document>\n")

    #01) Write Earthquakes >= 3.5:
    #####################################################################################

    fpkml.write("<Folder>\n")
    fpkml.write("   <name>SED_Earthquakes Ml ge 3.5</name>\n")
    fpkml.write("   <open>1</open>\n")
    fpkml.write("   <Style id=\"large_dot_red\">\n")
    fpkml.write("        <IconStyle>\n")
    fpkml.write("            <color>ff0000ff</color>\n")
    fpkml.write("            <colorMode>normal</colorMode>>\n")
    fpkml.write("            <scale>2.00</scale>               <!-- float -->\n")
    fpkml.write("             <heading>0</heading>             <!-- float -->\n")
    fpkml.write("             <Icon>\n")
    #http://maps.google.com/mapfiles/kml/pal2/icon18.png
    #fpkml.write("                <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n")
    fpkml.write("                 <href>http://maps.google.com/mapfiles/kml/shapes/dot.png</href>\n")
    fpkml.write("             </Icon>\n")
    fpkml.write("             <hotSpot x=\"0.5\"  y=\"0.5\"\n")
    fpkml.write("              xunits=\"fraction\" yunits=\"fraction\"/>    <!-- kml:vec2Type -->>\n")
    fpkml.write("        </IconStyle>\n")
    fpkml.write("        <LabelStyle>\n")
    fpkml.write("             <scale>0.75</scale>\n")
    fpkml.write("             <color>#ff0000ff</color>\n")
    fpkml.write("        </LabelStyle>\n")
    fpkml.write("        <BalloonStyle>\n")
    fpkml.write("             <!-- a background color for the balloon -->\n")
    fpkml.write("             <bgColor>#FFFFFFFF</bgColor>\n")
    fpkml.write("             <!-- styling of the balloon text -->\n")
    fpkml.write("             <text><![CDATA[\n")
    fpkml.write("             <b><font color=\"#CC0000\" size=\"+2\">$[name]</font></b>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <font face=\"Arial\">$[description]</font>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             SED Double-Difference Location: Earthquake\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <!-- insert the to/from hyperlinks -->\n")
    fpkml.write("             $[geDirections]\n")
    fpkml.write("             ]]></text>\n")
    fpkml.write("        </BalloonStyle>\n")
    fpkml.write("   </Style>\n")

    #yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus

    for i in range(lS):
        if((hyplist[i].fstatus == "NEW") or (hyplist[i].fstatus == "UPD") or (hyplist[i].fstatus == "REF")) and (hyplist[i].etype == "T") and (hyplist[i].mag >= 3.5):
           fpkml.write("<Placemark><styleUrl>#large_dot_red</styleUrl>\n")
           fpkml.write(" <name>%09d / %-s</name>\n" % (hyplist[i].evID2,hyplist[i].evID1))
           fpkml.write(" <description>\n")
           fpkml.write("   %04d/%02d/%02d %02d:%02d:%04.1f <br/> Depth: %5.1f km; Magnitude: %4.1f %-3s Type: %s <br/> %-s\n" % (hyplist[i].yy,hyplist[i].mm,hyplist[i].dd,hyplist[i].hh,hyplist[i].mi,hyplist[i].ss,hyplist[i].dep,hyplist[i].mag,hyplist[i].mtype,hyplist[i].etype,hyplist[i].author))
           fpkml.write(" </description>\n")
           fpkml.write("          <Point><coordinates>%010.4f,%010.4f</coordinates></Point>\n" % (hyplist[i].lon,hyplist[i].lat))
           fpkml.write("</Placemark>\n")

    fpkml.write("</Folder>\n")


    #02) Write Earthquakes 2.5 <= Ml < 3.5:
    #####################################################################################

    fpkml.write("<Folder>\n")
    fpkml.write("   <name>SED_Earthquakes 2.5 le Ml lt 3.5</name>\n")
    fpkml.write("   <open>1</open>\n")
    fpkml.write("   <Style id=\"medium_dot_red\">\n")
    fpkml.write("        <IconStyle>\n")
    fpkml.write("            <color>ff0000ff</color>\n")
    fpkml.write("            <colorMode>normal</colorMode>>\n")
    fpkml.write("            <scale>1.50</scale>               <!-- float -->\n")
    fpkml.write("             <heading>0</heading>             <!-- float -->\n")
    fpkml.write("             <Icon>\n")
    #http://maps.google.com/mapfiles/kml/pal2/icon18.png
    #fpkml.write("                <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n")
    fpkml.write("                 <href>http://maps.google.com/mapfiles/kml/shapes/dot.png</href>\n")
    fpkml.write("             </Icon>\n")
    fpkml.write("             <hotSpot x=\"0.5\"  y=\"0.5\"\n")
    fpkml.write("              xunits=\"fraction\" yunits=\"fraction\"/>    <!-- kml:vec2Type -->>\n")
    fpkml.write("        </IconStyle>\n")
    fpkml.write("        <LabelStyle>\n")
    fpkml.write("             <scale>0.75</scale>\n")
    fpkml.write("             <color>#ff0000ff</color>\n")
    fpkml.write("        </LabelStyle>\n")
    fpkml.write("        <BalloonStyle>\n")
    fpkml.write("             <!-- a background color for the balloon -->\n")
    fpkml.write("             <bgColor>#FFFFFFFF</bgColor>\n")
    fpkml.write("             <!-- styling of the balloon text -->\n")
    fpkml.write("             <text><![CDATA[\n")
    fpkml.write("             <b><font color=\"#CC0000\" size=\"+2\">$[name]</font></b>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <font face=\"Arial\">$[description]</font>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             SED Double-Difference Location: Earthquake\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <!-- insert the to/from hyperlinks -->\n")
    fpkml.write("             $[geDirections]\n")
    fpkml.write("             ]]></text>\n")
    fpkml.write("        </BalloonStyle>\n")
    fpkml.write("   </Style>\n")

    #yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus

    for i in range(lS):
        if((hyplist[i].fstatus == "NEW") or (hyplist[i].fstatus == "UPD") or (hyplist[i].fstatus == "REF")) and (hyplist[i].etype == "T") and (hyplist[i].mag >= 2.5) and (hyplist[i].mag < 3.5):
           fpkml.write("<Placemark><styleUrl>#medium_dot_red</styleUrl>\n")
           fpkml.write("     <name>%09d / %-s</name>\n" % (hyplist[i].evID2,hyplist[i].evID1))
           fpkml.write("     <description>\n")
           fpkml.write("        %04d/%02d/%02d %02d:%02d:%04.1f <br/> Depth: %5.1f km; Magnitude: %4.1f %-3s Type: %s <br/> %-s\n" % (hyplist[i].yy,hyplist[i].mm,hyplist[i].dd,hyplist[i].hh,hyplist[i].mi,hyplist[i].ss,hyplist[i].dep,hyplist[i].mag,hyplist[i].mtype,hyplist[i].etype,hyplist[i].author))
           fpkml.write("     </description>\n")
           fpkml.write("     <Point><coordinates>%010.4f,%010.4f</coordinates></Point>\n" % (hyplist[i].lon,hyplist[i].lat))
           fpkml.write("</Placemark>\n")

    fpkml.write("</Folder>\n")

    #03) Write Earthquakes Ml < 2.5:
    #####################################################################################

    fpkml.write("<Folder>\n")
    fpkml.write("   <name>SED_Earthquakes Ml lt 2.5</name>\n")
    fpkml.write("   <open>1</open>\n")
    fpkml.write("   <Style id=\"small_dot_red\">\n")
    fpkml.write("        <IconStyle>\n")
    fpkml.write("            <color>ff0000ff</color>\n")
    fpkml.write("            <colorMode>normal</colorMode>>\n")
    fpkml.write("            <scale>0.75</scale>               <!-- float -->\n")
    fpkml.write("             <heading>0</heading>             <!-- float -->\n")
    fpkml.write("             <Icon>\n")
    #http://maps.google.com/mapfiles/kml/pal2/icon18.png
    #fpkml.write("                <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n")
    fpkml.write("                 <href>http://maps.google.com/mapfiles/kml/shapes/dot.png</href>\n")
    fpkml.write("             </Icon>\n")
    fpkml.write("             <hotSpot x=\"0.5\"  y=\"0.5\"\n")
    fpkml.write("              xunits=\"fraction\" yunits=\"fraction\"/>    <!-- kml:vec2Type -->>\n")
    fpkml.write("        </IconStyle>\n")
    fpkml.write("        <LabelStyle>\n")
    fpkml.write("             <scale>0.00</scale>\n")
    fpkml.write("             <color>#ff0000ff</color>\n")
    fpkml.write("        </LabelStyle>\n")
    fpkml.write("        <BalloonStyle>\n")
    fpkml.write("             <!-- a background color for the balloon -->\n")
    fpkml.write("             <bgColor>#FFFFFFFF</bgColor>\n")
    fpkml.write("             <!-- styling of the balloon text -->\n")
    fpkml.write("             <text><![CDATA[\n")
    fpkml.write("             <b><font color=\"#CC0000\" size=\"+2\">$[name]</font></b>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <font face=\"Arial\">$[description]</font>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             SED Double-Difference Location: Earthquake\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <!-- insert the to/from hyperlinks -->\n")
    fpkml.write("             $[geDirections]\n")
    fpkml.write("             ]]></text>\n")
    fpkml.write("        </BalloonStyle>\n")
    fpkml.write("   </Style>\n")

    #yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus

    for i in range(lS):
        if((hyplist[i].fstatus == "NEW") or (hyplist[i].fstatus == "UPD") or (hyplist[i].fstatus == "REF")) and (hyplist[i].etype == "T") and (hyplist[i].mag < 2.5):
           fpkml.write("<Placemark><styleUrl>#small_dot_red</styleUrl>\n")
           fpkml.write("     <name>%09d / %-s</name>\n" % (hyplist[i].evID2,hyplist[i].evID1))
           fpkml.write("     <description>\n")
           fpkml.write("        %04d/%02d/%02d %02d:%02d:%04.1f <br/> Depth: %5.1f km; Magnitude: %4.1f %-3s Type: %s <br/> %-s\n" % (hyplist[i].yy,hyplist[i].mm,hyplist[i].dd,hyplist[i].hh,hyplist[i].mi,hyplist[i].ss,hyplist[i].dep,hyplist[i].mag,hyplist[i].mtype,hyplist[i].etype,hyplist[i].author))
           fpkml.write("     </description>\n")
           fpkml.write("     <Point><coordinates>%010.4f,%010.4f</coordinates></Point>\n" % (hyplist[i].lon,hyplist[i].lat))
           fpkml.write("</Placemark>\n")

    fpkml.write("</Folder>\n")

    #04) Write  all blasts:
    #####################################################################################
    fpkml.write("<Folder>\n")
    fpkml.write("   <name>SED_Blasts</name>\n")
    fpkml.write("   <open>1</open>\n")
    fpkml.write("   <Style id=\"small_dot_yellow\">\n")
    fpkml.write("        <IconStyle>\n")
    fpkml.write("            <color>ffffff00</color>\n")
    fpkml.write("            <colorMode>normal</colorMode>>\n")
    fpkml.write("            <scale>0.75</scale>               <!-- float -->\n")
    fpkml.write("             <heading>0</heading>             <!-- float -->\n")
    fpkml.write("             <Icon>\n")
    #http://maps.google.com/mapfiles/kml/pal2/icon18.png
    #fpkml.write("                <href>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</href>\n")
    fpkml.write("                 <href>http://maps.google.com/mapfiles/kml/shapes/dot.png</href>\n")
    fpkml.write("             </Icon>\n")
    fpkml.write("             <hotSpot x=\"0.5\"  y=\"0.5\"\n")
    fpkml.write("              xunits=\"fraction\" yunits=\"fraction\"/>    <!-- kml:vec2Type -->>\n")
    fpkml.write("        </IconStyle>\n")
    fpkml.write("        <LabelStyle>\n")
    fpkml.write("             <scale>0.00</scale>\n")
    fpkml.write("             <color>ffffff00</color>\n")
    fpkml.write("        </LabelStyle>\n")
    fpkml.write("        <BalloonStyle>\n")
    fpkml.write("             <!-- a background color for the balloon -->\n")
    fpkml.write("             <bgColor>#FFFFFFFF</bgColor>\n")
    fpkml.write("             <!-- styling of the balloon text -->\n")
    fpkml.write("             <text><![CDATA[\n")
    fpkml.write("             <b><font color=\"#CC0000\" size=\"+2\">$[name]</font></b>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <font face=\"Arial\">$[description]</font>\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             SED Double-Difference Location: Blast\n")
    fpkml.write("             <br/><br/>\n")
    fpkml.write("             <!-- insert the to/from hyperlinks -->\n")
    fpkml.write("             $[geDirections]\n")
    fpkml.write("             ]]></text>\n")
    fpkml.write("        </BalloonStyle>\n")
    fpkml.write("   </Style>\n")

    #yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus

    for i in range(lS):
        if((hyplist[i].fstatus == "NEW") or (hyplist[i].fstatus == "UPD") or (hyplist[i].fstatus == "REF")) and (hyplist[i].etype == "E"):
           fpkml.write("<Placemark><styleUrl>#small_dot_yellow</styleUrl>\n")
           fpkml.write("     <name>%09d / %-s</name>\n" % (hyplist[i].evID2,hyplist[i].evID1))
           fpkml.write("     <description>\n")
           fpkml.write("        %04d/%02d/%02d %02d:%02d:%04.1f <br/> Depth: %5.1f km; Magnitude: %4.1f %-3s Type: %s <br/> %-s\n" % (hyplist[i].yy,hyplist[i].mm,hyplist[i].dd,hyplist[i].hh,hyplist[i].mi,hyplist[i].ss,hyplist[i].dep,hyplist[i].mag,hyplist[i].mtype,hyplist[i].etype,hyplist[i].author))
           fpkml.write("     </description>\n")
           fpkml.write("     <Point><coordinates>%010.4f,%010.4f</coordinates></Point>\n" % (hyplist[i].lon,hyplist[i].lat))
           fpkml.write("</Placemark>\n")

    fpkml.write("</Folder>\n")

    #Close kml file:
    fpkml.write("</Document> </kml>\n")
    fpkml.close()

    print ""
    print "Events have been exported to kml format:",ofile

    return
#------------------------------------------------------------------------------------------

###########################################################################################
# Main code:
###########################################################################################

#Load two velest-cnv files. Compare evens based on KP
#event ID and write matching events to output (with both hypocenters + shift)

#Now from command line:
##################################################################################################################
#Get command line arguments:
oparser = optparse.OptionParser()
oparser.add_option('-P', '--PphaseFile',  action='store', dest='PphaseFile' ,  help='CNV file with P phases labeled as P')
oparser.add_option('-S', '--SphaseFile',  action='store', dest='SphaseFile',   help='CNV file with S phases labeled as P')
oparser.add_option('-H', '--HypoInfo',    action='store', dest='HypoInfo',     help='Defines source of hypocenter information; P: From PphaseFile (Default); S: From SphaseFile')
oparser.add_option('-o', '--output',      action='store', dest='OutputFile',   help='Filename of merged P+S output') 
oparser.add_option('-s', '--statfile',    action='store', dest='StationFile',  help='Velest station file')
oparser.add_option('-m', '--mergemode',    action='store', dest='mergemode',  help='Merge-Mode: MatchOnly (Default); ALL (Merge all events no matter if in both files); ALL+ (same as ALL, but events are complemented by phases in complementary file provided by -c option)')
oparser.add_option('-c', '--complementaryPhaseFile', action='store', dest='CompPhaseFile',  help='CNV File with complementary phases to be added in case -m ALL+ option is used (Default: None), additional S waves are added tp P-only events and vice verca.')
(options, arg) = oparser.parse_args()

#Check if mergemode is included in list of input arguments:
if(options.mergemode != None):
   mergemode = options.mergemode
else:
   mergemode = 'MatchOnly'
if(mergemode != 'MatchOnly') and (mergemode != 'ALL') and (mergemode != 'ALL+'):
   print "--> WARNING: mergemode option unknown ",mergemode," -> set mergemode to MatchOnly"
   mergemode = 'MatchOnly'

#Check if CompPhaseFile is included in list of input arguments:
if(options.CompPhaseFile != None):
   CompPhaseFile = options.CompPhaseFile
else:
   CompPhaseFile = 'None'

#Check if PphaseFile is included in list of input arguments:
if(options.PphaseFile != None):
   PphaseFile = options.PphaseFile
else:
   PphaseFile = '/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/P+S_Inversion/MergeModels/CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2.cnv'

#Check if SphaseFile is included in list of input arguments:
if(options.SphaseFile != None):
   SphaseFile = options.SphaseFile
else:
   SphaseFile = '/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/P+S_Inversion/MergeModels/CS_VMGralTrendMiddle_invertratio1_S_ratio_rm.cnv'

#Check if OutputFile is included in list of input arguments:
if(options.OutputFile != None):
   OutputFile = options.OutputFile
else:
   OutputFile = 'cnv_merge_P2S_Hengill.out.cnv'

#Check if StationFile is included in list of input arguments:
if(options.StationFile != None):
   StationFile = options.StationFile
else:
   StationFile = '/Users/alejandro/ICELAND_Project/Min1D_Model/VELEST/TESTS/P+S_Inversion/MergeModels/CS_VMGralTrendMiddle_invertratio1_P_CNVFeb2020_2.sta'

#Check if HypoInfo is included in list of input arguments:
if(options.HypoInfo != None):
   HypoInfo = options.HypoInfo
else:
   HypoInfo = 'P'

#Hardcoded-log file
LogFile1 = 'cnv_merge_PP2PS_Phases+Stations.log'
LogFile2 = 'cnv_merge_PP2PS_EventMerging.log'
##################################################################################################################

#Debug-Mode (0 or 1):
dbgmode = 0

#Load station file:
print ""
print "--> Load station file from:",StationFile
stations = loadVelestStations(StationFile)
print '    Number stations:',len(stations)

#Load P-phases:
print ""
print "--> Load P-phase CNV from: ",PphaseFile
Phypos,Pphases = loadCNV(PphaseFile,'P',dbgmode)
print '    Number events P:',len(Phypos)
print '    Number phases P:',len(Pphases)
print ""

#Load S-phases:
print "--> Load S-phase CNV from: ",SphaseFile
Shypos,Sphases = loadCNV(SphaseFile,'S',dbgmode)
print '    Number events S:',len(Shypos)
print '    Number phases S:',len(Sphases)
print ""

#Load complementary P+S file if requested:
if(CompPhaseFile != "None") and (os.path.isfile(CompPhaseFile)) and (mergemode == 'ALL+'):
   print "--> Load complementary P+S phases from CNV file: ",CompPhaseFile
   Chypos,Cphases = loadCNV(CompPhaseFile,'PS',dbgmode)
else:
   Chypos=[]
   Cphases=[]
print '    Number events C:',len(Chypos)
print '    Number phases C:',len(Cphases)
print ""

#Make sure S phases are labeled as "S" in merged file:
for j in range(len(Sphases)):
    Sphases[j].ptype = 'S' 

#Calculate some statistics:
error = np.zeros([3,5])

#Define errors - SED:
#error[0,0] = 0.025
#error[0,1] = 0.050
#error[0,2] = 0.100
#error[0,3] = 0.200
#error[0,4] = 0.400


#Define errors - COSEISMIQ:
error[0,0] = 0.0125
error[0,1] = 0.025
error[0,2] = 0.050
error[0,3] = 0.100
error[0,4] = 0.200



fp1 = open(LogFile1,'w')
fp2 = open(LogFile2,'w')

#Loop over P-Catalog to get statistics on error, station observations, distance...
for i in range(len(Phypos)):
    for j in range(len(Pphases)):
        if(Phypos[i].evID1 == Pphases[j].evID):

           k = getstatidx(stations,Pphases[j].stat)

           if(k>=0):
              dist = gps2dist_azimuth(Phypos[i].lat,Phypos[i].lon,stations[k].lat,stations[k].lon)
              fp1.write("%-5s %-8s %1d %5.3f %9.5f %10.5f %10.5f %5.1f\n" % (Pphases[j].stat,Pphases[j].ptype,Pphases[j].qual,1.0/(2.0**Pphases[j].qual),stations[k].lat,stations[k].lon,dist[0]/1000.0,dist[1]))
              error[1,Pphases[j].qual] += 1.0
              stations[k].nobsP += 1              
              stations[k].nobsPw = stations[k].nobsPw + (1/(2**Pphases[j].qual))
              Pphases[j].dist = dist[0]/1000.0

#Loop over S-Catalog to get statistics on error, station observations, distance...
for i in range(len(Shypos)):
    for j in range(len(Sphases)):
        if(Shypos[i].evID1 == Sphases[j].evID):

           k = getstatidx(stations,Sphases[j].stat)

           if(k>=0):
              dist = gps2dist_azimuth(Shypos[i].lat,Shypos[i].lon,stations[k].lat,stations[k].lon)
              fp1.write("%-5s %-8s %1d %5.3f %9.5f %10.5f %10.5f %5.1f\n" % (Sphases[j].stat,Sphases[j].ptype,Sphases[j].qual,1.0/(2.0**Sphases[j].qual),stations[k].lat,stations[k].lon,dist[0]/1000.0,dist[1]))
              error[2,Sphases[j].qual] += 1.0
              stations[k].nobsS += 1
              stations[k].nobsSw = stations[k].nobsSw + (1/(2**Sphases[j].qual))
              Sphases[j].dist = dist[0]/1000.0

for i in range(len(stations)):
    fp1.write("#STATION-STAT: %-5s %9.5f %10.5f %6d %6d %9.3f %9.3f\n" % (stations[i].stat,stations[i].lat,stations[i].lon,stations[i].nobsP,stations[i].nobsS,stations[i].nobsPw,stations[i].nobsSw))

#print error.shape[0]
#print error.shape[1]

#Calculate error statistsics:
errsumP = 0
errsumS = 0
averrorP = 0.0
averrorS = 0.0
for i in range(error.shape[1]):
    print "Error class",i,":",error[0,i],'Nobs P',int(error[1,i]),'Nobs S',int(error[2,i])
    errsumP = errsumP + int(error[1,i])
    errsumS = errsumS + int(error[2,i])
    averrorP = averrorP + (error[1,i]*error[0,i])
    averrorS = averrorS + (error[2,i]*error[0,i])


print ""
print "Average Picking error P (s):",averrorP/errsumP
print "Average Picking error S (s):",averrorS/errsumS

fp1.close()

#Find common events:
mergedHypo = []

#Loop over P-hypos:
for i in range(len(Phypos)):

    #Loop over S-hypos:
    for j in range(len(Shypos)):

        #Check if events match:
        if(Phypos[i].evID1 == Shypos[j].evID1):

           #Events match, add hypocenter to "mergedHypo" list and mark it as match in both hypolists:
           #Mark as match:
           Phypos[i].matchIDX = j
           Shypos[j].matchIDX = i

           #Log matching event:
           fp2.write("#Matching-ID %-20s P-cat-idx: %8d  S-cat-idx: %8d\n" % (Phypos[i].evID1,i,j))

           #Now decide on hypocenter and copy it to merged file:
           if(HypoInfo == 'P'):
              #Take hypocenter from P:
              mergedHypo.append(copy.deepcopy(Phypos[i]))
           else:
              #Take hypocenter from S:
              mergedHypo.append(copy.deepcopy(Shypos[j]))
           break

#Now check non-matching events:
   
#Add Non-Matching P-events:
#Loop over P-hypos:
print "Check non-Matching P-events (and add if ALL/ALL+ mergemode)..."
for i in range(len(Phypos)):
    #Check if event is non-matching:
    if(Phypos[i].matchIDX < 0):
       if(mergemode == 'ALL') or (mergemode == 'ALL+'):
          #Indicate that event is from P-data:
          Phypos[i].matchIDX = -1
          mergedHypo.append(copy.deepcopy(Phypos[i]))
          fp2.write("#NonMatching-P-Event-ADDED   %-20s\n" % (Phypos[i].evID1))
       else:
          fp2.write("#NonMatching-P-Event-SKIPPED %-20s\n" % (Phypos[i].evID1))

#Loop over S-hypos:
print "Check non-Matching S-events (and add if ALL/ALL+ mergemode)..."
for i in range(len(Shypos)):
    #Check if event is non-matching:
    if(Shypos[i].matchIDX < 0):
       if(mergemode == 'ALL') or (mergemode == 'ALL+'):
          #Indicate that event is from S-data:
          Shypos[i].matchIDX = -2
          mergedHypo.append(copy.deepcopy(Shypos[i]))
          fp2.write("#NonMatching-S-Event-ADDED   %-20s\n" % (Shypos[i].evID1))
       else:
          fp2.write("#NonMatching-S-Event-SKIPPED %-20s\n" % (Shypos[i].evID1))

#Now Merge:
print "--> Merge CNV-files of common events:",len(mergedHypo)

#Sort hypocenters according to OT-timestamp:
mergedHypo.sort(key=lambda x: (x.timestamp), reverse=False)

#Merge common events and write it to output
fp = open(OutputFile,'w')
for i in range(len(mergedHypo)):

    #Write hypocenter to file:   
    writeCNVHeader(fp,mergedHypo[i])

    #Merge phases and adjust traveltimes to common OT:
    tmpPhases = []

    #Add P-phases:
    if(mergedHypo[i].matchIDX >= 0) or (mergedHypo[i].matchIDX == -1):
       for j in range(len(Pphases)):
           if(mergedHypo[i].evID1 == Pphases[j].evID):
              #Update tt regarding new OT:
              Pphases[j].tt = Pphases[j].timestamp - mergedHypo[i].timestamp
              tmpPhases.append(Pphases[j])

    #Add S-phases:
    if(mergedHypo[i].matchIDX >= 0) or (mergedHypo[i].matchIDX == -2):
       for j in range(len(Sphases)):
           if(mergedHypo[i].evID1 == Sphases[j].evID):
              #Update tt regarding new OT:
              Sphases[j].tt = Sphases[j].timestamp - mergedHypo[i].timestamp
              tmpPhases.append(Sphases[j])

    #Now complement if requested (this step could be faster with indexing (adding position-range of phases to event info), but for CH-dimensions ok at the moment...):
    #Case 1): Find S-phases for a P-only event:
    if(mergemode == 'ALL+') and (mergedHypo[i].matchIDX == -1) and (len(Cphases)>0):
       for j in range(len(Cphases)):
           if(mergedHypo[i].evID1 == Cphases[j].evID) and (Cphases[j].ptype == 'S'):
              #Update tt regarding new OT:
              Cphases[j].tt = Cphases[j].timestamp - mergedHypo[i].timestamp
              fp2.write("#NonMatching-P-Event-PHASE-ADDED: %-20s %-9s %1s %8.3f\n" % (mergedHypo[i].evID1,Cphases[j].stat,Cphases[j].ptype,Cphases[j].tt))
              tmpPhases.append(Cphases[j])

    #Case 2): Find P-phases for a S-only event:
    if(mergemode == 'ALL+') and (mergedHypo[i].matchIDX == -2) and (len(Cphases)>0):
       for j in range(len(Cphases)):
           if(mergedHypo[i].evID1 == Cphases[j].evID) and (Cphases[j].ptype == 'P'):
              #Update tt regarding new OT:
              Cphases[j].tt = Cphases[j].timestamp - mergedHypo[i].timestamp
              fp2.write("#NonMatching-S-Event-PHASE-ADDED: %-20s %-9s %1s %8.3f\n" % (mergedHypo[i].evID1,Cphases[j].stat,Cphases[j].ptype,Cphases[j].tt))
              tmpPhases.append(Cphases[j])

    #Now sort according to distance and phase type - STILL MISSING - not necessary at the moment, problem is that distances could
    #be different because different hypocenters were used, sorting can be done later, not required for tomo codes anyway... 

    #Write phases to file:
    writeCNVPhases(fp,tmpPhases)
 
#    inotherfile = False
#    for j in range(len(velestlist)):
#        if(refevlist[i].evID1 == velestlist[j].evID1):
#           inotherfile = True
#           break
#    #Final check:
#    if(inotherfile == False):
#       fp.write("%10.5f %9.5f %7.3f %4.1f %-2s %04d %02d %02d %02d %02d %05.2f %-4s %6.3f %3d %5.1f %3d %1s %1s %-20s %9d | %4.0f | %4.0f | %-15s | %-25s | %s | %s\n" % (refevlist[i].lon,refevlist[i].lat,refevlist[i].dep,refevlist[i].mag,refevlist[i].mtype[0:2],refevlist[i].yy,refevlist[i].mm,refevlist[i].dd,refevlist[i].hh,refevlist[i].mi,refevlist[i].ss,refevlist[i].methodID,refevlist[i].rms,refevlist[i].gap,refevlist[i].mdist,refevlist[i].nobs,refevlist[i].etype,refevlist[i].lqual,refevlist[i].evID1,refevlist[i].evID2,refevlist[i].chx,refevlist[i].chy,refevlist[i].author,refevlist[i].region,refevlist[i].orID1,refevlist[i].evPID))        
fp.close()
fp2.close()

print ""
print "Hypocenters taken from file                         :",HypoInfo
print "File with merged P+S phases of common events        :",OutputFile
print "Number of                      common events        :",len(mergedHypo)
print "Statistics on phase data written to file            :",LogFile1
print "Log-file with event-merging information             :",LogFile2

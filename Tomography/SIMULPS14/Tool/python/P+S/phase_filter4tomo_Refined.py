#!/usr/bin/env python2

import re
import math
import time
import numpy as np
import cPickle
from obspy.core import UTCDateTime #Comented ADN
#from obspy.core.util import gps2DistAzimuth
from obspy.geodetics import gps2dist_azimuth  #Comented ADN
#import pg  #Comented ADN
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

      def __init__(self,yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,wgsum,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,stphaidx,enphaidx,matchidx,evscore,gid,useflg):
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
          self.wgsum = float(wgsum)               #SUm of weights
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
          self.catID  = int(catID)        #ID of catalog (integer number) / here: number of events in the same grid volume
          self.stphaidx = long(stphaidx)  #Position of first phase of event 
          self.enphaidx = long(enphaidx)  #Position of last  phase of event
          self.matchidx = long(matchidx)  #Position of matching event
          self.evscore = float(evscore)
          self.gid     = long(gid)        #ID of grid volume associated with this event
          self.useflg  = int(useflg)      #Use or don't use flage

      def __str__(self):
          return '%04d/%02d/%02d %02d:%02d:%06.3f %9.4f %8.4f %6.2f %4.1f %-5s %6.3f %3d %7.1f %3d %-2s %1s %3.0f %3.0f %-7s %-15s %-35s %-3s %15d %-15s %-s %-s' % (self.yy,self.mm,self.dd,self.hh,self.mi,self.ss,self.lon,self.lat,self.dep,self.mag,self.mtype,self.rms,self.gap,self.mdist,self.nobs,self.etype,self.lqual,self.chx,self.chy,self.agency,self.author,self.region,self.fstatus,self.evID2,self.methodID[0:15],self.evID1,self.orID1)

#------------------------------------------------------------------------------------------
class phase(object):
    
      def __init__(self,stat,ptype,qual,err,weig,tt,timestamp,edist,azi,baz,evID,slat,slon,stidx,oscore,phaidx):
          self.stat = str(stat)
          self.ptype = str(ptype)
          self.qual = int(qual)
          self.err  = float(err)
          self.weig = float(weig)
          self.tt   = float(tt)
          self.timestamp = float(timestamp)
          self.edist = float(edist)
          self.azi = float(azi)
          self.baz = float(baz)
          self.evID = str(evID)
          self.slat = float(slat)
          self.slon = float(slon)
          self.stidx = long(stidx)
          self.oscore = float(oscore)
          self.phaidx = long(phaidx)       #Position of corrssponding phase in phase-list -> Can be used later to adjust oscore in phase for final filtering
#------------------------------------------------------------------------------------------
class residual(object):
      def __init__(self,stat,ptype,qual,err,weig,res,dis,evID,AvResInQuad,NoResInQuad,Quad,evidx,stidx,auxf1,auxf2,auxf3,auxf4,auxf5):

          self.stat = str(stat)
          self.ptype = str(ptype)
          self.qual = int(qual)
          self.err  = float(err)
          self.weig = float(weig)
          self.res  = float(res)
          self.dis = float(dis)
          self.evID = str(evID)
          self.AvResInQuad = float(AvResInQuad)
          self.NoResInQuad = int(NoResInQuad)
          self.Quad = int(Quad)
          self.evidx = int(evidx)
          self.stidx = int(stidx)
          self.auxf1 = float(auxf1)
          self.auxf2 = float(auxf2)
          self.auxf3 = float(auxf3)          
          self.auxf4 = float(auxf4)
          self.auxf5 = float(auxf5)
#------------------------------------------------------------------------------------------
class outlier(object):
      def __init__(self,stat,statorg,pha,qua,err,picktimestamp,res,SCnob,SCval,SCasr,SCssr,SCasrQ,SCnsrQ,Quad,epi,azi,baz,slon,slat,asepi,ssepi,OTtimestamp,elon,elat,edep,emag,oscore,WFrevf,evID,phaidx):

          self.stat = str(stat)
          self.statorg = str(statorg)
          self.pha = str(pha)
          self.qua = int(qua)
          self.err = float(err)
          self.picktimestamp = float(picktimestamp)
          self.res = float(res)            #Travel-time residual
          self.SCnob = int(SCnob)          #Station-Correction: Total number of observation for station
          self.SCval = float(SCval)        #Station-Correction: Value of station-correction
          self.SCasr = float(SCasr)        #Station-Correction: Average Residual for station (all quadrants)
          self.SCssr = float(SCssr)        #Station-Correction: Standard deviation of residuals for station (all quadrants)
          self.SCasrQ = float(SCasrQ)      #Station-Correction: Average Residual for station for quadrant associated with baz of phase
          self.SCnsrQ = int(SCnsrQ)        #Station-Correction: Number of Residual for station for quadrant associated with baz of phase
          self.Quad = int(Quad)            #VELEST Quadrant (based on baz)
          self.epi = float(epi)            #Epicentral distance of phase
          self.azi = float(azi)            #Azimuth of phase
          self.baz = float(baz)            #Back-Azimuth of phase
          self.slon = float(slon)         
          self.slat = float(slat)
          self.asepi = float(asepi)        #Average epicentral distance of all phases at this station
          self.ssepi = float(ssepi)        #Standard deviation of epicentral distance of all phases at this station
          self.OTtimestamp = float(OTtimestamp)
          self.elon = float(elon)
          self.elat = float(elat)
          self.edep = float(edep)
          self.emag = float(emag)
          self.oscore = float(oscore)      #Score value to decide how likely phase is real outlier
          self.WFrevf = int(WFrevf)        #Flag indicating status of waveform review: -1: No review performed, 0: Review Performed, classiefied as no outlier, 1: Review Performed, classified as outlier, -2: Waveform not loaded/missingWFrevf
          self.evID = str(evID)
          self.phaidx = long(phaidx)       #Position of corrssponding phase in phase-list -> Can be used later to adjust oscore in phase for final filtering

#------------------------------------------------------------------------------------------
class wscheme(object):
      def __init__(self,qual,err,weig,nobsP,nobsS):
          self.qual = int(qual)
          self.err  = float(err)
          self.weig = float(weig)
          self.nobsP = int(nobsP)
          self.nobsS = int(nobsS)
#------------------------------------------------------------------------------------------
class station(object):

      def __init__(self,stat,statorg,netw,lat,lon,ele,nobsP,nobsS,nobsPw,nobsSw,SC_P_NobAll,SC_P_DelAll,SC_P_AREAll,SC_P_STDAll,SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4,SC_P_DISAve,SC_P_DISSTD,SC_P_BAZAve,SC_P_BAZSTD):
          self.stat = str(stat)                    #VELEST 4 character station name (could be alias name)
          self.statorg = str(statorg)              #Original (full) station name
          self.netw = str(netw)
          self.lat = float(lat)
          self.lon = float(lon)
          self.ele = float(ele)
          self.nobsP = int(nobsP)
          self.nobsS = int(nobsS)
          self.nobsPw = float(nobsPw)
          self.nobsSw = float(nobsSw)
          self.SC_P_NobAll = int(SC_P_NobAll)      #All P observations: Number of observations
          self.SC_P_DelAll = float(SC_P_DelAll)    #All P observations: Delay
          self.SC_P_AREAll = float(SC_P_AREAll)    #All P observations: Average (unweighted residual)
          self.SC_P_STDAll = float(SC_P_STDAll)    #All P observations: Standard deviation
          self.SC_P_DelQA1 = float(SC_P_DelQA1)    #P Quadrant 1: Delay
          self.SC_P_NobQA1 = int(SC_P_NobQA1)      #P Quadrant 1: Nobs
          self.SC_P_DelQA2 = float(SC_P_DelQA2)    #P Quadrant 2: Delay
          self.SC_P_NobQA2 = int(SC_P_NobQA2)      #P Quadrant 2: Nobs
          self.SC_P_DelQA3 = float(SC_P_DelQA3)    #P Quadrant 3: Delay
          self.SC_P_NobQA3 = int(SC_P_NobQA3)      #P Quadrant 3: Nobs
          self.SC_P_DelQA4 = float(SC_P_DelQA4)    #P Quadrant 4: Delay
          self.SC_P_NobQA4 = int(SC_P_NobQA4)      #P Quadrant 4: Nobs
          self.SC_S_NobAll = int(SC_S_NobAll)      #All S observations: Number of observations
          self.SC_S_DelAll = float(SC_S_DelAll)    #All S observations: Delay
          self.SC_S_AREAll = float(SC_S_AREAll)    #All S observations: Average (unweighted residual)
          self.SC_S_STDAll = float(SC_S_STDAll)    #All S observations: Standard deviation
          self.SC_S_DelQA1 = float(SC_S_DelQA1)    #S Quadrant 1: Delay
          self.SC_S_NobQA1 = int(SC_S_NobQA1)      #S Quadrant 1: Nobs
          self.SC_S_DelQA2 = float(SC_S_DelQA2)    #S Quadrant 2: Delay
          self.SC_S_NobQA2 = int(SC_S_NobQA2)      #S Quadrant 2: Nobs
          self.SC_S_DelQA3 = float(SC_S_DelQA3)    #S Quadrant 3: Delay
          self.SC_S_NobQA3 = int(SC_S_NobQA3)      #S Quadrant 3: Nobs
          self.SC_S_DelQA4 = float(SC_S_DelQA4)    #S Quadrant 4: Delay
          self.SC_S_NobQA4 = int(SC_S_NobQA4)      #S Quadrant 4: Nobs
          self.SC_P_DISAve = float(SC_P_DISAve)    #P distance of all observations: average
          self.SC_P_DISSTD = float(SC_P_DISSTD)    #P distance of all observations: STD
          self.SC_P_BAZAve = float(SC_P_BAZAve)    #P BAZ      of all observations: average
          self.SC_P_BAZSTD = float(SC_P_BAZSTD)    #P BAZ      of all observations: STD
#------------------------------------------------------------------------------------------
class stationDistVector(object):

      def __init__(self,stat,phas,dist):
          self.stat = str(stat)                    #VELEST 4 character station name (could be alias name)
          self.phas = str(phas)
          self.dist = float(dist)
#------------------------------------------------------------------------------------------
class gridv(object):
      def __init__(self,lomin,lomax,lamin,lamax,demin,demax,eqnto,eqncc,lqnto,lqncc,grdid):
          self.lomin = float(lomin)
          self.lomax = float(lomax)
          self.lamin = float(lamin)
          self.lamax = float(lamax)
          self.demin = float(demin)
          self.demax = float(demax)
          self.eqnto  = int(eqnto)
          self.eqncc  = int(eqncc)
          self.lqnto  = int(lqnto)
          self.lqncc  = int(lqncc)
          self.grdid  = float(grdid)
#------------------------------------------------------------------------------------------

###########################################################################################
#Subroutines/functions:

def writeSIMPhases(fp,PhasesIn):

    #Writes phase to a SIMULPS format file

    nop = len(PhasesIn)

    #SIMULPS accepts only P and S-P times, therefore we have to convert "PhasesIn" first into P and S-P "Phases"
    #Empty list for simulps format
    Phases = []

    #Sort phases corresponding distance and phase type: 
    PhasesIn.sort(key=lambda x: (x.edist, x.ptype), reverse=False)

    #Loop over all phases 
    for i in range(nop):
        if(PhasesIn[i].ptype == "P"):
           PhasesIn[i].ptype = "P-"
           Phases.append(copy.deepcopy(PhasesIn[i]))

        if(PhasesIn[i].ptype == "S"):  
           #Check if P-phase exist for same station to calculate S-P time:
           SmPtidx = -9
           for j in range(nop):
               if((PhasesIn[j].ptype == 'P') or (PhasesIn[j].ptype == 'P-')) and (PhasesIn[j].stat == PhasesIn[i].stat):
                 SmPtidx = j
                 break
           #If P-phase exists, calculate S-P time and adjust quality:
           #take the lower of the two
           if(SmPtidx >= 0):     
              PhasesIn[i].ptype = "SP"
              PhasesIn[i].tt = PhasesIn[i].tt - PhasesIn[SmPtidx].tt

              #Get the lower of the two qualities:
              if(PhasesIn[SmPtidx].qual > PhasesIn[i].qual):
                 PhasesIn[i].qual = PhasesIn[SmPtidx].qual
              
              #Final check:
              if(PhasesIn[i].tt > 0):
                 Phases.append(copy.deepcopy(PhasesIn[i]))

    #Calculate final number of obervations:
    nop = len(Phases)
    in6 = int(nop/6)
    rem = nop - (in6*6)

    #stat,ptype,qual,tt,timestamp,evID

    cnt = 0

    #Print full lines:
    #(a4,3(a1),i1,f6.2))') (pArUseC(h,1),'I',pArUseC(h,2),'-',pArUseI(h),pArUseR(h)
    for i in range(in6):
        fp.write("%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt,Phases[cnt+3].stat,Phases[cnt+3].ptype,Phases[cnt+3].qual,Phases[cnt+3].tt,Phases[cnt+4].stat,Phases[cnt+4].ptype,Phases[cnt+4].qual,Phases[cnt+4].tt,Phases[cnt+5].stat,Phases[cnt+5].ptype,Phases[cnt+5].qual,Phases[cnt+5].tt))
        cnt = cnt + 6

    #Print rest:
    if(rem == 1):
       fp.write("%-4sI%-2s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt))
       cnt = cnt + 1

    if(rem == 2):
       fp.write("%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt))
       cnt = cnt + 2

    if(rem == 3):
       fp.write("%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt))
       cnt = cnt + 3


    if(rem == 4):
       fp.write("%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt,Phases[cnt+3].stat,Phases[cnt+3].ptype,Phases[cnt+3].qual,Phases[cnt+3].tt))
       cnt = cnt + 4

    if(rem == 5):
       fp.write("%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f%-4sI%-2s%1d%6.2f\n" % (Phases[cnt+0].stat,Phases[cnt+0].ptype,Phases[cnt+0].qual,Phases[cnt+0].tt,Phases[cnt+1].stat,Phases[cnt+1].ptype,Phases[cnt+1].qual,Phases[cnt+1].tt,Phases[cnt+2].stat,Phases[cnt+2].ptype,Phases[cnt+2].qual,Phases[cnt+2].tt,Phases[cnt+3].stat,Phases[cnt+3].ptype,Phases[cnt+3].qual,Phases[cnt+3].tt,Phases[cnt+4].stat,Phases[cnt+4].ptype,Phases[cnt+4].qual,Phases[cnt+4].tt))
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
def writeSIMHeader(fp,Hypo):

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

    #Convert spherical coordinates into deg + minutes
    latd = int(lat)
    latm = (lat-latd)*60.0
    lond = int(lon)
    lonm = (lon-lond)*60.0 

    #NEEDS TO BE ADJUSTED FOR SIMULPS, STILL CNV:
    # write(12,'(3(i2.2),1x,2(i2.2),1x,f5.2,1x,i2.2,a1,f5.2,1x,i3.3,a1,f5.2,1x,f6.2,1x,f6.2)')
    #       yy,mm,dd,hh,mi,ss,lati,latc,latm,loni,lonc,lonm,rdep,mag
    fp.write("%02d%02d%02d %02d%02d %05.2f %02d%1s%05.2f %03d%1s%05.2f %6.2f %6.2f EVID: %-s\n" % (vyy,Hypo.mm,Hypo.dd,Hypo.hh,Hypo.mi,Hypo.ss,latd,latc,latm,lond,lonc,lonm,Hypo.dep,Hypo.mag,Hypo.evID1))

    return
#------------------------------------------------------------------------------------------
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
#lomin,lomax,lamin,lamax,demin,demax,ntot,ncc
def setupgrid(dx,dy,dz,gslat,gslon,gsdep,gelat,gelon,gedep):

    print ""
    print "--> Setup grid for spatial declustering"

    glist = []

    #Get mean latitude of area:
    mlat = (gslat+gelat)/2.0
    print "    Mean latitude:",mlat

    #Convert grid spacing in km to degrees:
    dxde=km2lon(dx,mlat)
    dyde=km2lat(dy)

    #Start values:
    clon=gslon
    clat=gslat
    cdep=gsdep
    gcnt=-1

    #For last depth interval use max depth:
    cdep_max = 100.0

    #get number of depth layers:
    ndepth = int((float(gedep)-float(gsdep))/float(dz))
    cntdep = 0

    #No setup grid volumes:
    #loop over depth
    while(cdep < gedep):
          
          #Count number of depth intervals:
          cntdep += 1

          #Check if final depth layer reach, if yes, adjust upper bound to max value to make
          #sure that deep seismicity is included in lowermost layer.
          if(cntdep == ndepth) and (cdep+dz < cdep_max):
             cdep_up = cdep_max
          else:
             cdep_up = cdep+dz

          print "--> Setup grid between depth",cdep,"and",cdep_up,"layer #",cntdep

          #loop over X:
          clon=gslon
          while(clon < gelon):
    
                #loop over Y:
                clat=gslat
                while(clat > gelat):

                      #Add grid to list here:
                      gcnt += 1
                      #print "Add volume",clon,clon+dxde,clat,clat+dyde,cdep,cdep+dz,0,0,0,0
                      glist.append(gridv(clon,clon+dxde,clat,clat+dyde,cdep,cdep_up,0,0,0,0,gcnt))

                      #Update Y
                      clat = clat-dyde

                #Update X:
                clon = clon+dxde

          #Update depth for next iteration:
          cdep = cdep+dz

    print "--> Number of generated  depth-layers                              :",cntdep
    print "--> Number of calculated layers (int(float(gedep)-float(gsdep)/dz)):",ndepth

    return glist
#------------------------------------------------------------------------------------------
def km2lon(km,lat):
    #Converts km in degree-longitude (at latitude lat)     
    #
    #T. Diehl, LDEO 2009/11
    #
    #Some fixed/initial values:
    pi=3.14159265358979323846

    return km/(111.1949*math.cos((pi/180)*lat))
#------------------------------------------------------------------------------------------
def km2lat(km):
    #Converts km in degree-latitude
    #
    #T. Diehl, LDEO 2009/11
    #
    #Some fixed/initial values:
    pi=3.14159265358979323846

    return km/(111.1949)
#------------------------------------------------------------------------------------------
def setupWeighting_ALPINE():
    wschemelist = []

    wschemelist.append(wscheme(0,0.050,1.0000,0,0))
    wschemelist.append(wscheme(1,0.100,0.5000,0,0))
    wschemelist.append(wscheme(2,0.200,0.2500,0,0))
    wschemelist.append(wscheme(3,0.400,0.1250,0,0))
    wschemelist.append(wscheme(4,9.990,0.0000,0,0))

    return wschemelist
#------------------------------------------------------------------------------------------
def setupWeighting_SED():
    wschemelist = []

    wschemelist.append(wscheme(0,0.025,1.0000,0,0))
    wschemelist.append(wscheme(1,0.050,0.5000,0,0))
    wschemelist.append(wscheme(2,0.100,0.2500,0,0))
    wschemelist.append(wscheme(3,0.200,0.1250,0,0))
    wschemelist.append(wscheme(4,0.400,0.0625,0,0))
    wschemelist.append(wscheme(5,9.990,0.0000,0,0))

    return wschemelist
#------------------------------------------------------------------------------------------
def getAvRes4Quad(stations,statidx,hypolist,evntidx,pha):

    avResQuad = 0.0
    noResQuad = 0
    quad = 0

    if(statidx >= 0) and (evntidx >= 0):

       #Calculate BAZ:
       distinf = gps2dist_azimuth(hypolist[evntidx].lat,hypolist[evntidx].lon,stations[statidx].lat,stations[statidx].lon)
       baz = distinf[2]

       #Now check quadrant:
       # VELEST Definition:
       # Residuals of the stations according to the azimuth
       # (in right-handed coordinate system):
       # (RES  = total average residual at station)
       # (RES1 = average residual of rays from 1st quadrant)
       # (1st quadrant, 0 - 90 deg --> N -> W)
       # (RES2 = average residual of rays from 2nd quadrant)
       # (2nd quadrant, 90 - 180 deg --> W -> S)
       # (RES3 = average residual of rays from 3rd quadrant)
       # (3rd quadrant, 180 - 270 deg --> S -> E)
       # (RES4 = average residual of rays from 4th quadrant)
       # (4th quadrant, 270 - 360 deg --> E -> N)

       #N to E:
       if(baz >= 0.0) and (baz < 90.0):
          quad = 4
          if(pha == 'P') or (pha == 'p'):
             avResQuad = stations[statidx].SC_P_DelQA4
             noResQuad = stations[statidx].SC_P_NobQA4 
          if(pha == 'S') or (pha == 's'):
             avResQuad = stations[statidx].SC_S_DelQA4
             noResQuad = stations[statidx].SC_S_NobQA4

       #E to S:
       if(baz >= 90.0) and (baz < 180.0):
          quad = 3
          if(pha == 'P') or (pha == 'p'):
             avResQuad = stations[statidx].SC_P_DelQA3
             noResQuad = stations[statidx].SC_P_NobQA3
          if(pha == 'S') or (pha == 's'):
             avResQuad = stations[statidx].SC_S_DelQA3
             noResQuad = stations[statidx].SC_S_NobQA3

       #S to W:
       if(baz >= 180.0) and (baz < 270.0):
          quad = 2
          if(pha == 'P') or (pha == 'p'):
             avResQuad = stations[statidx].SC_P_DelQA2
             noResQuad = stations[statidx].SC_P_NobQA2
          if(pha == 'S') or (pha == 's'):
             avResQuad = stations[statidx].SC_S_DelQA2
             noResQuad = stations[statidx].SC_S_NobQA2

       #W to N:
       if(baz >= 270.0) and (baz <= 360.0):
          quad = 1
          if(pha == 'P') or (pha == 'p'):
             avResQuad = stations[statidx].SC_P_DelQA1
             noResQuad = stations[statidx].SC_P_NobQA1
          if(pha == 'S') or (pha == 's'):
             avResQuad = stations[statidx].SC_S_DelQA1
             noResQuad = stations[statidx].SC_S_NobQA1

    return avResQuad,noResQuad,quad
#------------------------------------------------------------------------------------------
def getID(evID,evIDlist):
    idx = -1
    for n in range(len(evIDlist)):
        if(evIDlist[n] == evID):
           idx = n
           break
    return idx
#------------------------------------------------------------------------------------------
def getevntidx(hypolist,evID2find):

    idx = -1

    for n in range(len(hypolist)):
        if(hypolist[n].evID1 == evID2find):
           idx = n
           break

    return idx
#------------------------------------------------------------------------------------------
def getstatidx(stations,stat2find):

    idx = -1

    for n in range(len(stations)):
        if(stations[n].stat == stat2find):
           idx = n
           break

    return idx
#------------------------------------------------------------------------------------------
def getWeightFromError(err,weights):
    weight = 0.0  #By Default: lowest weight
    for i in range(len(weights)):
        if(err <= weights[i].err):
           weight = weights[i].weig
           break
            
    return weight
#------------------------------------------------------------------------------------------
def getQualityFromError(err,weights):
    qual = len(weights) - 1  #By Default: lowest quality
    for i in range(len(weights)):
        if(err <= weights[i].err):
           qual = weights[i].qual 
           break

    return qual
#------------------------------------------------------------------------------------------
def getGridID4Event(lat,lon,dep,grid):
    gridID = -9

    for k in range(len(grid)):
        if(lon >= grid[k].lomin) and (lon < grid[k].lomax) and (lat > grid[k].lamin) and (lat <= grid[k].lamax) and (dep > grid[k].demin) and (dep <= grid[k].demax):

           #Assign grid ID to event:
           gridID = grid[k].grdid

           #Update event counter for grid point:
           grid[k].eqnto += 1
  
           #finished:
           break

    return grid,gridID
#------------------------------------------------------------------------------------------
def getEventScore(wgap,wmdis,wwgsum,wnobs,cgap,cmdistz,cwgsum,cnobs):

    evscore = 0.0

    gap_sc = 1.0/(1.0+math.exp((wgap-cgap)*0.025))
    mdi_sc = 1.0/(1.0+math.exp((wmdis-(cmdistz*1.5))*0.20))
    wsu_sc = 1.0/(1.0+math.exp((wwgsum-cwgsum)*-0.25))
    nob_sc = 1.0/(1.0+math.exp((wnobs-cnobs)*-0.25))

    evscore = gap_sc + mdi_sc + wsu_sc + nob_sc

    return evscore
#------------------------------------------------------------------------------------------
def getDistAzi4phase(phase,elat,elon,stations,maxepidist,maxepidistS,DistVect):

    #default
    edist = -9.0
    azi = -9.0

    #Get position of station in station list:
    statidx = getstatidx(stations,phase.stat)

    if(statidx >= 0):
       distinf = gps2dist_azimuth(elat,elon,stations[statidx].lat,stations[statidx].lon)
       phase.edist = distinf[0]/1000.0
       phase.azi = distinf[1]
       phase.baz = distinf[2]
       phase.slat = stations[statidx].lat
       phase.slon = stations[statidx].lon
       phase.stidx = statidx

       #Calculate statistics for station corrections:
       stations[statidx].SC_P_DISAve = stations[statidx].SC_P_DISAve + distinf[0]/1000.0
       stations[statidx].SC_P_BAZAve = stations[statidx].SC_P_BAZAve + distinf[2]

       DistVect.append(stationDistVector(phase.stat,phase.ptype,distinf[0]/1000.0))

       #Check if distance is within allowed range and used for further processing:
       if(distinf[0]/1000.0 > maxepidist):
          #Disable it for further use:
          phase.qual = 9
          phase.weig = 0.0
          phase.oscore = 1.0
       #Check if distance is within allowed range and used for further processing for S-phases:
       if(distinf[0]/1000.0 > maxepidistS) and (phase.ptype.startswith("S")) and (phase.weig > 0.0):
          #print "--> Disable S-phase >",maxepidistS,distinf[0]/1000.0,phase.ptype
          #Disable it for further use:
          phase.qual = 9
          phase.weig = 0.0
          phase.oscore = 1.0
    else:
       phase.stidx = statidx   #Station not found (-1)
       #Disable it for further use:
       phase.qual = 9
       phase.weig = 0.0
       phase.oscore = 1.0

    return phase,DistVect
#------------------------------------------------------------------------------------------
def updateEventPhaseData(Phases,commonevID,weights):

    #Simplified version of getDistAziGap (assuming dist,epi,statinfo already linked to phse in previous step (getDistAzi4phase)

    azivec = []
    disvec = []
    weightsum = 0.0
    newObs = 0
    newpgap = 359.9
    mdist = 999.9
    phasesout = []

    for i in range(len(Phases)):

        #Adjust weight according to provided weighting scheme:
        Phases[i].weig = getWeightFromError(Phases[i].err,weights)
        Phases[i].qual = getQualityFromError(Phases[i].err,weights)

        #Use one common event ID (the one used in the hypocenter catalog) for
        #all the phases associated with this event:
        Phases[i].evID = commonevID

        #Update sum of weights: (now only used for known stations)
        #weightsum = weightsum + Phases[i].weig

        #Check if station info linked to phase, if not ignore it
        if(Phases[i].stidx >= 0):

           #Update sum of weights:
           weightsum = weightsum + Phases[i].weig

           #Use phase for gap only if weight is > 0
           if(Phases[i].weig > 0.0001):
              azivec.append(Phases[i].azi)
              disvec.append(Phases[i].edist)
              phasesout.append(Phases[i])
        else:
           print "--> WARNING: Station not included in station list:",Phases[i].stat,"--> not considered for gap calculation"

    #Calculate GAP:
    if(len(azivec) > 1):
       newpgap,newsgap = getGap(azivec)

    #Get minimum distance
    if(len(disvec) > 0):
       mdist = getMinimumDistance(disvec)

    #Number of observation with weight > 0
    if(len(disvec) > 0):
       newObs = len(disvec)

    return newpgap,mdist,weightsum,newObs,phasesout

#------------------------------------------------------------------------------------------
def getDistAziGap(Phases,commonevID,elat,elon,stations,weights):

    #Was split into two routines: getDistAzi4phase + updateEventPhaseData

    azivec = []
    disvec = []
    weightsum = 0.0
    newObs = 0
    newpgap = 359.9
    mdist = 999.9
    phasesout = []

    for i in range(len(Phases)):

        #Get position of station in station list:
        statidx = getstatidx(stations,Phases[i].stat)

        #print ""
        #print Phases[i].weig,Phases[i].qual,Phases[i].err

        #Adjust weight according to provided weighting scheme:
        Phases[i].weig = getWeightFromError(Phases[i].err,weights)
        Phases[i].qual = getQualityFromError(Phases[i].err,weights)

        #Use one common event ID (the one used in the hypocenter catalog) for
        #all the phases associated with this event:
        Phases[i].evID = commonevID

        #print Phases[i].weig,Phases[i].qual,Phases[i].err

        #Update sum of weights:
        weightsum = weightsum + Phases[i].weig

        if(statidx >= 0):
           distinf = gps2dist_azimuth(elat,elon,stations[statidx].lat,stations[statidx].lon)
           Phases[i].edist = distinf[0]/1000.0
           Phases[i].azi = distinf[1]
           Phases[i].slat = stations[statidx].lat
           Phases[i].slon = stations[statidx].lon
           Phases[i].stidx = statidx
           #Use phase only if weight is > 0
           if(Phases[i].weig > 0.0001):
              azivec.append(distinf[1])
              disvec.append(distinf[0]/1000.0)
              phasesout.append(Phases[i])
        else:
           print "--> WARNING: Station not included in station list:",Phases[i].stat,"--> not considered for gap calculation"
           Phases[i].edist = -9.0
           Phases[i].azi = -9.0

    #Calculate GAP:
    if(len(azivec) > 1):
       newpgap,newsgap = getGap(azivec)

    #Get minimum distance
    if(len(disvec) > 0):
       mdist = getMinimumDistance(disvec)

    #Number of observation with weight > 0
    if(len(disvec) > 0):
       newObs = len(disvec)
        
    return newpgap,mdist,weightsum,newObs,phasesout
#------------------------------------------------------------------------------------------
def getGap(azivec):

    gap_primary_max = -1.0
    gap_secondary_max = -1.0

    #Sort azimuth vector:
    azivec.sort()

    #Number of observations:
    naz = len(azivec)

    if(naz > 0):
       az_last2 = azivec[naz-2] - 360.0
       az_last = azivec[naz-1] - 360.0

       for i in range(naz):
           az = azivec[i]

           gap_primary = az - az_last
           if(gap_primary > gap_primary_max):
              gap_primary_max = gap_primary

           gap_secondary = az - az_last2
           if(gap_secondary > gap_secondary_max):
              gap_secondary_max = gap_secondary

           az_last2 = az_last
           az_last = az

    else:
       gap_primary_max = 359.9
       gap_secondary_max = 359.9

    return gap_primary_max,gap_secondary_max
#------------------------------------------------------------------------------------------
def getMinimumDistance(disvec):
    #print len(disvec)
    disvec.sort()
    return disvec[0]
#------------------------------------------------------------------------------------------
def loadOutliers(ifile,oscoreMAX,weights):

    outliers = []
    reslist = []

    #Open
    fp = open(ifile,"r")

    for line in fp:

        #Test if line is empty:
        if(len(line.split())<1):
           print "Line not considered:", line

        else:
           t = string.split(line,None,-1)

           #Check if line is commented:
           if(t[0].startswith("#")):
              print "Line not considered:", line
              continue

           #Append outlier
           #stat,statorg,pha,qua,err,picktimestamp,res,SCnob,SCval,SCasr,SCssr,SCasrQ,SCnsrQ,Quad,epi,azi,baz,slon,slat,asepi,ssepi,OTtimestamp,elon,elat,edep,emag,oscore,WFrevf,evID,phaidx
           #WFrevf == 1   ->   reviewed outlier
           #WFrevf == -1  -> Unreviewed outlier (if oscore > oscoreMAX)
           if(int(t[31]) > 0) or ((int(t[31]) < 0) and (float(t[30]) > oscoreMAX)):
              #Store outlier
              outliers.append(outlier(t[1],t[29],t[2],t[3],t[4],UTCDateTime(t[5]).timestamp,t[6],t[8],t[9],t[10],t[11],t[12],t[13],t[14],t[16],t[17],t[18],t[19],t[20],t[21],t[22],UTCDateTime(t[24]).timestamp,t[25],t[26],t[27],t[28],t[30],t[31],t[32],-1))
              #Store residual
              #stat,ptype,qual,err,weig,res,dis,evID,AvResInQuad,NoResInQuad,Quad
              reslist.append(residual(t[1],t[2],t[3],t[4],getWeightFromError(float(t[4]),weights),t[6],t[16],t[32],t[12],t[13],t[14],-1,-1,-9.0,-9.0,-9.0,-9.0,-9.0))

    fp.close()

    return outliers,reslist
#------------------------------------------------------------------------------------------
def loadRES4check(ResidualFile,checkquakes,wscheme,stations,hypolist):

    #Load all residual for events listed in checkquakes - list

    reslist = []

    if(os.path.isfile(ResidualFile)):

       #Open
       fp = open(ResidualFile,"r")

       for line in fp:

          if(line[0:1]==" "):
             print "    Ignore line:",line
             continue

          #Split string:
          s = string.split(line)

          #Now check if eventID is in list of events to check, if not check next line:
          if(getID(s[len(s)-1],checkquakes)<0):
             continue 

          if(s[1]=="p") or (s[1]=="P"):
             pha="P"
 
          if(s[1]=="s") or (s[1]=="S"):
             pha="S"

          #Station name:
          sta = str(s[0])

          #Get pick uncertainty:
          errph = wscheme[int(s[2])].err
          wgtph = wscheme[int(s[2])].weig

          #Get position of station in station list:
          statidx = getstatidx(stations,sta)

          #Get information on event:
          evntidx = getevntidx(hypolist,s[len(s)-1])

          #Calculate BAZ to get correct VELEST quadrant and Average Residual for that quadrant:
          if(statidx >= 0) and (evntidx >= 0):
             #Calculate BAZ to get correct quadrant:
             AvResInQuad,NoResInQuad,Quad = getAvRes4Quad(stations,statidx,hypolist,evntidx,pha)
          else:
             #Define default values:
             AvResInQuad = 0.0
             NoResInQuad = 0
             Quad = 0

          #Append residual no matter what:
          reslist.append(residual(s[0],pha,s[2],errph,wgtph,s[3],s[5],s[len(s)-1],AvResInQuad,NoResInQuad,Quad,-1,-1,-9.0,-9.0,-9.0,-9.0,-9.0)) 

       fp.close()
    else:
       print "--> Residual File not existing:",ResidualFile,'--> Skip residual filter'
       return reslist

    return reslist
#------------------------------------------------------------------------------------------
def loadRES(ResidualFile,co_maxRES,wscheme,stations,hypolist):

    #Different ways to calculate the residual thresholds dynamically:
    #1: residual: residual; restr = abs(errph)*abs(co_maxRES)
    #2: residual: |AverageStationResidual - IndividualResidual| >= picking error
    dynamicthrsmeth = 2

    reslist = []
    reslist_All = []

    if(os.path.isfile(ResidualFile)):

       #Open
       fp = open(ResidualFile,"r")

       for line in fp:
       
          if(line[0:1]==" "):
             print "    Ignore line:",line
             continue
          
          #Split string:
          s = string.split(line)

          if(s[1]=="p") or (s[1]=="P"):
             pha="P"
 
          if(s[1]=="s") or (s[1]=="S"):
             pha="S"

          #Station name:
          sta = str(s[0])

          #Get pick uncertainty:
          errph = wscheme[int(s[2])].err
          wgtph = wscheme[int(s[2])].weig

          #Get position of station in station list:
          statidx = getstatidx(stations,sta)

          #Get information on event:
          evntidx = getevntidx(hypolist,s[len(s)-1])

          #Calculate BAZ to get correct VELEST quadrant and Average Residual for that quadrant:
          if(statidx >= 0) and (evntidx >= 0):
             #Calculate BAZ to get correct quadrant:
             AvResInQuad,NoResInQuad,Quad = getAvRes4Quad(stations,statidx,hypolist,evntidx,pha)
          else:
             #Define default values:
             AvResInQuad = 0.0
             NoResInQuad = 0
             Quad = 0

          # Save basic information to complete list of residuals
          # stat,ptype,qual,err,weig,res,dis,evID,AvResInQuad,NoResInQuad,Quad
          reslist_All.append(residual(s[0],pha,s[2],errph,wgtph,s[3],s[5],s[len(s)-1],AvResInQuad,NoResInQuad,Quad,-1,-1,-9.0,-9.0,-9.0,-9.0,-9.0))
     
          #Get residual threshold (qual,err,weig,nobsP,nobsS):
          #AverageStationResidual - IndividualResidual >= IndividualObservationError
          if(co_maxRES < 0):
             #Calculate threshold dynamically:
             if(dynamicthrsmeth == 1): 
                resva = abs(float(s[3]))
                restr = abs(errph)*abs(co_maxRES)
             if(dynamicthrsmeth == 2):
                if(statidx >= 0) and (evntidx >= 0):
                   #SC_P_AREAll,SC_P_STDAll
                   if((pha=='P') or (pha=='p')):

                      #Check if Average residual AND STD is available (requires >= 2 observations)
                      if(stations[statidx].SC_P_NobAll >= 2):

                         #Use AverageStationResidual - IndividualResidual as residual - Edi's original proposal
                         #resva = abs(stations[statidx].SC_P_AREAll - float(s[3]))
                         #restr = abs(errph)

                         ########################
                         #Check if pick+uncertainty overlapps with predicted+uncertainty:
                         #resva = abs(float(s[3])) - stations[statidx].SC_P_AREAll - stations[statidx].SC_P_STDAll - abs(errph)
                         if(float(s[3]) < 0.0):
                            resva = abs(float(s[3])) + AvResInQuad - stations[statidx].SC_P_STDAll - abs(errph)
                            #resva = abs(float(s[3])) + AvResInQuad - 1.0*abs(errph)
                         else:
                            resva = float(s[3]) - AvResInQuad - stations[statidx].SC_P_STDAll - abs(errph)
                            #resva = float(s[3]) - AvResInQuad - 1.0*abs(errph)
                         restr = 0.0

                      #Check if Average residual is available (requires == 1 observations)
                      if(stations[statidx].SC_P_NobAll == 1):

                         #Use AverageStationResidual - IndividualResidual as residual - Edi's original proposal
                         #resva = abs(stations[statidx].SC_P_AREAll - float(s[3]))
                         #restr = abs(errph)

                         ########################
                         #Check if pick+uncertainty overlapps with predicted+uncertainty:
                         #If STD is not available use only average residual:
                         #resva = abs(float(s[3])) - stations[statidx].SC_P_AREAll - abs(errph)
                         if(float(s[3]) < 0.0):
                            resva = abs(float(s[3])) + AvResInQuad - 1.0*abs(errph)
                         else:
                            resva = float(s[3]) - AvResInQuad - 1.0*abs(errph)
                         restr = 0.0     

                      #No observation available:                      
                      if(stations[statidx].SC_P_NobAll < 1):
                         #If no observation:
                         AvResInQuad = 0.0
                         NoResInQuad = 0
                         Quad = 0
                         restr = 0.0
                         if(float(s[3]) < 0.0):
                            resva = abs(float(s[3])) - 1.0*abs(errph)
                         else:
                            resva = float(s[3]) - 1.0*abs(errph)
 
                   #if((pha=='S') or (pha=='s')) and (stations[statidx].SC_S_NobAll >= 1):
                   #   #Use AverageStationResidual - IndividualResidual as residual
                   #   resva = abs(stations[statidx].SC_S_AREAll - float(s[3]))
                   #else:
                   #   #If average residual not available use simply residual as residual-value:
                   #   resva = abs(float(s[3]))

                   #Debug:
                   #print "Residual          :",resva,sta,pha
                   #print "Threshold-Residual:",restr

                #If station is missing use residual and picking error to be conservative:
                else:            
                   AvResInQuad = 0.0
                   NoResInQuad = 0
                   Quad = 0    
                   resva = abs(float(s[3])) - 1.0*abs(errph)
                   restr = 0.0

          #Use static cut-off of residuals:    
          else:
             restr = co_maxRES
             resva = abs(float(s[3]))

          #          ROM_ p  1   0.182    6.50   35.22 960116 0700 20.80 19960116_0700
          #Reminder: stat,ptype,qual,err,weig,res,dis,evID
          if(resva >= restr):
             #print "--> Potential outlier-residual",s[0],pha,s[2],s[3],errph,wgtph,restr,s[len(s)-1]
             reslist.append(residual(s[0],pha,s[2],errph,wgtph,s[3],s[5],s[len(s)-1],AvResInQuad,NoResInQuad,Quad,-1,-1,-9.0,-9.0,-9.0,-9.0,-9.0)) 
        
       fp.close()
    else:
       print "--> Residual File not existing:",ResidualFile,'--> Skip residual filter'
       return reslist

    return reslist,reslist_All
#------------------------------------------------------------------------------------------
def check4residual(evID,stat,ptype,qual,reslist):

    resflg = 0
    resval = -99.0
    asrQ = 0.0
    nsrQ = 0
    Quad = 0
    residx = -1

    for i in range(len(reslist)):

        if(reslist[i].evID == evID) and (reslist[i].stat == stat) and (reslist[i].ptype == ptype) and (reslist[i].qual == qual):
           resflg = 1
           resval = reslist[i].res
           asrQ = reslist[i].AvResInQuad
           nsrQ = reslist[i].NoResInQuad
           Quad = reslist[i].Quad
           #print "--> Remove phase:",evID,stat,ptype,qual,"Residual:",reslist[i].res
           residx = i
           break

    return resflg,resval,residx,asrQ,nsrQ,Quad
#------------------------------------------------------------------------------------------
def loadstatcor(stations,ifile):

    #Check if file exists, otherwise skip...
    if(os.path.isfile(ifile)):

       #Open
       fp = open(ifile,"r")

       #Loop over lines:
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

           #Station correction P:
           delP = float(line[35:40])
           nobP = int(line[54:59])

           #Station correction S:
           delS = float(line[42:47])
           nobS = int(line[60:65])

           #Find station in station file:
           idx = getstatidx(stations,stat)

           #stat,netw,lat,lon,ele,nobsP,nobsS,nobsPw,nobsSw,SC_P_NobAll,SC_P_DelAll,SC_P_STDAll,SC_P_DelQA1,SC_P_DelQA2,SC_P_DelQA3,SC_P_DelQA4,SC_S_NobAll,SC_S_DelAll,SC_S_STDAll,SC_S_DelQA1,SC_S_DelQA2,SC_S_DelQA3,SC_S_DelQA4
           #Check is station exists in previous file:
           if(idx >= 0):
              stations[idx].SC_P_NobAll = nobP
              stations[idx].SC_P_DelAll = delP

              stations[idx].SC_S_NobAll = nobS
              stations[idx].SC_S_DelAll = delS

           else:
              print "WARNING: Station missing in original station file - add it from VELEST output:",stat

              #Some statistics related to station corretions (SC):
              SC_P_NobAll = 0
              SC_P_DelAll = 0.0
              SC_P_AREAll = 0.0
              SC_P_STDAll = 0.0
              SC_P_DelQA1 = 0.0
              SC_P_NobQA1 = 0
              SC_P_DelQA2 = 0.0
              SC_P_NobQA2 = 0
              SC_P_DelQA3 = 0.0
              SC_P_NobQA3 = 0
              SC_P_DelQA4 = 0.0
              SC_P_NobQA4 = 0
              SC_S_NobAll = 0
              SC_S_DelAll = 0.0
              SC_S_AREAll = 0.0
              SC_S_STDAll = 0.0
              SC_S_DelQA1 = 0.0
              SC_S_NobQA1 = 0
              SC_S_DelQA2 = 0.0
              SC_S_NobQA2 = 0
              SC_S_DelQA3 = 0.0
              SC_S_NobQA3 = 0
              SC_S_DelQA4 = 0.0
              SC_S_NobQA4 = 0
              SC_P_DISAve = 0.0
              SC_P_DISSTD = 0.0
              SC_P_BAZAve = 0.0
              SC_P_BAZSTD = 0.0

              #SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4
              statlist.append(station(stat,'XXXXX',netw,lat,lon,elev,0,0,0.0,0.0,SC_P_NobAll,SC_P_DelAll,SC_P_AREAll,SC_P_STDAll,SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4,SC_P_DISAve,SC_P_DISSTD,SC_P_BAZAve,SC_P_BAZSTD))

       #Close file
       fp.close()

    return stations
#------------------------------------------------------------------------------------------    
def loadstatcorstati(stations,ifile):

    #Check if file exists, otherwise skip...
    if(os.path.isfile(ifile)):

       #Open
       fp = open(ifile,"r")

       #Loop over lines:
       for line in fp:
          if(line[0:16]==" sta phase  nobs"):
              print "    Ignore line:",line
              continue
          if(len(line)!=61):
              print "    Ignore line:",line,len(line)
              continue

          #sta phase nobs avres  avwres    std    wsum    delay
          s = string.split(line)

          if(len(s) != 8):
             print "    Ignore line:",line
             continue  

          sta = string.strip(s[0])   
          phase = string.strip(s[1])
          nobs = int(string.strip(s[2]))
          avres = float(string.strip(s[3]))
          avwres = float(string.strip(s[4]))
          std = float(string.strip(s[5]))
          wsum = float(string.strip(s[6]))
          delay = float(string.strip(s[7]))

          #Find station in station file:
          idx = getstatidx(stations,sta)

          #Check is station exists in previous file:
          if(idx >= 0):
             if(phase=='P'):
                stations[idx].SC_P_STDAll = std
                stations[idx].SC_P_AREAll = avres
             if(phase=='S'):
                stations[idx].SC_S_STDAll = std
                stations[idx].SC_S_AREAll = avres
          else:
             print "WARNING: Station missing in original station file - ignore STD in VELEST output for:",sta

       #Close file
       fp.close()

    return stations
#------------------------------------------------------------------------------------------
def loadstatcorazi(stations,ifile):

    #Check if file exists, otherwise skip...
    if(os.path.isfile(ifile)):

       #Open
       fp = open(ifile,"r")

       #Loop over lines:
       for line in fp:

          if(line[0:18]==" Stn#  Stn     RES"):
              print "    Ignore line:",line
              continue
 
          if(len(line)!=78):
              if(len(line)!=16):
                 print "    Ignore line:",line,len(line)
              continue

          # Stn#  Stn     RES          RES1         RES2         RES3         RES4
          s = string.split(line)

          sta = string.strip(s[1])
          aral = float(line[13:19])
          noal = int(line[20:24])
          arQ1 = float(line[26:32])
          noQ1 = int(line[33:37])
          arQ2 = float(line[39:45])
          noQ2 = int(line[46:50])
          arQ3 = float(line[52:58])
          noQ3 = int(line[59:63])
          arQ4 = float(line[65:71])
          noQ4 = int(line[72:76])

          #Find station in station file:
          idx = getstatidx(stations,sta)

          #Check is station exists in previous file:
          if(idx >= 0):

             #Not sure how to distiguish S from P info, at the moment assume all is P...
             stations[idx].SC_P_DelQA1 = arQ1
             stations[idx].SC_P_NobQA1 = noQ1
             stations[idx].SC_P_DelQA2 = arQ2
             stations[idx].SC_P_NobQA2 = noQ2
             stations[idx].SC_P_DelQA3 = arQ3
             stations[idx].SC_P_NobQA3 = noQ3
             stations[idx].SC_P_DelQA4 = arQ4
             stations[idx].SC_P_NobQA4 = noQ4

             #stat,netw,lat,lon,ele,nobsP,nobsS,nobsPw,nobsSw,SC_P_NobAll,SC_P_DelAll,SC_P_AREAll,SC_P_STDAll,SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4
          
          else:
             print "WARNING: Station missing in original station file - ignore AZI-RES in VELEST output for:",sta

       #Close file
       fp.close()

    return stations
#------------------------------------------------------------------------------------------
def check_ST_avresperquadrant(stations,arpQthr,phaseTyp,fplog):

    #Minimum number of observations per azimuth bin:
    #nobsminperQ = 2
    nobsminperQ = 1

    #header for file:
    fplog.write("#Stat       lon       lat    ele        Phase     NobsTot AvReTot  STDTot |  AvReQ1 NobQ1 |   AvReQ2 NobQ2 |   AvReQ3 NobQ3 |   AvReQ4 NobQ4 |           AbsMaxDiff | Statcorr\n")

    stations2check = []

    for i in range(len(stations)):

        avreslist = []

        if(phaseTyp == 'P') and (stations[i].SC_P_NobAll > 1):

           if(stations[i].SC_P_NobQA1 >= nobsminperQ):
              avreslist.append(stations[i].SC_P_DelQA1)

           if(stations[i].SC_P_NobQA2 >= nobsminperQ):
              avreslist.append(stations[i].SC_P_DelQA2)

           if(stations[i].SC_P_NobQA3 >= nobsminperQ):
              avreslist.append(stations[i].SC_P_DelQA3)

           if(stations[i].SC_P_NobQA4 >= nobsminperQ):
              avreslist.append(stations[i].SC_P_DelQA4)

        if(phaseTyp == 'S') and (stations[i].SC_S_NobAll > 1):

           if(stations[i].SC_S_NobQA1 >= nobsminperQ):
              avreslist.append(stations[i].SC_S_DelQA1)

           if(stations[i].SC_S_NobQA2 >= nobsminperQ):
              avreslist.append(stations[i].SC_S_DelQA2)

           if(stations[i].SC_S_NobQA3 >= nobsminperQ):
              avreslist.append(stations[i].SC_S_DelQA3)

           if(stations[i].SC_S_NobQA4 >= nobsminperQ):
              avreslist.append(stations[i].SC_S_DelQA4)   

        #Now calculate differences if enough observations:
        if(len(avreslist) > 1):
           #Sort residuals
           avreslist.sort()

           #Calculate maximum difference (smallest and largest value)
           if(abs(avreslist[0]-avreslist[len(avreslist)-1]) > arpQthr):
              stations2check.append(stations[i])              
              print "WARNING: Abs. difference in average station residual between azi quadrants",abs(avreslist[0]-avreslist[len(avreslist)-1])," > ",arpQthr," for station ",stations[i].stat
              fplog.write("%-6s %9.4f %8.4f %6.3f Phase: %s NObsAll: %6d %7.3f %7.3f | %7.3f %5d |  %7.3f %5d |  %7.3f %5d |  %7.3f %5d | AbsMaxDiff: %8.3f | %8.3f\n" % (stations[i].stat,stations[i].lon,stations[i].lat,stations[i].ele,phaseTyp,stations[i].SC_P_NobAll,stations[i].SC_P_AREAll,stations[i].SC_P_STDAll,stations[i].SC_P_DelQA1,stations[i].SC_P_NobQA1,stations[i].SC_P_DelQA2,stations[i].SC_P_NobQA2,stations[i].SC_P_DelQA3,stations[i].SC_P_NobQA3,stations[i].SC_P_DelQA4,stations[i].SC_P_NobQA4,abs(avreslist[0]-avreslist[len(avreslist)-1]),stations[i].SC_P_DelAll))

    return stations2check
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

        #Some statistics related to station corretions (SC):
        SC_P_NobAll = 0
        SC_P_DelAll = 0.0
        SC_P_AREAll = 0.0
        SC_P_STDAll = 0.0
        SC_P_DelQA1 = 0.0
        SC_P_NobQA1 = 0
        SC_P_DelQA2 = 0.0
        SC_P_NobQA2 = 0
        SC_P_DelQA3 = 0.0
        SC_P_NobQA3 = 0
        SC_P_DelQA4 = 0.0
        SC_P_NobQA4 = 0
        SC_S_NobAll = 0
        SC_S_DelAll = 0.0
        SC_S_AREAll = 0.0
        SC_S_STDAll = 0.0
        SC_S_DelQA1 = 0.0
        SC_S_NobQA1 = 0
        SC_S_DelQA2 = 0.0
        SC_S_NobQA2 = 0
        SC_S_DelQA3 = 0.0
        SC_S_NobQA3 = 0
        SC_S_DelQA4 = 0.0
        SC_S_NobQA4 = 0
        SC_P_DISAve = 0.0
        SC_P_DISSTD = 0.0
        SC_P_BAZAve = 0.0
        SC_P_BAZSTD = 0.0

        #SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4
        statlist.append(station(stat,'XXXXX',netw,lat,lon,elev,0,0,0.0,0.0,SC_P_NobAll,SC_P_DelAll,SC_P_AREAll,SC_P_STDAll,SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4,SC_P_DISAve,SC_P_DISSTD,SC_P_BAZAve,SC_P_BAZSTD))

    fp.close()

    return statlist
#------------------------------------------------------------------------------------------
def loadAliasStations(ifile):

    statlist = []

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
              elev   = float(line[33:39])
              netw   = ""
              loc    = ""

              #Some statistics related to station corretions (SC):
              SC_P_NobAll = 0
              SC_P_DelAll = 0.0
              SC_P_AREAll = 0.0
              SC_P_STDAll = 0.0
              SC_P_DelQA1 = 0.0
              SC_P_NobQA1 = 0
              SC_P_DelQA2 = 0.0
              SC_P_NobQA2 = 0
              SC_P_DelQA3 = 0.0
              SC_P_NobQA3 = 0
              SC_P_DelQA4 = 0.0
              SC_P_NobQA4 = 0
              SC_S_NobAll = 0
              SC_S_DelAll = 0.0
              SC_S_AREAll = 0.0
              SC_S_STDAll = 0.0
              SC_S_DelQA1 = 0.0
              SC_S_NobQA1 = 0
              SC_S_DelQA2 = 0.0
              SC_S_NobQA2 = 0
              SC_S_DelQA3 = 0.0
              SC_S_NobQA3 = 0
              SC_S_DelQA4 = 0.0
              SC_S_NobQA4 = 0
              SC_P_DISAve = 0.0
              SC_P_DISSTD = 0.0
              SC_P_BAZAve = 0.0
              SC_P_BAZSTD = 0.0

              #SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4
              statlist.append(station(staali,staorg,netw,lat,lon,elev,0,0,0.0,0.0,SC_P_NobAll,SC_P_DelAll,SC_P_AREAll,SC_P_STDAll,SC_P_DelQA1,SC_P_NobQA1,SC_P_DelQA2,SC_P_NobQA2,SC_P_DelQA3,SC_P_NobQA3,SC_P_DelQA4,SC_P_NobQA4,SC_S_NobAll,SC_S_DelAll,SC_S_AREAll,SC_S_STDAll,SC_S_DelQA1,SC_S_NobQA1,SC_S_DelQA2,SC_S_NobQA2,SC_S_DelQA3,SC_S_NobQA3,SC_S_DelQA4,SC_S_NobQA4,SC_P_DISAve,SC_P_DISSTD,SC_P_BAZAve,SC_P_BAZSTD))
    
    #Close
    fp.close()

    return statlist
#------------------------------------------------------------------------------------------
def loadCNVhypos(ifile,catID,latmin,latmax,lonmin,lonmax):

    evlist = []

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
           wsu=-9.0
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

           #Check for formating problems in absolute time:
           corrsec = 0.0

           #Check second:
           if(ss >= 60.0):
              ss = ss - 60.0
              corrsec = corrsec + 60.0
              print '    Second >= 60.0 --> substract 60 s from ss and add 60 s          later for event:',ID1,corrsec
              #ss = ss - 0.00001

           #Check minute:
           if(mi >= 60.0):
              mi = mi - 60.0
              corrsec = corrsec + (60.0*60.0)
              print '    Minute >= 60.0 --> substract 60 m from mi and add 60 m (3600 s) later for event:',ID1,corrsec

           #print "%04d-%02d-%02d %02d:%02d:%09.6f %s" % (yy,mm,dd,hh,mi,ss,ID1)
           wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.999995 -> 60.0 Fix: increase precision

           #Now add possible corrections due to formating errors:
           if(corrsec > 0.0):
              print wdate,wdate.timestamp
           wdate = wdate + corrsec
           if(corrsec > 0.0):
              print wdate,wdate.timestamp

           #Check if event is within box:
           if(lat>=latmin) and (lat<=latmax) and (lon>=lonmin) and (lon<=lonmax):
              #BUG! -> This still contains potentially wrong minutes + seconds... get mi, ss from wdate instead (because that contains the corrected version)
              #evlist.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,rphIDX+1,-9,-9,0.0,-9,1))
              evlist.append(hypo(wdate.year,wdate.month,wdate.day,wdate.hour,wdate.minute,float(wdate.second)+float(wdate.microsecond/1000000.0),wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,-9,-9,-9,0.0,-9,1))

           else:
              #Blank ID1 will prevent the phases appended to list
              print '--> WARNING: Event outside of defined region, skip event:',ID1,lat,lon,mag
              ID1=''

           #next line:
           continue

        else:
           #next line:
           continue

    return evlist
#------------------------------------------------------------------------------------------
def loadCNV(ifile,catID,dbgmode,latmin,latmax,lonmin,lonmax,weights,hypo4replace):

    evlist = []
    phlist = []
    revIDX = -1
    rphIDX = -1  

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
           wsu=-9.0     
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

           #Check for formating problems in absolute time:
           corrsec = 0.0

           #Check second:
           if(ss >= 60.0):
              ss = ss - 60.0
              corrsec = corrsec + 60.0
              print '    Second >= 60.0 --> substract 60 s from ss and add 60 s          later for event:',ID1,corrsec
              #ss = ss - 0.00001

           #Check minute:
           if(mi >= 60.0):
              mi = mi - 60.0
              corrsec = corrsec + (60.0*60.0)
              print '    Minute >= 60.0 --> substract 60 m from mi and add 60 m (3600 s) later for event:',ID1,corrsec

           #print "%04d-%02d-%02d %02d:%02d:%09.6f %s" % (yy,mm,dd,hh,mi,ss,ID1)
           wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.999995 -> 60.0 Fix: increase precision

           #Now add possible corrections due to formating errors:
           if(corrsec > 0.0):
              print wdate,wdate.timestamp
           wdate = wdate + corrsec
           if(corrsec > 0.0):
              print wdate,wdate.timestamp

           #Check if alternative solutions exists:
           if(len(hypo4replace)>0):
              altID=getevntidx(hypo4replace,ID1)
           else:
              altID=-9  

           #Check if event should be considered:
           if(altID == -1):
              #Alternative solutions exist but this event is not included...
              #Blank ID1 will prevent the phases appended to list
              print '--> WARNING: Event not found in list of hypocenters for replacement:',ID1,'-> skip event'
              ID1=''

              #next line:
              continue

           #Check if hypocenter information should be replaced:
           #yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,wgsum,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,stphaidx,enphaidx,matchidx,evscore,gid,useflg
           if(altID >= 0):
              Fyy = hypo4replace[altID].yy
              Fmm = hypo4replace[altID].mm
              Fdd = hypo4replace[altID].dd
              Fhh = hypo4replace[altID].hh
              Fmi = hypo4replace[altID].mi
              Fss = hypo4replace[altID].ss
              Fts = hypo4replace[altID].timestamp
              FOTcorr = wdate.timestamp - hypo4replace[altID].timestamp #Important to correct traveltimes later...
              lon = hypo4replace[altID].lon
              lat = hypo4replace[altID].lat
              dep = hypo4replace[altID].dep
              rms = hypo4replace[altID].rms
              mag = hypo4replace[altID].mag
              gap = hypo4replace[altID].gap
              #Update swiss coordinates:
              chx,chy = celleb(lon,lat)
              if(chy < 62.0) or (chy > 302.0) or (chx < 480.0) or (chx > 847.5):
                 chy = 999
                 chx = 999
           else:
              Fyy = wdate.year
              Fmm = wdate.month
              Fdd = wdate.day
              Fhh = wdate.hour
              Fmi = wdate.minute
              Fss = float(wdate.second)+float(wdate.microsecond/1000000.0)
              Fts = wdate.timestamp
              FOTcorr = 0.0

           #Check if event is within box:
           if(lat>=latmin) and (lat<=latmax) and (lon>=lonmin) and (lon<=lonmax):
              #First check if phase-position of previous event needs to be stored in hypo-list:
              if(rphIDX >= 0) and (revIDX >= 0):
                 evlist[revIDX].enphaidx = rphIDX

              #Now new event starts
              revIDX += 1
              #BUG! -> This still contains potentially wrong minutes + seconds... get mi, ss from wdate instead (because that contains the corrected version)
              #evlist.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,rphIDX+1,-9,-9,0.0,-9,1))

              #Former working version
              #evlist.append(hypo(wdate.year,wdate.month,wdate.day,wdate.hour,wdate.minute,float(wdate.second)+float(wdate.microsecond/1000000.0),wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,rphIDX+1,-9,-9,0.0,-9,1))

              #New version allowing for replaced hypocenters:
              evlist.append(hypo(Fyy,Fmm,Fdd,Fhh,Fmi,Fss,Fts,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,rphIDX+1,-9,-9,0.0,-9,1))

              #if(dbgmode == 0):
                  #print hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,rphIDX+1,-9,-9,0.0,-9,1)
           else:
              #Blank ID1 will prevent the phases appended to list
              print '--> WARNING: Event outside of defined region, skip event:',ID1,lat,lon,mag
              ID1=''           

           #next line:
           continue

        else:

           #Check if line containes one phase (stat,ptype,qual,err,weig,tt,timestamp,edist,azi,baz,evID): weights[int(line[5:6])].err,weights[int(line[5:6])].weig
           if(len(ID1)>0) and (len(line.rstrip())==12):
              #print float(line[6:12]),FOTcorr,float(line[6:12])+FOTcorr
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),weights[int(line[5:6])].err,weights[int(line[5:6])].weig,float(line[6:12])+FOTcorr,Fts+float(line[6:12])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1)) 
              #Update phase index counter:
              rphIDX += 1

           #Check if line containes two phase:
           if(len(ID1)>0) and (len(line.rstrip())==24):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),weights[int(line[5:6])].err,weights[int(line[5:6])].weig,float(line[6:12])+FOTcorr,Fts+float(line[6:12])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),weights[int(line[17:18])].err,weights[int(line[17:18])].weig,float(line[18:24])+FOTcorr,Fts+float(line[18:24])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              #Update phase index counter:
              rphIDX += 2

           #Check if line containes three phase:
           if(len(ID1)>0) and (len(line.rstrip())==36):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),weights[int(line[5:6])].err,weights[int(line[5:6])].weig,float(line[6:12])+FOTcorr,Fts+float(line[6:12])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),weights[int(line[17:18])].err,weights[int(line[17:18])].weig,float(line[18:24])+FOTcorr,Fts+float(line[18:24])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),weights[int(line[29:30])].err,weights[int(line[29:30])].weig,float(line[30:36])+FOTcorr,Fts+float(line[30:36])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              #Update phase index counter:
              rphIDX += 3

           #Check if line containes four phase:
           if(len(ID1)>0) and (len(line.rstrip())==48):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),weights[int(line[5:6])].err,weights[int(line[5:6])].weig,float(line[6:12])+FOTcorr,Fts+float(line[6:12])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),weights[int(line[17:18])].err,weights[int(line[17:18])].weig,float(line[18:24])+FOTcorr,Fts+float(line[18:24])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),weights[int(line[29:30])].err,weights[int(line[29:30])].weig,float(line[30:36])+FOTcorr,Fts+float(line[30:36])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[36:40],line[40:41],int(line[41:42]),weights[int(line[41:42])].err,weights[int(line[41:42])].weig,float(line[42:48])+FOTcorr,Fts+float(line[42:48])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              #Update phase index counter:
              rphIDX += 4 

           #Check if line containes five phase:
           if(len(ID1)>0) and (len(line.rstrip())==60):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),weights[int(line[5:6])].err,weights[int(line[5:6])].weig,float(line[6:12])+FOTcorr,Fts+float(line[6:12])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),weights[int(line[17:18])].err,weights[int(line[17:18])].weig,float(line[18:24])+FOTcorr,Fts+float(line[18:24])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),weights[int(line[29:30])].err,weights[int(line[29:30])].weig,float(line[30:36])+FOTcorr,Fts+float(line[30:36])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[36:40],line[40:41],int(line[41:42]),weights[int(line[41:42])].err,weights[int(line[41:42])].weig,float(line[42:48])+FOTcorr,Fts+float(line[42:48])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[48:52],line[52:53],int(line[53:54]),weights[int(line[53:54])].err,weights[int(line[53:54])].weig,float(line[54:60])+FOTcorr,Fts+float(line[54:60])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              #Update phase index counter:
              rphIDX += 5

           #Check if line containes six phase:
           if(len(ID1)>0) and (len(line.rstrip())==72):
              phlist.append(phase(line[0:4],line[4:5],int(line[5:6]),weights[int(line[5:6])].err,weights[int(line[5:6])].weig,float(line[6:12])+FOTcorr,Fts+float(line[6:12])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[12:16],line[16:17],int(line[17:18]),weights[int(line[17:18])].err,weights[int(line[17:18])].weig,float(line[18:24])+FOTcorr,Fts+float(line[18:24])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[24:28],line[28:29],int(line[29:30]),weights[int(line[29:30])].err,weights[int(line[29:30])].weig,float(line[30:36])+FOTcorr,Fts+float(line[30:36])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[36:40],line[40:41],int(line[41:42]),weights[int(line[41:42])].err,weights[int(line[41:42])].weig,float(line[42:48])+FOTcorr,Fts+float(line[42:48])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[48:52],line[52:53],int(line[53:54]),weights[int(line[53:54])].err,weights[int(line[53:54])].weig,float(line[54:60])+FOTcorr,Fts+float(line[54:60])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              phlist.append(phase(line[60:64],line[64:65],int(line[65:66]),weights[int(line[65:66])].err,weights[int(line[65:66])].weig,float(line[66:72])+FOTcorr,Fts+float(line[66:72])+FOTcorr,-9.0,-9.0,-9.0,ID1,-99.0,-999.0,-1,0.0,-1))
              #Update phase index counter:
              rphIDX += 6

    #Close
    fp.close()

    #Finally add last phase-position of previous event to hypo-list:
    evlist[revIDX].enphaidx = rphIDX

    #Finally, index the entire phase list:
    for l in range(len(phlist)):
        phlist[l].phaidx = l

    return evlist,phlist

#------------------------------------------------------------------------------------------
def loadIDlist(ifile):

    IDlist = []

    #Open
    fp = open(ifile,"r")

    for line in fp:
        #Test if line is empty:
        if(len(line.split())<1):
           print "Line not considered:", line
        else:
           #Test if line is comment:
           #t=line.split()
           t = string.split(line,None,-1)
           s = string.split(line,"|",-1)

           if(t[0].startswith("#") == False):
              IDlist.append(t[0])

    #Close
    fp.close()

    return IDlist
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
                 wsu = -9

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
                 wsu = -9


              #Check if event was extracted from SC3 DB structure:
              #STILL MISSING 
              #else:

              #Check if format is correct:
              if(ss < 60.0) and (mi < 60) and (hh < 24) and (dd < 32) and (mm < 13):
                 wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss)) #Problems with seconds close to 60.0 e.g.: 59.99995 -> 60.0 Fix: increase precision
                 #wdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss))

                 #print ID1
                 list.append(hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,rphIDX+1,-9,-9,0.0,-9,1))
                 if(dbgmode == 1):
                    print hypo(yy,mm,dd,hh,mi,ss,wdate.timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mat,mnobs,merr,magmeth,rms,gap,mdi,nob,wsu,ety,lqa,chx,chy,agy,ID1,ID2,ID4,ID3,loc,modlID,author,region,fstat,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,rphIDX+1,-9,-9,0.0,-9,1)
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
def hypoFilter(hypos,quakes2delete):
    filteredlist = []

    print ""
    print "--> Events to delete loaded from file quakes2delete"

    #Loop over all events:
    for i in range(len(hypos)):

        todelete = 0

        #Check if event is listed to be removed:
        for j in range(len(quakes2delete)):
            if(hypos[i].evID1 == quakes2delete[j]):
               todelete = 1
               break

        #Append event:
        if(todelete == 0):
           filteredlist.append(hypos[i])
        else:
           print "--> Event marked to be deleted, remove event",hypos[i].evID1

    return filteredlist
#------------------------------------------------------------------------------------------

###########################################################################################
# Main code:
###########################################################################################

#Loads a phase file, calculates all geometrical location parameters from scratch and filters the phase
#data based on VELEST output, some basic filters are implemented:
#1) GAP cut-off
#2) RMS cut-off
#3) Station-Residual cut-off (e.g. if residual is X times larger than picking error)
#4) Remove events specified in a list
#5) Remove phases specified in a list

#Recalculate all parameters after filtering to make sure quality is still sufficient

#Now from command line:
##################################################################################################################
#Get command line arguments:
oparser = optparse.OptionParser()
oparser.add_option('-F', '--phaseFile1',  action='store', dest='PhaseFile1' ,  help='Phase file to filter (default: /Users/tdiehl/data/LOTO/SWISS_TOMO/CH_Update/VELEST/Data3/simul_01_SHA01/SH_A_001.cnv)')
oparser.add_option('-I', '--phaseForm1',  action='store', dest='PhaseForm1' ,  help='Format of Phase file (default: cnv)')
oparser.add_option('-o', '--output',      action='store', dest='OutputFile',   help='Filename of filtered output (default: phase_filter4tomo.out )') 
oparser.add_option('-e', '--exportformat',action='store', dest='OutputForm',   help='Format of merged output: cnv, simulps, MANUPICK, MANUPICK_alias2org (replaces station names in input with original name, use for NLL) (default: cnv)')
oparser.add_option('-s', '--statfile',    action='store', dest='StationFile',  help='Alias station file (default: /Users/tdiehl/lib/sed_stations.GSE_SED.alias)')
oparser.add_option('-v', '--velestoutputroot',action='store', dest='VelestRoot', help='Root of VELEST output files including additional information on minimum 1D-statistics (station-corrections, etc., default: SH_A_001)')
oparser.add_option('-r', '--residualoutlier',action='store', dest='outlierfile', help='File with residuals identified as outliers, e.g. output of filter4tomo_Waveform.py, default: None')
oparser.add_option('-q', '--quakes2delete',action='store', dest='quakes2deletefile', help='File with earthquakes to delete (one eventID per line), default: None')
oparser.add_option('-c', '--checkquakes',action='store', dest='checkquakesfile', help='File with earthquakes to check (one eventID per line), markes all phases of specified events as potential outliers, default: None')
oparser.add_option('-H', '--hypocenter4replacement',action='store', dest='hypo4replace', help='Hypocenters in specified cnv file will replace the original ones, events not included in this file will be removed, no further filtering, default: None')
oparser.add_option('-C', '--ConsistencyCheckResidual',action='store_true', dest='RCCflg', help='If set, consitency check of |residuals| > threshold is performed (default: false); station2source+event2receiver')
oparser.add_option('-d', '--DistanceSelection', action='store', dest='DistSelFlg', help='Distance selection options; All: no restriction applied (default); Pg_Dyn: Dynamic distance range for Pg dominated selections (distance range defined by parameters maxepidist1 and maxepidist2 harcoded in program (e.g. 80-140 km); Pg_Dyn_Reloc: same as Pg_Dyn but with loose event quality criteria to accept almost every event for relocation')
oparser.add_option('-E', '--ExportOnly',action='store_true', dest='Exportflg', help='If set, no filters will be applied, only export will be performed, default: False -- NOT FULLY OPERATIONAL/TESTED YET')
oparser.add_option('-b', '--boxrange',           action='store', dest='boxrange',    help='This parameter defines a box for event selection by latmin/latmax/lonmin/lonmax; default: 45.0/48.5/5.25/11.5')
(options, arg) = oparser.parse_args()

#Box definition:
if(options.boxrange != None):
   boxrange = options.boxrange
else:
   boxrange = '45.0/48.5/5.25/11.5'
sp = string.split(boxrange,'/',-1)              #Split string:
latmin = float(sp[0])
latmax = float(sp[1])
lonmin = float(sp[2])
lonmax = float(sp[3])

#Check if file with quakes to check exists:
if(options.DistSelFlg != None):
   DistSelFlg = options.DistSelFlg
else:
   DistSelFlg = 'All'
if(DistSelFlg!='All') and (DistSelFlg!='Pg_Dyn') and (DistSelFlg!='Pg_Dyn_Reloc'):
   print "  -> WARNING: Unknown Distance Selection Option",DistSelFlg,'use All instead'
   DistSelFlg = 'All'

#Check if Export-Only is requested:
if(options.Exportflg != None):
   Exportflg = True
else:
   Exportflg = False

#Check if consistency check is performed:
if(options.RCCflg != None):
   RCCflg = True
else:
   RCCflg = False

#Check if file with hypocenters for replacement exists:
if(options.hypo4replace != None):
   hypo4replacefile = options.hypo4replace
   rphypos = 1
else:
   hypo4replacefile = 'None'
   rphypos = 0

#Check if file with quakes to check exists:
if(options.checkquakesfile != None):
   checkquakesfile = options.checkquakesfile
   ckQuakes = 1
else:
   checkquakesfile = 'None'
   ckQuakes = 0

#Check if file with quakes to delete exists:
if(options.quakes2deletefile != None):
   quakes2deletefile = options.quakes2deletefile
   rmQuakes = 1
else:
   quakes2deletefile = 'None'
   rmQuakes = 0

#Check if root of outlierfile was specified:
if(options.outlierfile != None):
   outlierfile = options.outlierfile
   rmOutliers = 1
else:
   outlierfile = 'None'
   rmOutliers = 0

#Check if root of VELEST output file was specified:
if(options.VelestRoot != None):
   VelestRoot = options.VelestRoot
else:
   VelestRoot = 'SH_A_001' 

#Check if PhaseFile 1 is included in list of input arguments:
if(options.PhaseFile1 != None):
   PhaseFile1 = options.PhaseFile1
else:
    #PhaseFile1 = '/Users/tdiehl/data/LOTO/SWISS_TOMO/CH_Update/VELEST/Data3/simul_01_SHA01/SH_A_001.cnv'
    PhaseFile1 = '/Users/alejandro/ICELAND_Project/Tomography/SIMULPS14/Tools/python/P+S/CS_invertratio1_PS.cnv' #ADN

#Check if PhaseForm1 is included in list of input arguments:
if(options.PhaseForm1 != None):
   PhaseForm1 = options.PhaseForm1
else:
   PhaseForm1 = 'cnv'

#Check if OutputFile is included in list of input arguments:
if(options.OutputFile != None):
   OutputFile = options.OutputFile
else:
   if(OutputForm != 'MANUPICK') and (OutputForm != 'MANUPICK_alias2org'):
      OutputFile = 'phase_filter4tomo.out'
   else:
      OutputFile = 'events'

#Check if OutputForm is included in list of input arguments:
if(options.OutputForm != None):
   OutputForm = options.OutputForm
else:
   OutputForm = 'cnv'

#Check if StationFile is included in list of input arguments:
if(options.StationFile != None):
   StationFile = options.StationFile
else:
    #StationFile = '/Users/tdiehl/lib/sed_stations.GSE_SED.alias'
   #StationFile = '/Users/tdiehl/data/LOTO/SWISS_TOMO/CH_Update/VELEST/Data3/simul_01_SHA01/SH_A_001.sta'
   StationFile = '/Users/alejandro/ICELAND_Project/Tomography/SIMULPS14/Tools/python/P+S/CS_invertratio1_PS.sta' #ADN


#Hardcoded-log file
LogFileFilter = 'phase_filter4tomo.filter.log'
LogFileEvents = 'phase_filter4tomo.events.log'
LogFilePhases = 'phase_filter4tomo.phases.log'
LogFileStations = 'phase_filter4tomo.stations.log'
LogFileStatCorr = 'phase_filter4tomo.StatCorrPosOutlier.log'
if(os.path.isfile(outlierfile)) and (rmOutliers == 1):
   LogFilePhas2Del = 'phase_filter4tomo.Outliers2Delete.log'
else:
   LogFilePhas2Del = 'phase_filter4tomo.Phases2Delete.log'
LogFileEven2Del = 'phase_filter4tomo.Events2Delete.log'
LogFileRCC = 'phase_filter4tomo.RCC.log'

#Setup velest statistics output file names:
VelestOut_StatCor_AlSta = VelestRoot + '.sta'
VelestOut_StatCor_Stati = VelestRoot + '.statistic.sta'
VelestOut_StatCor_AzRes = VelestRoot + '.azres'
ResidualFile = VelestRoot + '.stares'

#Hardcoded applications:
#Distance selection (All: no distance restrictions; Pg_Dyn: Dynamic selection):
#NOW: Defined as command-line option 
#DistSelFlg = 'All'
#DistSelFlg = 'Pg_Dyn'

#Maximum epicentral distance (dynamic definition for Pg tomo):
#Use all phases upto distances of maxepidist1. Use distances
#maxepidist1 <= dist <= only until gap is <180 and nobs>=8
#For Pg+Pn (All, default):
maxepidist1 = 999999.9
maxepidist2 = 999999.9
maxepidist1S = 999999.9      #Maximum epicentral distance to include S-phases: HQ
maxepidist2S = 999999.9      #Maximum epicentral distance to include S-phases: Maximum
#for Pg dominated tomo:
if(DistSelFlg=='Pg_Dyn'):
   maxepidist1 = 80.0
   maxepidist2 = 140.0
   maxepidist1S = 80.0
   maxepidist2S = 140.0
if(DistSelFlg=='Pg_Dyn_Reloc'):
   maxepidist1 = 80.0
   maxepidist2 = 140.0
   maxepidist1S = 25.0
   maxepidist2S = 60.0
print '  -> Distance selction flag',DistSelFlg,'Set maxepidist1 to maxepidist2 to:',maxepidist1,maxepidist2,' | Maximum distance-range to use S-phases:',maxepidist1S,maxepidist2S

#For dynamic distance cut-off (e.g. Pg-tomo):
#Default values for tomography (Pg_Dyn)
#Minimum number of observations to finish phase compilation at distances > maxepidist1:
dyndist_minnob = 10
dyndist_minnobS = 5   #Stop including S-phases if distance > maxepidist1S and dyndist_minnobS - observations (P+S) are available
#Maximum gap to finish phase compilation at distances > maxepidist1:
dyndist_maxgap = 170  #Used for tomo - maybe we should have used a lower value here...
dyndist_maxgapS = 360 #Stop including S-phases if distance > maxepidist1S and dyndist_maxgapS is reached

#Default values for relocation studies:
if(DistSelFlg=='Pg_Dyn_Reloc'):
   dyndist_minnob = 10
   dyndist_maxgap = 130  #Used for relocation CaseStudies

#Cut-off values for quality criteria:
#Pure removal of phases (often done with initial CNV file) -> Skip test of gap:
if((os.path.isfile(outlierfile)) and (rmOutliers == 1)) or (rmQuakes == 1) or (ckQuakes == 1) or (rphypos == 1) or (Exportflg and ((OutputForm == "MANUPICK") or (OutputForm == "MANUPICK_alias2org"))):
   co_maxGap = 360.0
else:
   #co_maxGap = 180.0    #For VELEST filtering
   co_maxGap = 183.0    #Convert CNV to final EQKS
   #co_maxGap = 190.0   #For Sg select same criteria as for Pg

#Minimum number of observation per event:
#co_minNob = 8 #Conservative
if(Exportflg and ((OutputForm == "MANUPICK") or (OutputForm == "MANUPICK_alias2org"))):
   co_minNob = 4 #Minimum for initial selection - For export to MANUPICK take everything!
else:
   co_minNob = 6 #Minimum for initial selection
if(DistSelFlg=='Pg_Dyn'):
   co_minNob = 6      #For Pg tomo -> To make sure most events will have a Pg location...
   co_maxGap = 190.0  #Maybe Pg - location changes and gap becomes <180 after joint inversion...
if(DistSelFlg=='Pg_Dyn_Reloc'):
   co_minNob = 4      #For relocation CaseStudies
   co_maxGap = 360.0  #For relocation CaseStudies
print '  -> Minimum number of observation per event set to:',co_minNob
print '  -> Maximum GAP    of event                 set to:',co_maxGap

#Hardcoded filter parameters:
#Maximum event RMS, if larger remove (Test #1) Proposal of EK: threshold = 2 * Total RMS residual of inversion 
#should be not too strict before large residuals haven't been removed
#Pure removal of phases (often done with initial CNV file) -> Skip test of RMS:
if(os.path.isfile(outlierfile)) and (rmOutliers == 1):
   #co_maxRMS = 2.0  #For initial selections 
   co_maxRMS = 99999999.9  #After Minimum 1D, e.g. if convert to SIMULPS - Alpine 1D (unfiltered): RMS=0.3 --> 0.6 corresponds to 2*RMS
else:
   co_maxRMS = 7.0  #For initial selections 
   #co_maxRMS = 0.6  #After Minimum 1D, e.g. if convert to SIMULPS - Alpine 1D (unfiltered): RMS=0.3 --> 0.6 corresponds to 2*RMS
print '  -> Maximum event-RMS set to                      :',co_maxRMS

#Maximum difference in average station reidual between azimuthal quadrants (Test #2):
arpQthr = 0.4

#Maximum phase/station residual, if negative, it is calculated from picking error: maxRES = abs(co_maxRES)*pickuncertaintu
#                                if negative, an outlier is defined if AverageStationResidual - IndividualResidual >= IndividualObservationError
#                                if positive, it represents the absolute cut-off value
#co_maxRES = -2.0 
#co_maxRES = 0.6
co_maxRES = 0.4
print '  -> Maximum phase residual threshold              :',co_maxRES

# Perform outlier statistic test: Residual Consitency Check (RCC)
# For potential outliers (|RES|>co_maxRES), check:
if(RCCflg):
   RCC_Check_Flg = 1
   print ""
   print "  -> Perform Residual Consistency Check..."
else:
   RCC_Check_Flg = 0
   print ""
   print "  -> Skip    Residual Consistency Check..."
# 1) Station-check: check all residuals from the same source region for the given station.
#    check if residual is consistent with the rest of the residuals from this source region.
#    if significantly smaller residuals exists fro this source region, the residual is likely
#    an mispick or related to GPS timing issues... If similar residuals exist, make
#    sure that they occure at a different time (e.g. seperated by 4 months) to be sure
#    that it isn't related to the same GPS/Timing problem...
# Search radius in source region (epicenter, km):
RCC_Stat_rad_Epi = 10.0  #5.0
# Search radius in source region (depth, km):
RCC_Stat_rad_Dep = 10.0  #7.5

# 2) Event check: check residuals of SAME event and similar receiver-region:
#  Search sector in receiver region (epicenter, km):
RCC_Even_rad_Sta = 10.0  #+/- Distance in km to consider stations
RCC_Even_bin_Azi = 10.0  #+/- Azimuth in deg to consider stations 

# Parameters for automatic decision:
#Use dynamic tolerance based on standard deviation and mean (False/True):
RCC_ResDynamicTolerance = True
# Static tolerance in sec (use +/- average residual of minimum 1D inversion), also used if less than 2 obervations for std:
RCC_ResStaticTolerance = 0.1
# Temporal gap between two residuals in months (to make sure possible GPS issues were fixed)
RCC_TemporalGap = 6
#Minimum number of matches for statistical-criteria:
RCC_Min_RSR = 1  #Similar source  
RCC_Min_RRR = 2  #Similar receiver (at least two stations are required?)
#Minimum number of matches for count-criteria:
RCC_MinCount = 3

#Minimum obervation class (qualities > co_minQual will be removed):
co_minQual = 4

#Maximum outlier-score to keep a phase, if outlier-score > this value for a phase, phase will be deleted:
oscoreMAX = 0.9

#Lat/Lon range of events to include  45.0/48.5/5.25/11.5
#Now from input argument...
#latmin = 45.0
#latmax = 48.5
#lonmin = 5.25
#lonmax = 11.5

#---ADN 12.03.2019---
#Limits COSEISMIQ Region
latmin = 60
latmax = 70
lonmin = -25
lonmax = -15

##################################################################################################################

#Debug-Mode (0 or 1):
dbgmode = 0

#Check if export-only flag was set:
if(Exportflg):
   print "--> INFO: Export only flag is set, do not apply any filter"
   #Set all other flags to 0:
   rmQuakes = 0
   ckQuakes = 0
   rmOutliers = 0

fplog = open(LogFileFilter,'w')
fplogEv = open(LogFileEvents,'w')
fplogPh = open(LogFilePhases,'w')
fplogSt = open(LogFileStations,'w')
fplogSC = open(LogFileStatCorr,'w')
fplogP2D = open(LogFilePhas2Del,'w')
fplogE2D = open(LogFileEven2Del,'w')

if(RCC_Check_Flg==1):
   fplogRCC = open(LogFileRCC,'w')

#Setup weighting schemes for file:
wscheme1 = setupWeighting_SED()

#Load station file:
print ""
print "--> Load station file from:",StationFile
stations = loadVelestStations(StationFile)
#stations = loadAliasStations(StationFile)   #ADN
print '    Number stations:',len(stations)

#Load hypocenters for replacement:
if(rphypos == 1):
   print ""
   print "--> Load hypocenters for replacement from CNV-file: ",hypo4replacefile
   hypo4replace = loadCNVhypos(hypo4replacefile,2,latmin,latmax,lonmin,lonmax)
   print "--> Number of hypocenters loaded for replacement:",len(hypo4replace)
else:
   hypo4replace = []
   
#Load phase file #1:
print ""
if(PhaseForm1 == 'cnv'):
   if(rphypos == 1):
      print "--> Load Phases #1 from CNV-file: ",PhaseFile1,"for which a location is provided in",hypo4replacefile
   else:
      print "--> Load Phases #1 from CNV-file: ",PhaseFile1
   hypos1,phases1 = loadCNV(PhaseFile1,1,dbgmode,latmin,latmax,lonmin,lonmax,wscheme1,hypo4replace)

print ""
print '    Number events phasefile1:',len(hypos1)
print '    Number phases phasefile1:',len(phases1)
print ""
org_no_phases = len(phases1)
org_no_events = len(hypos1)

#Load events to check:
if(ckQuakes == 1):
   print "--> Load events to check  from  : ",checkquakesfile
   checkquakes = loadIDlist(checkquakesfile)
else:
   checkquakes = []

#Load potential events to delete:
if(rmQuakes==1):
   print "--> Load events to delete from  : ",quakes2deletefile
   quakes2delete = loadIDlist(quakes2deletefile)
else:
   quakes2delete = []

#Load station corrections:
print ""
print "--> Load statcor      from:",VelestOut_StatCor_AlSta
stations = loadstatcor(stations,VelestOut_StatCor_AlSta)

#Add statistical information (STD) for station corrections
print ""
print "--> Load statcor-stat from:",VelestOut_StatCor_Stati
stations = loadstatcorstati(stations,VelestOut_StatCor_Stati)

#Add azimuthal information Q1-Q4 for station corrections
print ""
print "--> Load statcor-azim from:",VelestOut_StatCor_AzRes
stations = loadstatcorazi(stations,VelestOut_StatCor_AzRes)

#Compare average residuals per quadrants and check if difference is > arpQthr
stations2check = check_ST_avresperquadrant(stations,arpQthr,'P',fplogSC) 

#Sort both hypocenter lists with respect to origin time:
hypos1.sort(key=lambda x: (x.timestamp), reverse=False)

#Now copy the hypocenter list #1 as the merged hypocenter list:
#If you want to delete entire events specificed in a list, you could do it here. Only
#keep events not filtered out...
if(rmQuakes==1) and (len(quakes2delete)>0):
   filteredHypo = hypoFilter(hypos1,quakes2delete)
else:
   filteredHypo = hypos1

#Now sort the merged catalog (with respect to origin time):
#filteredHypo.sort(key=lambda x: (x.timestamp), reverse=False)

#Counters:
no_removed_events = 0
no_removed_phases = 0

#Write some header:
fplogP2D.write("#Action    Stat Pha Q   Error Pick-Time                     resid | SCNob StatCor AvResid ResiSTD AvResQa NobQa Qa | EpiDist  Azim   BAZ   statlon  statlat  AvDist STDDist | Hypo-Origin-Time                evlon    evlat    evdep   mag StaOrg O-SCORE WFRew EV-ID\n")
fplogE2D.write("#Action    Reason       RMS    GAP   Nobs Hypo-Origin-Time                evlon    evlat    evdep   mag EV-ID\n")

######################################################################################################################################
#Loop #1: Loop ove all phases of selected phases to calculate distance and azimuth/baz of each phase (also adds additional information to phase-fields):
#         This loop compiles information also about station statistics, such as average distance etc.

#Check if first first cut-off should be applied:
#Should not be applied if only phases removed, events checked, events filtered, hypocenters updated:
if(rmQuakes == 0) and (rmOutliers == 0) and (ckQuakes == 0) and (rphypos == 0):
   maxepidistTMP = maxepidist2
   maxepidistTMPS = maxepidist2S
else:
   maxepidistTMP = 9999999.9
   maxepidistTMPS = 9999999.9

print ""
print "--> calculate distance & azimuth & baz for each phase and ignore all phases at distances >",maxepidistTMP,"km and S-phases >",maxepidistTMPS,"km"

#Vector containing information on distances for phases
DistVect = []

#Event loop:
for i in range(len(filteredHypo)):

    #Phase loop:
    for j in range(filteredHypo[i].stphaidx, filteredHypo[i].enphaidx+1):

        #Adjust travel time to correct origin time in filtered catalog
        phases1[j].tt = phases1[j].timestamp - filteredHypo[i].timestamp

        #Final check of event IDs:
        if(phases1[j].evID == filteredHypo[i].evID1):

           #ADDITIONAL PHASE FILTERS COULD BE ADDED HERE LATER (residuals, quality, Auto-Manual, etc...

           #Get distance and azimuth of each phase (also adds information to phase-fields):
           #In this function we could add a distance filter to seperate Pg from Pn...
           #We use maxepidistTMP as absolute maximum for distance
           phases1[j],DistVect = getDistAzi4phase(phases1[j],filteredHypo[i].lat,filteredHypo[i].lon,stations,maxepidistTMP,maxepidistTMPS,DistVect)

######################################################################################################################################
#Loop #2: Loop over stations to calculate average distance/std and average BAZ/std:

print ""
print "--> calculate mean distance & standard-deviation for each station"

#Station loop:
for i in range(len(stations)):

    nobs = 0

    #Calculate statistsics for P-wave corrections:
    if(stations[i].SC_P_NobAll > 0):

       #Calculate average:
       stations[i].SC_P_DISAve = stations[i].SC_P_DISAve/stations[i].SC_P_NobAll

       #Standard deviation:
       for j in range(len(DistVect)):
           if(DistVect[j].stat == stations[i].stat) and (DistVect[j].phas == 'P'):
              nobs += 1
              stations[i].SC_P_DISSTD =  stations[i].SC_P_DISSTD + ((DistVect[j].dist - stations[i].SC_P_DISAve)**2)

       if(nobs == stations[i].SC_P_NobAll) and (nobs > 1):
          stations[i].SC_P_DISSTD = math.sqrt(stations[i].SC_P_DISSTD/(nobs-1))
       else:
          print "--> WARNING: Cannot calculate the standard deviation:",stations[i].stat,stations[i].SC_P_NobAll,nobs
          stations[i].SC_P_DISSTD = -1
    
######################################################################################################################################
#Read list of outliers if existing and requested, otherwise read residuals and decide what is outlier
if(os.path.isfile(outlierfile)) and (rmOutliers == 1):
   print ""
   print "--> Load outliers     from:",outlierfile 
   phase2delete,reslist = loadOutliers(outlierfile,oscoreMAX,wscheme1)
   print '    Number of phases to remove as defined in outlierlist             :',len(reslist)
else:
   #If file not existing switch-off rmOutliers and load residuals and select potential outliers based on residual, and station-correction statistics...
   rmOutliers = 0

   #Load residual lists only if rphypos == 0 (in case only hypocenters will be updated, residual check is not necessary):
   if(rphypos == 1) or (rmQuakes == 1) or (Exportflg):
      reslist = []
      reslist_All = []
      phase2delete = []
   else:
      #Get all residuals for certain events or identify outliers:
      if(ckQuakes == 0):
         print ""
         print "--> Load residuals    from:",ResidualFile
         reslist,reslist_All = loadRES(ResidualFile,co_maxRES,wscheme1,stations,hypos1)
         #reslist    : Contains all potentential residuals exceeding thresholds
         #reslist_All: Contains ALL residuals for later statistsics
         print '    Number of phases to remove due to residuals larger than threshold:',len(reslist)
         print '    Number of phases in complete residual list (ALL phases)          :',len(reslist_All)        
         phase2delete = []   #Empty list which contains possible outliers
      else:
         print ""
         print "--> Load residuals    from:",ResidualFile,"for all events to be checked in file",checkquakesfile
         reslist = loadRES4check(ResidualFile,checkquakes,wscheme1,stations,hypos1)
         print '    Number of phases to check for events (and remove in this run)    :',len(reslist)
         phase2delete = []   #Empty list which contains possible outliers
         reslist_All = []    #Dummy empty list, just to make sure it's not missing later

######################################################################################################################################
#Loop #3: Loop over all events: collect information of potential outlier phases and mark corresponding phases in phases-list
#         -> Link Residual vector with phaselist

print ""
print "--> Associate residuals with phases in phase-list and link with statistical information: -> List with potential outliers to delete"

#Loop over all selected events:
for i in range(len(filteredHypo)):

    #Loop over all phases associated with this event:
    for j in range(filteredHypo[i].stphaidx, filteredHypo[i].enphaidx+1):

        #Final check of event IDs and if phase is still used:
        if(phases1[j].evID == filteredHypo[i].evID1) and (phases1[j].qual<wscheme1[len(wscheme1)-1].qual) and (phases1[j].weig>0.00001) and (phases1[j].oscore <= oscoreMAX):

           #Check if phase is listed in residuals2remove:
           #stat,ptype,qual,err,weig,res,dis,evID
           outlierflg,resid,residx,SCasrQ,SCnsrQ,Quad = check4residual(phases1[j].evID,phases1[j].stat,phases1[j].ptype,phases1[j].qual,reslist)

           #Check if phase is listed in residual vector: if outlierflg != 0 phase is listed...
           if(outlierflg != 0):

              #Default outlier-score:
              oscore = 1.0
              phases1[j].oscore = oscore

              #Store phase as potential outlier (not if outliers are read already):
              #stat,statorg,pha,qua,err,picktimestamp,res,SCnob,SCval,SCasr,SCssr,SCasrQ,SCnsrQ,Quad,epi,azi,baz,slon,slat,asepi,ssepi,OTtimestamp,elon,elat,edep,emag,oscore,evID,phaidx
              if(rmOutliers == 0):
                 phase2delete.append(outlier(phases1[j].stat,stations[phases1[j].stidx].statorg,phases1[j].ptype,phases1[j].qual,phases1[j].err,phases1[j].timestamp,resid,stations[phases1[j].stidx].SC_P_NobAll,stations[phases1[j].stidx].SC_P_DelAll,stations[phases1[j].stidx].SC_P_AREAll,stations[phases1[j].stidx].SC_P_STDAll,SCasrQ,SCnsrQ,Quad,phases1[j].edist,phases1[j].azi,phases1[j].baz,stations[phases1[j].stidx].lon,stations[phases1[j].stidx].lat,stations[phases1[j].stidx].SC_P_DISAve,stations[phases1[j].stidx].SC_P_DISSTD,filteredHypo[i].timestamp,filteredHypo[i].lon,filteredHypo[i].lat,filteredHypo[i].dep,filteredHypo[i].mag,oscore,-1,phases1[j].evID,j))
              else:
                 #Outlier was already stored, in this case only update position of phase:
                 phase2delete[residx].phaidx = j

######################################################################################################################################
#Loop #4: Now go through list of potential outliers and perform additional analysis of outlier candidates:
#         1) Define cut-off value, if |res| > cut-off -> oscore = 1.0 - cut off value e.g. 1.5 s
#         2) Station-wise analysis of outliers: Does same station have similar outliers from same distance/baz range within years (-> GPS problem unlikely, mispick unlikely unless systematically wrong phase)
#         3) Event-wise  analysis of outliers : Does same station have similar outliers at stations in similar distance/azi range (-> GPS problem unlikely, mispick unlikely unless systematically wrong phase)
#         4) Check if residual is positive, more likely mispick if pick is delayed...
#         CHECK/Likelihood  STILL MISSING

#Number of suspicious residuals, which passed the Residual Consistency check
RCC_Check_CNT = 0
RCC_Check_CNT_RRR_Match = 0

if(rmOutliers == 0) and (ckQuakes == 0) and (rphypos == 0) and (rmQuakes == 0):
   print ""
   if(RCC_Check_Flg == 1):
      #Link hypocenter information with ALL residuals:
      print "--> Prepare Consistency check for potential outliers - link hypocenter/station parameters with ALL residuals..."
      #Loop over ALL residuals
      for j in range(len(reslist_All)):

          #Find position of hypocenter:
          for k in range(len(hypos1)):
              if(reslist_All[j].evID == hypos1[k].evID1):
                 reslist_All[j].evidx = k
                 break

          #Find position of station:
          for k in range(len(stations)):
              if(reslist_All[j].stat == stations[k].stat):
                 reslist_All[j].stidx = k
                 break

      print "--> Consistency check for potential outliers: station to source-region AND event to receiver-region..."
      #Loop over all potential outliers:
      print "--> Number of residuals to check:",len(phase2delete)
      for i in range(len(phase2delete)):

          residuals_SourceRegion = []
          residuals_ReceiverRegion = []
          
          #stat,statorg,pha,qua,err,picktimestamp,res,SCnob,SCval,SCasr,SCssr,SCasrQ,SCnsrQ,Quad,epi,azi,baz,slon,slat,asepi,ssepi,OTtimestamp,elon,elat,edep,emag,oscore,WFrevf,evID,phaidx

          #Loop over ALL residuals of the same station, for which hypocenter is close to the actual one OR same event, and receiver close to stations :
          for j in range(len(reslist_All)):

              #Check if hypocenter information exists for residual:
              if(reslist_All[j].evidx < 0):
                 print "--> WARNING: No hypocentral information available for residual of event, ignore residual of event  ",reslist_All[j].evID
                 continue

              #Check if station information exists for residual:
              if(reslist_All[j].stidx < 0):
                 print "--> WARNING: No station     information available for residual of event, ignore residual of station",reslist_All[j].stat,reslist_All[j].evID
                 continue

              #1) Now test if residual matches SAME STATION, phase, and depth:
              if(phase2delete[i].stat == reslist_All[j].stat) and (phase2delete[i].pha == reslist_All[j].ptype) and (phase2delete[i].evID != reslist_All[j].evID) and (abs(phase2delete[i].edep-hypos1[reslist_All[j].evidx].dep)<=RCC_Stat_rad_Dep):
                 #First test passed, now calculate difference in epicenter between two events:
                 epidiff = gps2dist_azimuth(phase2delete[i].elat,phase2delete[i].elon,hypos1[reslist_All[j].evidx].lat,hypos1[reslist_All[j].evidx].lon)
              
                 #Check if event is within radius:
                 if(epidiff[0]/1000.0 <= RCC_Stat_rad_Epi):
                    #Calculate distance,Azi,Baz between same station and similar event:
                    evepidiff = gps2dist_azimuth(hypos1[reslist_All[j].evidx].lat,hypos1[reslist_All[j].evidx].lon,phase2delete[i].slat,phase2delete[i].slon)

                    #store epicentral distance:
                    reslist_All[j].auxf1 = epidiff[0]/1000.0
                    reslist_All[j].auxf2 = abs(phase2delete[i].edep-hypos1[reslist_All[j].evidx].dep)
                    reslist_All[j].auxf3 = abs(phase2delete[i].OTtimestamp-hypos1[reslist_All[j].evidx].timestamp)/(60.0*60.0*24.0)
                    reslist_All[j].auxf4 = evepidiff[1] #AZI
                    reslist_All[j].auxf5 = evepidiff[2] #BAZ
                    residuals_SourceRegion.append(reslist_All[j])

              #2) Now test if residual matches SAME EVENT, phase and is from different station:
              if(phase2delete[i].stat != reslist_All[j].stat) and (phase2delete[i].pha == reslist_All[j].ptype) and (phase2delete[i].evID == reslist_All[j].evID):  

                 #Default flag:
                 insector = False

                 #Calculate azimuth between reslist_All-epicenter and reslist_All-station
                 epidiffepi = gps2dist_azimuth(phase2delete[i].elat,phase2delete[i].elon,stations[reslist_All[j].stidx].lat,stations[reslist_All[j].stidx].lon)

                 azi2test = epidiffepi[1]
                 dis2test = epidiffepi[0]/1000.0
                 baz2test = epidiffepi[2]

                 #Setup azimuth range defining the sector:
                 azi_min = phase2delete[i].azi - RCC_Even_bin_Azi
                 azi_max = phase2delete[i].azi + RCC_Even_bin_Azi

                 #Setup distance range defining the sector:
                 dis_min = phase2delete[i].epi - RCC_Even_rad_Sta
                 dis_max = phase2delete[i].epi + RCC_Even_rad_Sta

                 #Now test if within 0-359.9:
                 if(azi_min < 0):
                    azi_min = 360.0 + azi_min
                 if(azi_max >= 360.0):
                    azi_max = azi_max - 360.0

                 #Now check if zero-crossing occurs:
                 if(azi_min>azi_max):
                    #Zero-crossing case, e.g. 350-10 deg
                    if((azi2test <= azi_max) or (azi2test >= azi_min)) and (dis2test>=dis_min) and (dis2test<=dis_max):
                       insector = True
                       #print "Azi-Residual ZERO CROSSING in sector:",phase2delete[i].azi,azi_min,azi_max,dis_min,dis_max,azi2test,dis2test
                 else:
                    #Normal case, e.g. 10-30 deg:
                    if(azi2test>=azi_min) and (azi2test<=azi_max) and (dis2test>=dis_min) and (dis2test<=dis_max):
                       insector = True
                       #print "Azi-Residual               in sector:",phase2delete[i].azi,azi_min,azi_max,dis_min,dis_max,azi2test,dis2test

                 #calculate difference in station coordinates (previous test based on radius around station... replaced by segment-area, therefore not used anymore):
                 #epidiff = gps2dist_azimuth(phase2delete[i].slat,phase2delete[i].slon,stations[reslist_All[j].stidx].lat,stations[reslist_All[j].stidx].lon)

                 #Check if station is within sector:
                 #Old version radius based:
                 #if(epidiff[0]/1000.0 <= RCC_Even_rad_Sta):
                 #New version: Consider sector defined by azithuth range and distance:
                 if(insector):
                    #Calculate distance,Azi,Baz between similar station and same event (done above):
                    #evepidiff = gps2dist_azimuth(phase2delete[i].elat,phase2delete[i].elon,stations[reslist_All[j].stidx].lat,stations[reslist_All[j].stidx].lon)                         

                    #store epicentral distance/azi/baz:
                    reslist_All[j].auxf1 = dis2test  #DIS
                    reslist_All[j].auxf4 = azi2test  #AZI
                    reslist_All[j].auxf5 = baz2test  #BAZ

                    #reslist_All[j].auxf1 = epidiff[0]/1000.0
                    #reslist_All[j].auxf4 = evepidiff[1] #AZI
                    #reslist_All[j].auxf5 = evepidiff[2] #BAZ

                    #Update counter:
                    RCC_Check_CNT_RRR_Match += 1

                    #Append residual
                    residuals_ReceiverRegion.append(reslist_All[j])
                 
          #Sort results according to distance from reference station/event:
          residuals_SourceRegion.sort(key=lambda x: (x.auxf1), reverse=False)
          residuals_ReceiverRegion.sort(key=lambda x: (x.auxf1), reverse=False)

          #Calculate mean and standard deviation of residuals_SourceRegion to calculate threshold:
          mean_RSR = 0.0
          stde_RSR = 0.0
          nobs_RSR = len(residuals_SourceRegion)
          nobs_RSRTGap = 0 #Number of observations outside of emporal gap
          #Case #1: No observations: -> Nothing to do...
          #Case #2: Only one observation: STD does't exists, use static tolerance instead
          if(nobs_RSR > 0):
             for k in range(nobs_RSR):
                 if(residuals_SourceRegion[k].auxf3/30.0 >= RCC_TemporalGap):
                    mean_RSR += residuals_SourceRegion[k].res
                    nobs_RSRTGap += 1
             if(nobs_RSRTGap > 0):
                mean_RSR = mean_RSR/nobs_RSRTGap
             for k in range(nobs_RSR):
                 if(residuals_SourceRegion[k].auxf3/30.0 >= RCC_TemporalGap):
                    stde_RSR += ((residuals_SourceRegion[k].res - mean_RSR)**2)
             if(nobs_RSRTGap > 1):
                stde_RSR = math.sqrt(stde_RSR/(nobs_RSR-1))
             else:
                #If only one residual exists, use static value for tolerance:
                stde_RSR = RCC_ResStaticTolerance

          #Calculate mean and standard deviation of residuals_ReceiverRegion to calculate threshold:
          mean_RRR = 0.0
          stde_RRR = 0.0
          nobs_RRR = len(residuals_ReceiverRegion)
          #Case #1: No observations: -> Nothing to do...
          #Case #2: Only one observation: STD does't exists, use static tolerance instead
          if(nobs_RRR > 0):
             for k in range(nobs_RRR):
                 mean_RRR += residuals_ReceiverRegion[k].res
             mean_RRR = mean_RRR/nobs_RRR
             for k in range(nobs_RRR):
                 stde_RRR += ((residuals_ReceiverRegion[k].res - mean_RRR)**2)
             if(nobs_RRR > 1):       
                stde_RRR = math.sqrt(stde_RRR/(nobs_RRR-1))
             else:
                #If only one residual exists, use static value for tolerance:
                stde_RRR = RCC_ResStaticTolerance

          #Count number of similar residuals:
          RCC_DeletFlagI = 1      #By default remove potential outlier...
          RCC_DeletFlagS = 'DELETE'
          cnt_SourceRegion = 0
          cnt_ReceiverRegion = 0
          #Define tolerance level for counting similar residuals
          if(RCC_ResDynamicTolerance):
             RCC_ResStaticTolerance_RSR_Cnt = stde_RSR
             RCC_ResStaticTolerance_RRR_Cnt = stde_RRR
          else:
             RCC_ResStaticTolerance_RSR_Cnt = RCC_ResStaticTolerance
             RCC_ResStaticTolerance_RRR_Cnt = RCC_ResStaticTolerance

          for k in range (len(residuals_SourceRegion)):
              #if(residuals_SourceRegion[k].res >= RCC_RES_low) and (residuals_SourceRegion[k].res <= RCC_RES_upp) and (residuals_SourceRegion[k].auxf3/30.0 >= RCC_TemporalGap):
              if(abs(phase2delete[i].res-residuals_SourceRegion[k].res) <= RCC_ResStaticTolerance_RSR_Cnt) and (residuals_SourceRegion[k].auxf3/30.0 >= RCC_TemporalGap):
                 cnt_SourceRegion += 1
         
          for k in range (len(residuals_ReceiverRegion)):
              #if(residuals_ReceiverRegion[k].res >= RCC_RES_low) and (residuals_ReceiverRegion[k].res <= RCC_RES_upp):
              if(abs(phase2delete[i].res-residuals_ReceiverRegion[k].res) <= RCC_ResStaticTolerance_RRR_Cnt):
                 cnt_ReceiverRegion += 1       

          #Final decision:
          #Dynamic threshold:
          if(RCC_ResDynamicTolerance):
             #print "--> Apply dynamic RCC"
             if(((nobs_RSRTGap>=RCC_Min_RSR) and (phase2delete[i].res>=mean_RSR-stde_RSR) and (phase2delete[i].res<=mean_RSR+stde_RSR)) or (cnt_SourceRegion>=RCC_MinCount) or ((nobs_RRR>=RCC_Min_RRR) and (phase2delete[i].res>=mean_RRR-stde_RRR) and (phase2delete[i].res<=mean_RRR+stde_RRR)) or (cnt_ReceiverRegion>=RCC_MinCount)):
                RCC_DeletFlagI = 0
                RCC_DeletFlagS = 'KEEP'
                RCC_Check_CNT += 1                                               #Update counter for residuals which passed RCC
                #Update score / flag for use...
                phase2delete[i].oscore = 0.0                                     #Lower oscore from 1 to 0.0, should be used if threshold is 0.9
                phases1[phase2delete[i].phaidx].oscore = phase2delete[i].oscore  #use same oscore also for original phase
          else:
             #print "--> Apply static RCC"
             if(cnt_SourceRegion>=RCC_MinCount) or (cnt_ReceiverRegion>=RCC_MinCount):
                RCC_DeletFlagI = 0
                RCC_DeletFlagS = 'KEEP'
                RCC_Check_CNT += 1                                               #Update counter for residuals which passed RCC
                #Update score / flag for use...
                phase2delete[i].oscore = 0.0                                     #Lower oscore from 1 to 0.0, should be used if threshold is 0.9
                phases1[phase2delete[i].phaidx].oscore = phase2delete[i].oscore  #use same oscore also for original phase

          #document residual analysis: residual: stat,ptype,qual,err,weig,res,dis,evID,AvResInQuad,NoResInQuad,Quad,evidx,stidx
          #Write basic information related to residual to check:
          #Total number of observations for station/phase:
          fplogRCC.write("#                  Stat     Phase   Residual Q   Error   Dist AZI BAZ       elat       elon   edep eMag Event-ID                   slat       slon StatOrg  NobsTot DeleteFlg noOTGap-SS av-SS sd-SS cntSS   | no-ER av-ER sd-ER cntER\n")
          #               RCC-FOR-RESIDUAL : LKB1     P          0.685 4   0.400   13.4  68 249   46.34280    7.46510   13.1  1.6 KP199601042354         46.38700    7.62710 LKBD
          #               RCC-FOR-RESIDUAL : PLOS     P         -0.421 0   0.025   32.7 351 171   46.75900    9.44700    7.0  4.0 KP200801211640         47.04920    9.38080 PLONS        360 DELETE            8  0.032 0.059 0000   0 |  0.000 0.100 0000
          fplogRCC.write("RCC-FOR-RESIDUAL : %-8s %-8s %7.3f %1d %7.3f %6.1f %3.0f %3.0f %10.5f %10.5f %6.1f %4.1f %-20s %10.5f %10.5f %-8s %7d %-8s       %4d %6.3f %5.3f %4d | %4d %6.3f %5.3f %4d\n" % (phase2delete[i].stat,phase2delete[i].pha,phase2delete[i].res,phase2delete[i].qua,phase2delete[i].err,phase2delete[i].epi,phase2delete[i].azi,phase2delete[i].baz,phase2delete[i].elat,phase2delete[i].elon,phase2delete[i].edep,phase2delete[i].emag,phase2delete[i].evID,phase2delete[i].slat,phase2delete[i].slon,phase2delete[i].statorg,phase2delete[i].SCnob,RCC_DeletFlagS,nobs_RSRTGap,mean_RSR,stde_RSR,cnt_SourceRegion,nobs_RRR,mean_RRR,stde_RRR,cnt_ReceiverRegion))

          #Document all residuals associated to similar sources for same station:
          if(len(residuals_SourceRegion)>0):
             fplogRCC.write("#                  Stat     Phase   Residual Q   Error   Dist AZI BAZ       elat       elon   edep eMag Event-ID             DiffEpi DiffDep DiffOT (Months)\n")
             #               Station2Source   : LKB1     P         -0.487 4   0.400   49.0 136 317   46.70880    7.18610    7.1  2.9 KP199606272239           0.6    1.9    260
          for k in range(len(residuals_SourceRegion)):
              fplogRCC.write("Station2Source   : %-8s %-8s %7.3f %1d %7.3f %6.1f %3.0f %3.0f %10.5f %10.5f %6.1f %4.1f %-20s %7.1f %6.1f %6d\n" % (residuals_SourceRegion[k].stat,residuals_SourceRegion[k].ptype,residuals_SourceRegion[k].res,residuals_SourceRegion[k].qual,residuals_SourceRegion[k].err,residuals_SourceRegion[k].dis,residuals_SourceRegion[k].auxf4,residuals_SourceRegion[k].auxf5,hypos1[residuals_SourceRegion[k].evidx].lat,hypos1[residuals_SourceRegion[k].evidx].lon,hypos1[residuals_SourceRegion[k].evidx].dep,hypos1[residuals_SourceRegion[k].evidx].mag,residuals_SourceRegion[k].evID,residuals_SourceRegion[k].auxf1,residuals_SourceRegion[k].auxf2,residuals_SourceRegion[k].auxf3/30.0))
          fplogRCC.write("#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")

          #Document all residuals associated to similar receiver for same event:
          for k in range(len(residuals_ReceiverRegion)):
              fplogRCC.write("Event2Receiver   : %-8s %-8s %7.3f %1d %7.3f %6.1f %3.0f %3.0f %10.5f %10.5f %6.1f %4.1f %-20s %7.1f\n" % (residuals_ReceiverRegion[k].stat,residuals_ReceiverRegion[k].ptype,residuals_ReceiverRegion[k].res,residuals_ReceiverRegion[k].qual,residuals_ReceiverRegion[k].err,residuals_ReceiverRegion[k].dis,residuals_ReceiverRegion[k].auxf4,residuals_ReceiverRegion[k].auxf5,hypos1[residuals_ReceiverRegion[k].evidx].lat,hypos1[residuals_ReceiverRegion[k].evidx].lon,hypos1[residuals_ReceiverRegion[k].evidx].dep,hypos1[residuals_ReceiverRegion[k].evidx].mag,residuals_ReceiverRegion[k].evID,residuals_ReceiverRegion[k].auxf1))

          #Blank line
          fplogRCC.write("\n")

######################################################################################################################################
#Loop #5: Loop over all events: 1) Dynamic distance filter for e.g. Pg tomography:
#                                  --> Extract all used phases for an event
#                                  --> Sort by distance
#                                  --> Loop over all phases, change nothing for distances <= maxepidist1, use phases between maxepidist1 and maxepidist2 only
#                                      if they improve nobs or gap.

if(maxepidist2 <= 600.0) and (ckQuakes == 0) and (rmOutliers == 0) and (rphypos == 0) and (rmQuakes == 0):
   print ""
   print "--> Dynamic distance selection for distances between",maxepidist1,"and",maxepidist2,"and for S:",maxepidist1S,"and",maxepidist2S

   #Loop over all selected events:
   for i in range(len(filteredHypo)):

       #Extract phases associated with current event:
       tmpPhases1 = []

       #Loop over all phases associated with this event:
       for j in range(filteredHypo[i].stphaidx, filteredHypo[i].enphaidx+1):

           #Final check of event IDs and if phase is still used:
           if(phases1[j].evID == filteredHypo[i].evID1) and (phases1[j].qual<wscheme1[len(wscheme1)-1].qual) and (phases1[j].weig>0.00001) and (phases1[j].oscore <= oscoreMAX):

              #Add phase to temp. phase-list for event:
              tmpPhases1.append(phases1[j])

       #Now sort phases in tmp-array according to:
       # 1) Distance 
       # 2) Phase-Name P, then S. / alternatively tt
       tmpPhases1.sort(key=lambda x: (x.edist, x.ptype), reverse=False)

       tmpNob = 0
       tmpGAP = 360
       tmpAzi = []

       #stat,ptype,qual,err,weig,tt,timestamp,edist,azi,baz,evID,slat,slon,stidx,oscore,phaidx
       #Now loop starting from smallest distance and test if nobs and gap are sufficient...
       for j in range(len(tmpPhases1)):

           #Comment TD: Problem in previous version: phase-j was counted for nobs & gap, but then possibly removed...
           #            the improved algorithm checks now nobs and gap from previous loop before consider the phase used/unused
           #            In addition S-phases can be added only for certain distance ranges (until quality is sufficient)         
           #2019/03/27

           #Check if phase still needed:
           #Distance is above maxepidist1 (e.g. 80 km)
           #Total number of phases is already > 10
           #Gap is already < 170
           #-> phase not needed anymore
           if(tmpPhases1[j].edist > maxepidist1)  and (tmpNob > dyndist_minnob) and (tmpGAP < dyndist_maxgap):
              #Phase is beyond maxepidist1, Number of phases and gap are already sufficient, we don't need that phase anymore
              phases1[tmpPhases1[j].phaidx].oscore = 9.9    #Mark original phase as not used... 
              phases1[tmpPhases1[j].phaidx].weig = 0.0
           else:
              #Phase is potentially still needed to get high-quality solutions
              #Check if phase is S. Only use S if 1) distance <= maxepidist1S
              #                                or 2) Minimum quality not yet reached
              if(tmpPhases1[j].ptype.startswith("S")) and (tmpPhases1[j].edist > maxepidist1S) and (tmpNob > dyndist_minnobS) and (tmpGAP < dyndist_maxgapS):
                 #S-phase beyond 25 km AND minimum quality for location reached, do not use that phase anymore
                 #print "--> Phase removed: S-phase >",maxepidist1S,"and minimum quality reached...",tmpPhases1[j].edist,tmpNob,tmpGAP
                 phases1[tmpPhases1[j].phaidx].oscore = 9.9    #Mark original phase as not used... 
                 phases1[tmpPhases1[j].phaidx].weig = 0.0
              else:
                 #Phase is still needed -> Update nobs and gap:
                 #Update number of phases nobs:
                 tmpNob += 1

                 #Update GAP:
                 tmpAzi.append(tmpPhases1[j].azi)
                 if(len(tmpAzi) > 1):
                    tmpGAP,newsgap = getGap(tmpAzi)
           
######################################################################################################################################
#Loop #6: Loop over all events: 1) Decide if events passe minimum criteria (number of observations, gap, RMS) after potential outliers have been removved
#                                   --> only events with "mergedHypo[i].useflg == 1" will be considered for final 

print ""
print "--> Final check of event quality after removal of outlier phases..."

#Loop over all selected events:
for i in range(len(filteredHypo)):

    #Extract phases associated with current event:
    tmpPhases1 = []

    #Loop over all phases associated with this event:
    for j in range(filteredHypo[i].stphaidx, filteredHypo[i].enphaidx+1):

        #Final check of event IDs and if phase is still used:
        if(phases1[j].evID == filteredHypo[i].evID1) and (phases1[j].qual<wscheme1[len(wscheme1)-1].qual) and (phases1[j].weig>0.00001) and (phases1[j].oscore <= oscoreMAX):
        
           #Add phase to temp. phase-list for event:
           tmpPhases1.append(phases1[j])
           
    #New version (updateEventPhaseData): Use info which was already added by getDistAzi4phase and only update the event info + weights...
    newgap,newmdist,newweightsum,newNobs,tmpPhases2 = updateEventPhaseData(tmpPhases1,filteredHypo[i].evID1,wscheme1)

    #Sort phases according to:
    # 1) Distance 
    # 2) Phase-Name P, then S. / alternatively tt
    #tmpPhases2.sort(key=lambda x: (x.edist, x.ptype), reverse=False)
    
    #Document significant changes in gap:
    if(abs(filteredHypo[i].gap-newgap) > 5.0):
       fplog.write("# GAPDIFF %-20s %6.1f %5.1f %5.1f\n" % (filteredHypo[i].evID1,filteredHypo[i].gap-newgap,filteredHypo[i].gap,newgap))

    #Update hypocenter information
    filteredHypo[i].gap = newgap
    filteredHypo[i].mdist = newmdist
    filteredHypo[i].nobs = newNobs        #Includes only observation with weight > 0    
    filteredHypo[i].wgsum = newweightsum  #Sum of weights
    filteredHypo[i].matchidx = i          #Save position in original list

    #Check if updated GAP is acceptable:
    if(filteredHypo[i].gap > co_maxGap) and (filteredHypo[i].useflg == 1):
        print "--> WARNING: Updated GAP >",co_maxGap," --> After VELEST RELOCATION                            --> skip event for final selection:",filteredHypo[i].evID1,filteredHypo[i].lat,filteredHypo[i].lon,filteredHypo[i].gap,filteredHypo[i].mag
        fplog.write("# GAPTOOBIG   : %-s %5.1f %5.1f\n" % (filteredHypo[i].evID1,filteredHypo[i].gap,filteredHypo[i].mag))   
        fplogE2D.write("DEL-EVENT: GAP-TOO-BIG %6.3f %5.1f %5d %s %9.4f %8.4f %8.3f %5.1f %s\n" % (filteredHypo[i].rms,filteredHypo[i].gap,filteredHypo[i].nobs,UTCDateTime(filteredHypo[i].timestamp),filteredHypo[i].lon,filteredHypo[i].lat,filteredHypo[i].dep,filteredHypo[i].mag,filteredHypo[i].evID1))
        filteredHypo[i].useflg = 0
        no_removed_events += 1
    
    #Check if updated Nobs is acceptable:
    if(filteredHypo[i].nobs < co_minNob) and (filteredHypo[i].useflg == 1):
        print "--> WARNING: Updated Nobs <",co_minNob," --> After residual filtering...                        --> skip event for final selection:",filteredHypo[i].evID1,filteredHypo[i].lat,filteredHypo[i].lon,filteredHypo[i].nobs,filteredHypo[i].mag
        fplog.write("# NOBSTOOSMALL: %-s %5d %5.1f\n" % (filteredHypo[i].evID1,filteredHypo[i].nobs,filteredHypo[i].mag))
        fplogE2D.write("DEL-EVENT: NOB-TOO-LOW %6.3f %5.1f %5d %s %9.4f %8.4f %8.3f %5.1f %s\n" % (filteredHypo[i].rms,filteredHypo[i].gap,filteredHypo[i].nobs,UTCDateTime(filteredHypo[i].timestamp),filteredHypo[i].lon,filteredHypo[i].lat,filteredHypo[i].dep,filteredHypo[i].mag,filteredHypo[i].evID1))
        filteredHypo[i].useflg = 0
        no_removed_events += 1

    #Check if RMS after simult, inversion is still acceptable:
    if(filteredHypo[i].rms > co_maxRMS) and (filteredHypo[i].useflg == 1):
        print "--> WARNING:         RMS  >",co_maxRMS," --> After VELEST RELOCATION                            --> skip event for final selection:",filteredHypo[i].evID1,filteredHypo[i].lat,filteredHypo[i].lon,filteredHypo[i].rms,filteredHypo[i].mag
        fplog.write("# RMSTOOBIG   : %-s %5.3f %5.1f\n" % (filteredHypo[i].evID1,filteredHypo[i].rms,filteredHypo[i].mag))
        fplogE2D.write("DEL-EVENT: RMS-TOO-BIG %6.3f %5.1f %5d %s %9.4f %8.4f %8.3f %5.1f %s\n" % (filteredHypo[i].rms,filteredHypo[i].gap,filteredHypo[i].nobs,UTCDateTime(filteredHypo[i].timestamp),filteredHypo[i].lon,filteredHypo[i].lat,filteredHypo[i].dep,filteredHypo[i].mag,filteredHypo[i].evID1))
        filteredHypo[i].useflg = 0
        no_removed_events += 1

######################################################################################################################################
#Write output: loop over all events selected for final declustered catalog:

#Open output phase file:
if(OutputForm != 'MANUPICK') and (OutputForm != 'MANUPICK_alias2org'):
   fpcnv = open(OutputFile,'w')
else:
   #Check if output directory exists, if not create:
   if((os.path.isdir(OutputFile)) is False):
      print "mkdir",OutputFile
      os.makedirs(OutputFile)

#Write header of log files
fplogEv.write("#lon      lat      dep  mag  yy   mm dd hh mi ss    rms   gap mdist  nobs  wsum  used evscore gridID TotEventsPerGrid EventID\n")
fplogPh.write("#Stat  Phase QA  err     wgt    traveltime   Picktime                  dist     azi    evID          elon     elat      slon      slat\n")
fplogSt.write("#Stat  lon lat elev NPobs Pwsum NSobs Swsum\n")

cntfilter_phase = 0
cntfilter_event = 0

#Loop over all selected events + update station counter!
for i in range(len(filteredHypo)):

    #Check is event was selected for declustered catalog
    if(filteredHypo[i].useflg > 0):

       #Extract phases associated with current event:
       tmpPhases1 = []

       cntfilter_event += 1

       #Collect all phases for event:
       #Reminder: stat,ptype,qual,tt,timestamp,evID
       for j in range(filteredHypo[i].stphaidx, filteredHypo[i].enphaidx+1):
           #Adjust travel time to correct origin time in merged catalog
           phases1[j].tt = phases1[j].timestamp - filteredHypo[i].timestamp

           #Final check of event IDs and if weight is > 0 / phase is not marked at outlier:
           if(phases1[j].evID == filteredHypo[i].evID1) and (phases1[j].qual <= co_minQual) and (phases1[j].oscore <= oscoreMAX):
              tmpPhases1.append(phases1[j])

       #Write final events to log file with all quality parameters necessary for declustering:
       fplogEv.write("%9.4f %8.4f %5.1f %4.1f %4d %02d %02d %02d %02d %04.1f %6.3f %3.0f %5.1f %3d %5.1f %3d %6.3f %9d %5d %-20s\n" % (filteredHypo[i].lon,filteredHypo[i].lat,filteredHypo[i].dep,filteredHypo[i].mag,filteredHypo[i].yy,filteredHypo[i].mm,filteredHypo[i].dd,filteredHypo[i].hh,filteredHypo[i].mi,filteredHypo[i].ss,filteredHypo[i].rms,filteredHypo[i].gap,filteredHypo[i].mdist,filteredHypo[i].nobs,filteredHypo[i].wgsum,filteredHypo[i].useflg,filteredHypo[i].evscore,filteredHypo[i].gid,filteredHypo[i].catID,filteredHypo[i].evID1))

       #Write phases to log file for statistics and visualization: stat,ptype,qual,err,weig,tt,timestamp,edist,azi,evID --> Problem: the eventID in the phase data should be adjusted for merged events
       for m in range(len(tmpPhases1)):
           cntfilter_phase += 1
           fplogPh.write("%-6s %-5s %2d %8.5f %7.5f %8.3f %-s %7.2f %5.1f %-20s %9.4f %8.4f %9.4f %8.4f\n" % (tmpPhases1[m].stat,tmpPhases1[m].ptype,tmpPhases1[m].qual,tmpPhases1[m].err,tmpPhases1[m].weig,tmpPhases1[m].tt,UTCDateTime(tmpPhases1[m].timestamp),tmpPhases1[m].edist,tmpPhases1[m].azi,tmpPhases1[m].evID,filteredHypo[i].lon,filteredHypo[i].lat,tmpPhases1[m].slon,tmpPhases1[m].slat))

           #Update station  and phase counter:
           if(tmpPhases1[m].ptype[0:1] == 'P'):
              stations[tmpPhases1[m].stidx].nobsP += 1
              stations[tmpPhases1[m].stidx].nobsPw = stations[tmpPhases1[m].stidx].nobsPw + tmpPhases1[m].weig
              wscheme1[tmpPhases1[m].qual].nobsP += 1             

           if(tmpPhases1[m].ptype[0:1] == 'S'):
              stations[tmpPhases1[m].stidx].nobsS += 1
              stations[tmpPhases1[m].stidx].nobsSw = stations[tmpPhases1[m].stidx].nobsSw + tmpPhases1[m].weig
              wscheme1[tmpPhases1[m].qual].nobsS += 1

       #Write hypocenter to file:   
       if(OutputForm == 'cnv'):
          #Write hypocenter to file:
          writeCNVHeader(fpcnv,filteredHypo[i])

          #Write phases to file:
          writeCNVPhases(fpcnv,tmpPhases1)

       if(OutputForm == 'simulps'):
          #Write hypocenter to file:
          writeSIMHeader(fpcnv,filteredHypo[i])

          #Write phases to file:
          writeSIMPhases(fpcnv,tmpPhases1)

       if(OutputForm == 'MANUPICK') or (OutputForm == 'MANUPICK_alias2org'):
          #Setup filename:
          if(filteredHypo[i].evID1[0:2]=='KP'):
             diryy = filteredHypo[i].evID1[2:6]
             dirmm = filteredHypo[i].evID1[6:8]
          else:
             diryy = filteredHypo[i].evID1[0:4]
             dirmm = filteredHypo[i].evID1[4:6]
          pathdir = OutputFile + '/' + diryy + '/' + dirmm

          #Check if output directory exists, if not create:
          if((os.path.isdir(pathdir)) is False):
             print "mkdir",pathdir
             os.makedirs(pathdir)

          MPfile = pathdir + '/' + filteredHypo[i].evID1 + '.MANUPICK'
          fpMP  = open(MPfile,'w')
 
          #Write MANUPICK header:
          #yy,mm,dd,hh,mi,ss,timestamp,lon,lat,dep,lonerr,laterr,deperr,mag,mtype,mnobs,merr,mmeth,rms,gap,mdist,nobs,wgsum,etype,lqual,chx,chy,agency,evID1,evID2,evPID,orID1,methodID,earthmodelID,author,region,fstatus,difflon,difflat,diffdep,diffepi,diffoti,diffazi,latalt,lonalt,depalt,catID,stphaidx,enphaidx,matchidx,evscore,gid,useflg)
          disrange = "Loca"
          OTime = UTCDateTime(filteredHypo[i].timestamp)
          rdate = OTime - 60.0
          #"2009-12-31T12:23:34.5"
          #UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (yy,mm,dd,hh,mi,ss))
          #Final reference data is OT - one minute and second = 0.0
          rdate = UTCDateTime("%04d-%02d-%02d %02d:%02d:%09.6f" % (rdate.year,rdate.month,rdate.day,rdate.hour,rdate.minute,0.0))
          fpMP.write(" 1 1  200.00  300.00   1.710 2                        INST\n")
          fpMP.write("   0.000     0.000  10.00 0000 00 00 00 00  0.00 0    TRIAL\n")
          fpMP.write("%-53s AUTOR\n" % ('CNV2KP'))
          fpMP.write("%04d/%02d/%02d %02d:%02d                                      %-s\n" % (rdate.year,rdate.month,rdate.day,rdate.hour,rdate.minute,disrange))

          #Loop over phases:
          #stat,ptype,qual,err,weig,tt,timestamp,edist,azi,baz,evID,slat,slon,stidx,oscore,phaidx
          for m in range(len(tmpPhases1)):
              tt=tmpPhases1[m].timestamp - rdate.timestamp
              if(tt >= 1000.0):
                 print "WARNING: traveltime exceeds 999.9 seconds, skip phase in MANUPICK-export: ",tmpPhases1[m].stat,tmpPhases1[m].ptype,tt,filteredHypo[i].evID1
                 continue
              #Write standard MANUPICK (use station-name provided as in input phase file):
              if(OutputForm == 'MANUPICK'):
                 fpMP.write("%-7s %-7s %-1s%-1s %7.3f %1d %4.2f %8d %-9s %1d Pick %-3s %7.3f %7.3f %7.3f %1d\n" % (tmpPhases1[m].stat,tmpPhases1[m].ptype,' ',' ',tt,1,0,0,'None',1,'???',0,tmpPhases1[m].err*-1.0,tmpPhases1[m].err,tmpPhases1[m].qual))
              #Write MANUPICK_alias2org -> Replace name in input phase file with ORIGINAL station name provided in alias list (useful if you want to use file for NLL):
              else:
                 fpMP.write("%-7s %-7s %-1s%-1s %7.3f %1d %4.2f %8d %-9s %1d Pick %-3s %7.3f %7.3f %7.3f %1d\n" % (stations[tmpPhases1[m].stidx].statorg,tmpPhases1[m].ptype,' ',' ',tt,1,0,0,'None',1,'???',0,tmpPhases1[m].err*-1.0,tmpPhases1[m].err,tmpPhases1[m].qual))

          #Finish MANUPICK file:
          fpMP.write("                                                      SKIP\n")

          fpMP.close()

####################################################################################################################
#Write station statistics: stat,netw,lat,lon,ele,nobsP,nobsS,nobsPw,nobsSw
for l in range(len(stations)):
    fplogSt.write("%-6s %9.4f %8.4f %6.3f %5d %7.2f %5d %7.2f\n" % (stations[l].stat,stations[l].lon,stations[l].lat,stations[l].ele,stations[l].nobsP,stations[l].nobsPw,stations[l].nobsS,stations[l].nobsSw))

####################################################################################################################
#Write pick statistics: qual,err,weig,nobsP,nobsS
averrP = 0
npicksP = 0
averrS = 0
npicksS = 0
print ""
print "--> Class error weight NobsP NobsS"
for i in range(len(wscheme1)):
    print "    %1d %6.3f %7.4f %9d %9d" % (wscheme1[i].qual,wscheme1[i].err,wscheme1[i].weig,wscheme1[i].nobsP,wscheme1[i].nobsS)
    fplog.write("#PICK-STAT: %1d %6.3f %7.4f %9d %9d\n" % (wscheme1[i].qual,wscheme1[i].err,wscheme1[i].weig,wscheme1[i].nobsP,wscheme1[i].nobsS))
    if(wscheme1[i].weig > 0):
       averrP = averrP + (wscheme1[i].nobsP*wscheme1[i].err)
       npicksP = npicksP + wscheme1[i].nobsP
       averrS = averrS + (wscheme1[i].nobsS*wscheme1[i].err)
       npicksS = npicksS + wscheme1[i].nobsS

if(npicksP > 0):
   averrP = averrP/npicksP
   print "    Average P-Pick error: %6.3f" % (averrP)
   fplog.write("#PICK-STAT-AVERAGE-P-Err: %6.3f" % (averrP))

if(npicksS > 0):
   averrS = averrS/npicksS
   print "    Average S-Pick error: %6.3f" % (averrS) 
   fplog.write("#PICK-STAT-AVERAGE-S-Err: %6.3f" % (averrS))

####################################################################################################################
#Write Outlier information:

#Sort list of outliers:
if(ckQuakes == 0):
   phase2delete.sort(key=lambda x: (x.stat, x.epi), reverse=False)

for i in range(len(phase2delete)):
    fplogP2D.write("DEL-PHASE: %-6s %1s %1d %7.3f %s %7.3f | %5d %7.3f %7.3f %7.3f %7.3f %5d %2d | %7.1f %5.1f %5.1f %9.4f %8.4f %7.1f %7.1f | %s %9.4f %8.4f %8.3f %5.1f %-6s %7.3f %5d %s\n" % (phase2delete[i].stat,phase2delete[i].pha,phase2delete[i].qua,phase2delete[i].err,UTCDateTime(phase2delete[i].picktimestamp),phase2delete[i].res,phase2delete[i].SCnob,phase2delete[i].SCval,phase2delete[i].SCasr,phase2delete[i].SCssr,phase2delete[i].SCasrQ,phase2delete[i].SCnsrQ,phase2delete[i].Quad,phase2delete[i].epi,phase2delete[i].azi,phase2delete[i].baz,phase2delete[i].slon,phase2delete[i].slat,phase2delete[i].asepi,phase2delete[i].ssepi,UTCDateTime(phase2delete[i].OTtimestamp),phase2delete[i].elon,phase2delete[i].elat,phase2delete[i].edep,phase2delete[i].emag,phase2delete[i].statorg,phase2delete[i].oscore,phase2delete[i].WFrevf,phase2delete[i].evID)) 

####################################################################################################################

#Close output files:       
fplog.close()
if(OutputForm != 'MANUPICK') and (OutputForm != 'MANUPICK_alias2org'):
   fpcnv.close()
fplogEv.close()
fplogPh.close()
fplogSt.close()
fplogSC.close()
fplogP2D.close()
fplogE2D.close()
if(RCC_Check_Flg==1):
   fplogRCC.close()

print ""
print "Phase file                                                         :",PhaseFile1
print "Station file                                                       :",StationFile
print "Residual file                                                      :",ResidualFile
print "Original number of phases loaded                                   :",org_no_phases
print "         number of phases to remove based on residual-analysis     :",len(reslist)
print "File with filtered    phases of accepted events                    :",OutputFile
print "Number of events loaded from input phase-file                      :",len(hypos1)
print "Number of events to delete in quakes2delete file                   :",len(quakes2delete)
print "Number of events removed because marked in quakes2delete file      :",len(hypos1)-len(filteredHypo) 
print "Number of events prior filtering                                   :",len(filteredHypo)  
print "Number of events after filtering                                   :",cntfilter_event
print "Number of phases after filtering                                   :",cntfilter_phase
print "Statistics merging       written to file                           :",LogFileFilter
print "Log file   of events passing quality criteria writen to            :",LogFileEvents
print "Log file   of phases passing quality criteria writen to            :",LogFilePhases
print "Log file   with stations used in final cnv file + pick statistics  :",LogFileStations
print "Log file   with stations of sign. diff. in av. residual per azimuth:",LogFileStatCorr
print "Log file   with phases which have been removed                     :",LogFilePhas2Del
print "Log file   with events which have been removed                     :",LogFileEven2Del
if(RCC_Check_Flg==1):
   print "Log file  with statistics and results of Residual Consistency Check:",LogFileRCC
   print "Number of suspicious residuals which passed Residual Consist. Check:",RCC_Check_CNT
   print "Number of station in Event2Receiver sectors Residual Consist. Check:",RCC_Check_CNT_RRR_Match

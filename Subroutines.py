# -*- coding: utf-8 -*-
"""
Created on Tue May  3 08:07:12 2016

@author: justinpringle
Subroutines for SBEACH
"""
import numpy as np
import datetime as dt

def allocate(sb):
    '''
    allocate arrays - set ot zero
    '''
    ndx = sb.arrsize.ndx
    sb.memorymain.dx = np.zeros(ndx)
    sb.memorymain.prix = np.zeros(ndx)
    sb.memorymain.prid = np.zeros(ndx)
    sb.memorymain.hbx = np.zeros(ndx)
    sb.memorymain.hbd = np.zeros(ndx)
    sb.memorymain.thmx = np.zeros(ndx)
    sb.memorymain.temx = np.zeros(ndx)
    sb.memorymain.tdpmx = np.zeros(ndx)
    sb.memorymain.tdmx = np.zeros(ndx)
    sb.memorymain.tdmn = np.zeros(ndx)
    sb.memorymain.ihb = np.zeros(ndx)
    sb.memorymain.di = np.zeros(ndx)
    sb.memorymain.dp = np.zeros(ndx)
    sb.memorymain.d = np.zeros(ndx)
    sb.memorymain.x = np.zeros(ndx)
    sb.memorymain.diss = np.zeros(ndx)
    sb.memorymain.xsv = np.zeros(ndx)
    sb.memorymain.dsv = np.zeros(ndx)
    sb.memorymain.teta = np.zeros(ndx)
    sb.memorymain.ateta = np.zeros((ndx,4))
    sb.memorymain.dmx = np.zeros(ndx)
    sb.memorymain.dmn = np.zeros(ndx)
    sb.memorymain.hmx = np.zeros(ndx)
    sb.memorymain.emx = np.zeros(ndx)
    sb.memorymain.dpmx = np.zeros(ndx)
    sb.memorymain.volch = np.zeros(ndx)
    sb.memorymain.brok = np.zeros(ndx)
    sb.memorymain.dsmo = np.zeros(ndx)
    sb.memorymain.edhdx = np.zeros(ndx)
    sb.memorymain.dbot = np.zeros(ndx)
    sb.memorymain.da = np.zeros(ndx)
    
    sb.memorymain.q = np.zeros((ndx+1))
    sb.memorymain.qp = np.zeros((ndx+1))
    sb.memorymain.h =np.zeros((ndx+1))    
    sb.memorymain.e = np.zeros((ndx+1))
    sb.memorymain.zw = np.zeros((ndx+1))
    
def temporalInfo(H,D,T,WL,Wi,globVars):
    '''
    reads in temporal info used to extrapolate between indices.
    '''
    
    if globVars.iwave == 1:
        globVars.twavl = 0
        globVars.twavh = globVars.dtwav
        globVars.hinl = H[0]
        globVars.hinh = H[1]
        globVars.tl = T[0]
        globVars.th = T[1]
        
    if globVars.iang == 1:
        globVars.tangl = 0
        globVars.tangh = globVars.dtang
        globVars.zinl = D[0]
        globVars.zinh = D[1]
        
    if globVars.ielev == 1:
        globVars.telvl = 0
        globVars.telvh = globVars.dtelv
        globVars.elvl = WL[0]
        globVars.elvh = WL[1]
        
    if globVars.iwind == 1:
        globVars.twndl = 0
        globVars.twndh = globVars.dtwnd
        globVars.wl = Wi[0,0]
        globVars.wh = Wi[1,0]
        globVars.zwndl = Wi[0,1]
        globVars.zwndh = Wi[1,1]
        
def readWaves(time,H,T,D,WL,Wi,globVars,sbeach):
    '''
    reads temporal info at each time step and interpolates where neccessary.
    '''
    rad = 0.01745 #convert degrees to radians
    if globVars.iwave == 1:
        if time < globVars.twavh:
            tratio = (time-globVars.twavl)/globVars.dtwav
            globVars.hin = globVars.hinl+ (globVars.hinh-globVars.hinl)*tratio
            globVars.tin = globVars.tl+ (globVars.th-globVars.tl)*tratio
        else:
            globVars.hin = globVars.hinh
            globVars.tin = globVars.th
            globVars.twavl = globVars.twavh
            globVars.twavh = globVars.twavl+globVars.dtwav
            
            globVars.tl = globVars.th
            globVars.hinl = globVars.hinh
            
            globVars.hinh = H[int(globVars.twavh/globVars.dtwav)]
            globVars.th = T[int(globVars.twavh/globVars.dtwav)]
        globVars.hin*=0.706 #convert to rms wave height
        sbeach.memory.lo = 1.5613*globVars.tin**2    
    
    if globVars.iang == 1:
        if time < globVars.tangh:
            tratio = (time-globVars.tangl)/globVars.dtang
            globVars.zin = globVars.zinl+ (globVars.zinh-globVars.zinl)*tratio
        else:
            globVars.zin = globVars.zinh            
            globVars.tangl = globVars.tangh
            globVars.tangh = globVars.tangl+globVars.dtang
            
            
            globVars.zinl = globVars.zinh            
            globVars.zinh = D[int(globVars.tangh/globVars.dtang)]
        globVars.zin*rad
        
    if globVars.ielev == 1:
        if time <globVars.telvh:
            tratio = (time-globVars.telvl)/globVars.dtelv
            globVars.wlin = globVars.elvl+ (globVars.elvh-globVars.elvl)*tratio
        else:
            globVars.wlin = globVars.elvh            
            globVars.telvl = globVars.telvh
            globVars.telvh = globVars.telvl+globVars.dtelv         
            globVars.elvl = globVars.elvh            
            globVars.elvh = WL[int(globVars.telvh/globVars.dtelv)]
        sbeach.memory.dsurge = globVars.wlin    
    if globVars.iwind == 1:
        if time <globVars.twndh:
            tratio = (time-globVars.twndl)/globVars.dtwnd
            globVars.w = globVars.wl+ (globVars.wh-globVars.wl)*tratio
            globVars.zwind = globVars.zwndl+ (globVars.zwndh-globVars.zwndl)*tratio
        else:
            globVars.w = globVars.wh            
            globVars.twndl = globVars.twndh
            globVars.twndh = globVars.twndl+globVars.dtwnd
            
            globVars.elvl = globVars.elvh
            globVars.elvh = Wi[int(globVars.twndh/globVars.dtwnd),0]
            globVars.zwind = globVars.zwndh
            globVars.zwndl = globVars.zwndh
            globVars.zwndh = Wi[int(globVars.twndh/globVars.dtwnd),1]
        globVars.zwind*rad
            
def hbelev(sb,xB,yB):
    '''
    interpolates the harbottom elevations.
    '''
    for i in range(sb.arrsize.ndx):
        for j in range(len(xB)-1):
            if sb.memorymain.x[i]>=xB[j] and sb.memorymain.x[i]<xB[j+1]:
                break
            
        sb.memorymain.dbot[i] = yB[j] + (yB[j+1]-yB[j])*(sb.memorymain.x[i]-xB[j])/(xB[j+1]-xB[j])
#        print(i,j,sb.memorymain.dbot[i])
    #check not above profile
#    print('check')
    for i in range(sb.arrsize.ndx):
        if sb.memorymain.dbot[i]<sb.memorymain.di[i]:
            sb.memorymain.dbot[i] = sb.memorymain.di[i]
     
            
            
def initialize(sb):
    '''
    initialize arrays
    '''
    sb.memorymain.teta = np.zeros(sb.arrsize.ndx)
    sb.memorymain.ateta = np.zeros((sb.arrsize.ndx,4))
    
def input(sb,file,globVars):
    '''
    reads control file for SBEACH containing all input param values.
    sb is the SBEACH module
    '''
    rad = 0.01745 #convert degrees to radians
    
    with open(file,'r') as f:
        lineOld = ''
        idx=1
        for line in f.readlines():
            split = line.split()
#            print(split)
            if split[0][0] == 'A' or split[0][0] == 'B' or split[0][0] == 'C' or split[0][0] == 'D' or split[0][0] == 'E' or split[0][0] == 'F':
                lineOld = split[0]
                continue
            if lineOld == 'A.3':
                sb.arrsize.ndx = np.int(split[0])
                sb.memory.xstart = np.float(split[-1])
            if lineOld == 'A.4':
                idx = np.int(split[0])
            if idx == 1 and lineOld in ['A.5','A.6','A.7']:
                continue
            elif idx == 0 and lineOld in ['A.8','A.9']:
                continue
            if lineOld == 'A.8':
                ngrid = np.int(split[0])
                sb.memorymain.nrndxv = ngrid
                sb.memorymain.nrdxv = ngrid
            if lineOld == 'A.9':
                sb.memorymain.dxv = np.asarray(split[::2],dtype=np.float)
                sb.memorymain.ndxv = np.asarray(split[1::2],dtype=np.int)
            if lineOld == 'A.10':
                globVars.ndt = np.int(split[0])
                sb.memory.dt = np.float(split[-1])*60
            if lineOld == 'A.11':
                continue
            if lineOld == 'A.12':
                sb.memorymain.wri = np.asarray(split,dtype=np.int)
            if lineOld == 'A.13':
                continue
            if lineOld == 'A.14':
                globVars.elv1 = np.asarray(split,dtype=np.float)
#                elv2 = np.float(split[1])
#                elv3 = np.float(split[2])
                continue
            if lineOld == 'A.15':
                globVars.refelev = np.float(split[-1])
                continue
            if lineOld == 'A.16':
                sb.memory.k = np.float(split[0])
            if lineOld == 'A.17':
                sb.memory.eps = np.float(split[0])
            if lineOld == 'A.18':
                sb.memory.lamm = np.float(split[0])
            if lineOld == 'A.19':
                sb.tempc = np.float(split[0])
                
            if lineOld == 'B.2':
                globVars.iwave = np.int(split[0])
            if lineOld in ['B.3','B.4','B.5']:
                if globVars.iwave == 1:
                    continue
            elif lineOld in ['B.6']:
                 if globVars.iwave == 0:
                     continue
            if lineOld == 'B.6':
                globVars.dtwav = np.float(split[0])*60
            if lineOld == 'B.7':
                globVars.iang = np.float(split[0])
            if lineOld in ['B.8','B.9','B.10']:
                if globVars.iang ==1:
                    continue
            if lineOld in ['B.11']:
                if globVars.iang == 0:
#                    dtang=0
                    continue
            if lineOld == 'B.11':
                globVars.dtang = np.float(split[0])*60
            if lineOld == 'B.12':
                globVars.dmeas = np.float(split[0])
            if lineOld == 'B.13':
                globVars.irand = np.int(split[0])
            if lineOld in ['B.15']:
                if globVars.irand == 0:
                    continue
            if lineOld == 'B.15':
                globVars.iseed = np.int(split[0])
                globVars.rperc = np.int(split[1])/100
            if lineOld == 'B.16':
                globVars.ielev = np.int(split[0])
            if lineOld in ['B.17','B.18','B.19']:
                if globVars.ielev == 1:
                    continue
            if lineOld == 'B.20':
                globVars.dtelv = np.float(split[0])*60
            if lineOld == 'B.21':
                globVars.iwind = np.float(split[0])
            if lineOld in ['B.22','B.23','B.24']:
                if globVars.iwind == 1:
                    continue
            if lineOld in ['B.25']:
                if globVars.iwind == 0:
#                    globdtwnd=0
                    continue
            if lineOld == 'B.23':
                globVars.w = np.float(split[0])
                globVars.zwind = np.float(split[1])*rad
            if lineOld == 'B.25':
                globVars.dtwnd = np.float(split[0])*60
            if lineOld == 'C.1':
                globVars.tpin = np.int(split[0])
            if lineOld in ['C.2','C.3']:
                if globVars.tpin == 1:
                    continue
            if lineOld == 'C.4':
                globVars.dfs = np.float(split[0])
                
            if lineOld == 'C.5':
                sb.memory.d50 = np.float(split[0])/1000
            if lineOld == 'C.6':
#                print(split)
                sb.memory.bmax = np.float(split[0])*rad
                sb.memory.bi = 0.1
                
            if lineOld == 'D.1':
                ibchfill = np.int(split[0])
            if lineOld in ['D.2','D.3','D.4','D.5']:
                if ibchfill == 0:
                    continue
            
            if lineOld == 'E.1':
                iswall = np.int(split[0])
            if lineOld in ['E.2','E.3','E.4','E.5','E.6']:
                if iswall == 1:
                    continue
            
            if lineOld == 'F.1':
                globVars.ihbot = np.float(split[0])
            if lineOld in ['F.2']:
                if globVars.ihbot == 1:
                    globVars.scf = np.float(split[0])
            
#    return iwave,iang,ielev,iwind,irand,ndt,dtwav,dtang,dtelv,dtwnd, dmeas, iseed, rperc, w, zwind, dfs, elv1, refelev
                
                
                
                
                
                
                
                
                
                
                

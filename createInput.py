# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:48:40 2016

@author: justinpringle

Creates the input file for sbeach
parameters defined in SBEACH/Docs/CERC_93-2 User Manual report 3
"""
import numpy as np
file = 'input.txt'

'''gridParams'''
ndx = 252
xstart = -320/3.2808
idx = 1
dxc = 0 #constant grid
ngrid = 3 #variable grid number of regions
dxv = [5/3.2808,10/3.2808,20/3.2808] #dx widths converting from feet to metres
ndxv = [96,84,72] #number of grid cells in each region
res = [None]*(len(dxv)+len(ndxv))
res[::2] = dxv
res[1::2] = ndxv 
ndt = 1440
dt = 3
nwr = 3
wri = [360,720,960]
icomp=1
#track contours
elv1 = 12/3.2808
elv2 = 5/3.2808
elv3 = -5/3.2808
epd1 = 5/3.2808
epd2 = 3/3.2808
epd3 = 1/3.2808
refelev = 0
k = 1.5*10**-6#1.5*10**-6
eps = 0.002#0.002
lamm = 0.5
tempc = 15

'''waves and water levels'''
wvtype = 2
iwave = 1
hin = 0
t = 0
dtwav = 60
iang = 0
zin = 0
dtang = 0
dmeas = 33.8
irand = 1
iseed = 7878
rperc = 20
ielev = 1
telev = 0
dtelv = 60
iwind = 0
w = 0
zwind = 0
dtwind = 0

'''Beach'''
tpin = 1
dfs = 0.304
d50 = 0.35
bmax = 17.2

'''beach fill'''
ibchfill = 0
xbfs = 0
xbfe = 0
nfill = 0
xf = np.zeros(nfill)
efill = np.zeros(nfill)

'''seawall/ revetment'''
iswall = 0
xswall = 0
iswfail = 0
pefail = 0
wefail = 0
hfail = 0
ihbot = 1
scf=2


with open(file,'w') as f:
    f.write('A_________________________SBEACH model setup_____________________A\n')
    f.write('A.3 Total number of calc cells and position of landward boundary: ndx, xstart\n')
    f.write('\t%d\t%d\n'%(ndx,xstart))
    
    f.write('A.4 Grid type (0 = constant, 1 = variable): idx\n')
    f.write('\t%d\n'%idx)
    f.write('A.5 if grid is variable continue to A.8\n')
    
    f.write('A.6 Constant grid cell width: dxc\n')
    f.write('\t%d\n'%dxc)    
    f.write('A.7 if grid is constant continue to A.10\n')
    
    f.write('A.8 Number of different grid cell regions: ngrid\n')
    f.write('\t%d\n'%ngrid)
    
    f.write('A.9 Grid cell widths and number of cells in each region (landward to seaward) dxv[i],ndxv[i],i=1,ngrid\n')
    f.write(' '.join(map(str,res)))
    f.write('\n')
    
    f.write('A.10 Number of time steps and value of time steps in minutes: ndt,dt\n')
    f.write('\t%d %d\n'%(ndt,dt))
    
    f.write('A.11 Number of time steps for intermediate output: nwr\n')
    f.write('\t%d\n'%nwr)
    
    f.write('A.12 Time steps of intermediate output: wri[i],i=1,nwr\n')
    f.write(' '.join(map(str,wri)))
    f.write('\n')
    
    f.write('A.13 Is a measured profile available for comparison? 1=Yes, 2=No: icomp\n')
    f.write('\t%d\n'%icomp)
    
    f.write('A.14 Three profile elevation contours (maximum horizontal recession of each will be determined): elv1,elv2,elv3\n')
    f.write('\t%d %d %d\n'%(elv1,elv2,elv3))
    
    f.write('A.15 Three erosion profile depths and ref elevation: epd1, epd2, epd3, refelev\n')
    f.write('\t%d %d %d %d\n'%(epd1,epd2,epd3,refelev))
    
    f.write('A.16 Transport rate coefficient: k (m^4/N)\n')
    f.write('\t%g\n'%k)
    
    f.write('A.17 Coefficietn for slope dependent term: eps (m^2/s)\n')
    f.write('\t%g\n'%eps)
    
    f.write('A.18 Transport rate decay coefficient multiplier: lamm\n')
    f.write('\t%g\n'%lamm)
    
    f.write('A.19 Water temp in degrees c: tempc\n')
    f.write('\t%d\n'%tempc)
    
    f.write('B__________________________Waves,Water,Wind______________________B\n')
    
    f.write('B.1 Wave type (Monochromatic =1, var = 2): wvtype\n')
    f.write('\t%d\n'%wvtype)
    
    f.write('B.2 Wave height and period input (constant=0, variable=1) iwave\n')
    f.write('\t%d\n'%iwave)
    
    f.write('B.6 Time step of variable wave height and period (minutes)\n')
    f.write('\t%g\n'%dtwav)
    
    f.write('B.7 Wave anlge input (constant = 0, variable = 1): iang\n')
    f.write('\t%g\n'%iang)
    
    f.write('B.8 Constant wave angle: zin\n')
    f.write('\t%g\n'%zin)
    
    f.write('B.11 Time step of variable wave angle input in minutes:dtang\n')
    f.write('\t%g\n'%dtang)
    
    f.write('B.12 Water depth of input waves (deepwater=0): dmeas\n')
    f.write('\t%g\n'%dmeas)
    
    f.write('B.13 Is randomization of input wave height desired? irand\n')
    f.write('\t%d\n'%irand)
    
    f.write('B.15 Seed value for randomization and percent of variability: iseed, rperc\n')
    f.write('\t%d %d\n'%(iseed, rperc))
    
    f.write('B.16 Total water elevation input (constant =0, variable = 1) ielev\n')
    f.write('\t%d\n'%ielev)
    
    f.write('B.20 Time step of variale total water elevation input in minutes: dtelev\n')
    f.write('\t%g\n'%dtelv)
    
    f.write('B.21 Wind speed and angle input (constant = 0, variable = 1) iwind\n')
    f.write('\t%d\n'%iwind)
    
    f.write('B.23 Constant wind speed and angle: w, zwind\n')
    f.write('\t%g %g\n'%(w,zwind))
    
    f.write('B.25 Time step variable wind speed and angle\n')
    f.write('\t%g\n'%dtwind)
    
    f.write('C________________________Beach___________________________________C\n')
    
    f.write('C.1 Type of input profile (arbitrary =1, schematized =2): tpin\n')
    f.write('\t%d\n'%tpin)
    
    f.write('C.4 Depth corresponding to the landward end of the surfzone: dfs\n')
    f.write('\t%g\n'%dfs)
    
    f.write('C.5 Effective grain size diameter in mm: d50\n')
    f.write('\t%g\n'%d50)
    
    f.write('C.6 Maximum profile slope prior to avalanching in degrees: bmax\n')
    f.write('\t%g\n'%bmax)
    
    f.write('D________________________Beach Fill______________________________D\n')
    f.write('D.1 Is a beach fill present? (No=0, Yes=1): ibchfill\n')
    f.write('\t%d\n'%ibchfill)
    
    f.write('E________________________Seawall/Revetment_______________________E\n')
    f.write('E.1 Is a seawall present? (No=0, Yes = 1): iswall\n')
    f.write('\t%d\n'%iswall)
    
    f.write('F________________________Hard Bottom_______________________F\n')
    f.write('F.1 Is a hard bottom present? (No=0, Yes = 1): ihbot\n')
    f.write('\t%d\n'%ihbot)
    f.write('F.2 Dondrift scale factor: scf\n')
    f.write('\t%d\n'%scf)
    
    
    

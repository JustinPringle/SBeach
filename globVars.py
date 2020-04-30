# -*- coding: utf-8 -*-
"""
Created on Mon May  9 12:19:35 2016

@author: justinpringle
"""
import numpy as np

twavl, twavh, tangl, tangh, telvl, telvh, twndl,twndh = 0,0,0,0,0,0,0,0
hinl, hinh, tl,th, zinl,zinh, elvl,elvh, wl,wh,zwndl,zwndh = 0,0,0,0,0,0,0,0,0,0,0,0

iwave =0
iang =0
ielev = 0
iwind =0
irand =0
ndt =0
dtwav =0
dtang =0 
dtelv =0
dtwnd =0
dmeas =0
iseed =0
rperc =0
w =0
zwind =0 
dfs = 0
elv1 = np.zeros(3)
xrefelv = np.zeros(3)
xelvi = np.zeros(3)
xelvr = np.zeros(3)
ielv = np.zeros(3)
refelev = 0
factl=1
iswlim = np.asarray([0])
#seawall
nswall = np.asarray([0])

iswall = 0
iswfail=0
swfail=0
xswall = 0
pefail = 0
wefail=0
hfail=0
#Hard bottom
ihbot = np.array([0])
scf = np.array([0])
#waves
hin = 0
tin = 0
L0 = 0
zin = 0
wlin  = 0
w = 0
zwind = 0
istop = np.asarray([0])
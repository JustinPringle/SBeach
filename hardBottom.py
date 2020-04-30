# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:23:01 2016

@author: justinpringle

create hardBottom file
"""
gPath = '/Volumes/NO NAME/SBEACH/'

import numpy as np

hx = [50,100]
hy = [2,-2]
with open('%sTest/reach3.hb'%gPath,'w') as f:
    for x,y in zip(hx,hy):
        f.write('%f\t%f\n'%(x,y))

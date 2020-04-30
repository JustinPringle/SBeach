# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 10:39:14 2016

@author: justinpringle
note: profile is positive below MSW
"""
gPath = ''
import model as sb
import numpy as np
import Subroutines as sub
import globVars
import matplotlib.pyplot as plt
import matplotlib.animation as animate
import os

def init():
    line.set_data([],[])
    line2.set_data([],[])
    line3.set_data([],[])
    line4.set_data([],[])
    return line,line2,line3,line4,
    
def aniFunc(i):
    plt.title('Hour = %s'%i)
    ind = np.where(hDict[i]>0)[0][0]
    ind2 = np.where(np.abs(wlDict[i][1::]-wlDict[i][0:-1])>0)[0][0]
#    ind3 = np.where(qDict[i][:]>0)[0][0]
#    print(ind2)
#    ind = 0
    line.set_data(x,sDict[i])
    line2.set_data(x[ind::],hDict[i][ind:-1])#+nDict[i][ind-1::]+shoreDict[i][ind-1::])
##    
    line3.set_data(x[ind2+1::],wlDict[i][ind2+1:-1])#+WL[i])
    line4.set_data(x,qDict[i][0:-1])#dDict[i]/100)
    return line,line2,line3,line4,
    
if __name__ == '__main__':
    Plot=True
    '''Initialize'''
#read control file and allocate params to memory
    file = 'input.txt'
    sub.input(sb,file,globVars)
#allocate main arrays to memory
#sub.initialize(sb)
    sb.memory.initmemory()
    sb.memorymain.allocatearr()
    sub.allocate(sb)
    sb.memoryq_u.rrm = sb.memory.k*100
#
##initial profile
    F = '%sTest/'%gPath
    xIn = np.asarray(np.loadtxt('%sreach3.xz'%F)[:,0],dtype=np.float64,order='F')/3.2808
    yIn = np.asarray(np.loadtxt('%sreach3.xz'%F)[:,1],dtype=np.float64,order='F')/3.2808
    xF = np.loadtxt('%sfinal.xz'%F)[:,0]/3.2808
    yF = np.loadtxt('%sfinal.xz'%F)[:,1]/3.2808
    xS = np.loadtxt('%sreach3_Final.xz'%F)[:,0]/3.2808
    yS = np.loadtxt('%sreach3_Final.xz'%F)[:,1]/3.2808
    
# Hard Bottoms
    if globVars.ihbot == 1:
        Fhb = '%sTest/'%gPath
        xHB = np.asarray(np.loadtxt('%sreach3.hb'%Fhb)[:,0],dtype = np.float64,order='F')
        yHB = np.asarray(np.loadtxt('%sreach3.hb'%Fhb)[:,1],dtype = np.float64,order='F')

#Load waves, and water levels etc
    H = np.asarray(np.loadtxt('%sreach3.hs'%F)[:,1],dtype=np.float64,order='F')/3.2808 #4.65#
#    H = np.asarray([10 for i in range(len(H2))])
    T = np.asarray(np.loadtxt('%sreach3.per'%F)[:,1],dtype=np.float64,order='F')
#    D = np.asarray(np.loadtxt('%sreach3.hs'%F)[:,1],dtype=np.float64,order='F')
    WL = np.asarray(np.loadtxt('%sreach3.wl'%F)[:,1],dtype=np.float64,order='F')/3.2808
#    WL = np.asarray([0 for i in range(len(WL2))])
    if globVars.iwind == 0:
        Wi = np.asarray([globVars.w,globVars.zwind])
    else:
        #load variable wave angles
        Wi = 0
    if globVars.iang == 0:
        D = np.asarray([0])
    else:
        D = np.asarray(np.loadtxt('%sreach3.hs'%F)[:,1],dtype=np.float64,order='F')
#Assign grid cell widths to each cell        
    sb.memory.dx = sb.vargrid(sb.memorymain.dxv,sb.memorymain.ndxv,sb.arrsize.ndx)
    
#Initial Profile - + below SWL and - Above SWl
    sb.initialprofile(globVars.tpin,-1*yIn,xIn,sb.memory.dx,sb.arrsize.ndx,len(xIn))
    if globVars.ihbot == 1:
        sub.hbelev(sb,xHB,-yHB)
        
#Contour elvs to track - choose 3
    globVars.ielv,globVars.xelvi,globVars.xelvr,globVars.xrefelv,nrmax,erroc = sb.initialcontourlocations(globVars.elv1,globVars.refelev)
#Initialize som arrays
    sb.memorymain.d = sb.memorymain.di
    sb.memorymain.dp = sb.memorymain.di
#    print(sb.memorymain.d[0],sb.memorymain.di[0])
##calculate settling velocity
    sb.kvisc(sb.memorymain.teta,sb.memorymain.ateta)
    vf = sb.falvel(sb.tempc,sb.memorymain.teta,sb.memorymain.ateta)
#
##calc hnsratio
    hnsr = np.zeros(249)
    sb.hnsratio(hnsr[12:249])
    
    #Loop through time NOTE: dt and ndt input are in minutes and then converted to seconds
    ndt=globVars.ndt
    sub.temporalInfo(H,D,T,WL,Wi,globVars)
    
    hDict = []
    qDict=[]
    sDict=[]
    wlDict =[]
    x = np.asarray([i for i in sb.memorymain.x])
    nru =np.asarray([0])
    nfs = np.asarray([0])
    imin = np.asarray([0])
    blr=np.asarray([0])
    izero=np.asarray([0])
    lru=np.asarray([0])
    drumax=np.asarray([0],dtype=np.float)
    dru=np.asarray([0],dtype=np.float)
    idrumax=np.asarray([0])
    overwash = np.asarray([0])
#    sb.memory.ceq=0
    
    for i in range(ndt):
#        if i >440 and i<450:
#            print(i)
        time = i*sb.memory.dt
        '''read wave arrays and interpolate between time points'''
        sub.readWaves(time,H,T,D,WL,Wi,globVars,sb)
        
        '''preprocess wave information (deep water boundary conditions)'''
#        globVars.hin*=np.random.uniform(low=0.8,high=1.2)
        hbeg,h0,zbeg = sb.preprocesswaves(globVars.iwave,globVars.irand,
                                         0,globVars.tin,globVars.dmeas,
                                         globVars.rperc,globVars.hin,sb.memorymain.d[-1],
                                         globVars.zin)
#        print('python',i,hbeg)                                 
        '''add water level to depth'''
#        if i == 0:
#            sb.memorymain.dp = sb.memorymain.d
        sb.memorymain.d+=sb.memory.dsurge
        sb.memorymain.dp+=sb.memory.psurge
        
        '''smooth profile'''                                
        domsmo = sb.domsmo(sb.memory.dx,sb.memorymain.d,h0,globVars.nswall,sb.arrsize.ndx)
#        print(domsmo[0])  
#        '''Transform waves from deep water to shore'''
##        sb.memorymain.h = np.zeros((sb.arrsize.ndx+1))
        sb.wavran(sb.memory.dx,globVars.nswall,hbeg,globVars.tin,zbeg,domsmo,
                  sb.memorymain.h,sb.memorymain.e,globVars.istop,sb.memorymain.zw,sb.memorymain.diss,
                  sb.memorymain.brok,sb.memorymain.x,globVars.iswlim,hnsr[12:249],globVars.w,globVars.zwind,
                  sb.memorymain.d,sb.arrsize.ndx)
        
#        print('Completed')    
        '''Add the setup to depth before calculating transport rates'''
        sb.addsetuptodepth(globVars.istop,domsmo,sb.memorymain.e,sb.arrsize.ndx)
#        
#        '''calculae the slope term in transport rate calcs'''
        sb.memorymain.edhdx = np.zeros(sb.arrsize.ndx)
        sb.slopetransport(globVars.istop,globVars.nswall,sb.memorymain.edhdx,sb.memorymain.d,
                          sb.memory.dx,sb.arrsize.ndx)
#                          
#        '''Cross shore transport rates in the surf zone and offshore zone'''
        sb.trran(h0,globVars.zin,globVars.tin,vf,sb.memorymain.x,domsmo,globVars.istop,
                 sb.memorymain.zw,sb.memorymain.diss,sb.memorymain.brok,
                 sb.memorymain.edhdx,sb.memorymain.q,globVars.iswlim,
                 sb.arrsize.ndx)
#                 
#        '''Cross shore transport rates in the swash zone'''
#        
#        
##        istop = 0
        sb.swshdom(sb.memorymain.x,sb.memorymain.d,sb.memorymain.e,h0,globVars.nswall,
                   globVars.dfs,nru,nfs,imin,blr,globVars.iswall,izero,lru,globVars.zin,
                   drumax,i,idrumax,sb.memorymain.brok,dru,sb.arrsize.ndx)
#        if nru<nfs and nru<imin and sb.memorymain.q[nfs]>0:
#            print(i,nru,imin,x[nru],x[imin],sb.memorymain.d[nru],sb.memorymain.d[imin])
#        if i>700 and i<750:
#            print(i,nru,-dru+sb.memory.dsurge)
#        '''Overwash is calculated here'''
        sb.exrun(nru,nfs,sb.memorymain.q,imin,sb.memorymain.x,sb.memorymain.d,
                 overwash,izero,lru,blr,globVars.dfs,dru,sb.arrsize.ndx)
#                 
#        '''Call conservation of sand and update shoreline'''
        if i ==0:
            sb.memorymain.qp=sb.memorymain.q
        
        sb.consan(nru,sb.memorymain.d,sb.memorymain.q,sb.memorymain.dp,
                  sb.memorymain.qp,globVars.nswall,globVars.pefail,
                  globVars.iswfail,globVars.swfail,time,globVars.istop,globVars.iswall,
                  globVars.factl,'m/s',globVars.ihbot,sb.memorymain.ihb,
                  sb.memorymain.dbot,globVars.scf,sb.arrsize.ndx)
        
#        sb.memorymain.d -=sb.memory.dsurge
#        sb.memorymain.dp-=sb.memory.psurge
#        '''store data'''
        if i%10==0:
            with open('outFiles/Wave/%d.txt'%i,'w') as f:
                f.write('\n'.join(list(map(str,sb.memorymain.h))))
            with open('outFiles/Shore/%d.txt'%i,'w') as f:
                f.write('\n'.join(list(map(str,-sb.memorymain.d))))
            with open('outFiles/Q/%d.txt'%i,'w') as f:
                f.write('\n'.join(list(map(str,sb.memorymain.q*10000))))
            with open('outFiles/Water/%d.txt'%i,'w') as f:
                f.write('\n'.join(list(map(str,sb.memorymain.e + sb.memory.dsurge))))
#            hDict.append(np.array([i for i in sb.memorymain.h]))
#            qDict.append(np.array([i*10000 for i in sb.memorymain.q]))
#            sDict.append(np.array([i for i in sb.memorymain.d]))
#            wlDict.append(np.array([i+0*sb.memory.dsurge for i in sb.memorymain.e]))
        
        if np.isnan(sb.memorymain.h.any()):
            print(i)
#        print(i,np.isnan(sb.memorymain.h.any()))
          
        sb.memorymain.q=np.zeros(len(sb.memorymain.q))
        sb.memory.psurge=sb.memory.dsurge
                   
#        hDict.append(np.ma.masked_array([i for i in sb.memorymain.h],0))
    
#    Clear memory    
    if not Plot:
        sb.memorymain.donememorymain()
        sb.memory.donememory()
        sb.memoryq_u.donememoryq_u()

    '''Plot and animate'''
    if Plot:
        hDict = [np.loadtxt('outFiles/Wave/%d.txt'%(i)) for i in range(0,ndt)  if i%10==0]
        sDict = [np.loadtxt('outFiles/Shore/%d.txt'%(i)) for i in range(0,ndt)  if i%10==0]
        qDict = [np.loadtxt('outFiles/Q/%d.txt'%(i)) for i in range(0,ndt) if i%10==0]
        wlDict = [np.loadtxt('outFiles/Water/%d.txt'%(i)) for i in range(0,ndt)  if i%10==0]
        fig = plt.figure(figsize=(15,10))#,dpi=300)
        ax = fig.add_subplot(111)
        
        ax.plot(x,sDict[0],color='r',alpha = 0.5)
        ax.plot(xF,yF,color='b',lw=3)
        ax.plot(xS,yS,color='r',lw=3)
        ax.plot(sb.memorymain.x,-sb.memorymain.dbot,color='gray')
        ax.set_xlim([-50,600])#[60,1000])
        ax.set_ylim([-10,5])
#        if Test:
#            ax.plot(xF,yF,'b',lw=3)
#            ax.plot(xS,yS,'r',lw=3)
#        ax.set_aspect('equal')
        line, = ax.plot([],[],lw=1,color='k')
        line2, = ax.plot([],[],lw=1,color='b')
        line3, = ax.plot([],[],lw=1,color = 'g')
        line4, = ax.plot([],[],lw=1,color='r',alpha=0.5)
        
        anim = animate.FuncAnimation(fig,aniFunc,init_func=init,frames=len(hDict),interval = 1, blit=True)
        plt.rcParams['animation.ffmpeg_path']='/Volumes/NO NAME/SBEACH/ffmpeg'
#        ffwriter = animate.FFMpegWriter(fps=10)
#        anim.save('test.avi',writer=ffwriter,bitrate=1000,extra_args=['-vcodec','libx24'])
        sb.memorymain.donememorymain()
        sb.memory.donememory()
        sb.memoryq_u.donememoryq_u()


################################################################################
## Call memory modules
#model.memory.initmemory
#
#model.memorymain.nrdxv = 3
#model.memorymain.ndxv =np.asarray([96,84,72],dtype = np.int32,order='F')
#model.memorymain.dxv = np.asarray([5,10,20],dtype = np.float64,order='F')
#model.arrsize.ndx = int(np.sum(model.memorymain.ndxv))
#
#F = '%sTest/'%gPath
#xIn = np.asarray(np.loadtxt('%sreach3.xz'%F)[:,0],dtype=np.float64,order='F')
#yIn = np.asarray(np.loadtxt('%sreach3.xz'%F)[:,1],dtype=np.float64,order='F')
#
#model.memory.xstart = -320
#assert model.memory.xstart >= xIn[0]
#model.memory.dx = model.vargrid(model.memorymain.dxv,model.memorymain.ndxv,model.arrsize.ndx)
#
#Y,X = model.initialprofile(1,-1*yIn,xIn,model.arrsize.ndx,len(xIn)-1)
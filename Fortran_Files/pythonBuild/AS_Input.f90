# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:12:56 2016

@author: justinpringle
"""

!***************************************************************************************
     Subroutine ReadTemporalInfo(IANG,IWAVE,IWIND,OFFSHORE,FACTH,DTR,TWAVL,TWAVH,&
    	HINL,TL,HINH,TH,TANGL,TANGH,ZINL,ZINH,TELVL,TELVH,DSURGEL,DSURGEH,TWNDL,&
    	TWNDH,WL,ZWINDL,WH,ZWINDH,TRATIO,T,TIME,ZIN,W,ZWIND,HIN,DTWAV,DTANG,DTELV,DTWND,FACTL)
    	Use Memory, ONLY: DSURGE,IELEV,LO,PI
    	Use IOHandles, ONLY: fhiWAV,fhiANG,fhiWND,fhiELV
    	Implicit NONE
    	Integer,intent(inout):: IANG,IWAVE,IWIND
    	Logical,intent(inout):: OFFSHORE
    	Real(8),intent(inout):: FACTH,DTR,TWAVL,TWAVH,HINL,TL,HINH,TH,TANGL,TANGH,ZINL,&
    		ZINH,TELVL,TELVH,DSURGEL,DSURGEH,TWNDL,TWNDH,WL,ZWINDL,WH,ZWINDH,TRATIO,&
    		T,TIME,ZIN,W,ZWIND,HIN,DTWAV,DTANG,DTELV,DTWND,FACTL
    	Logical:: CycleLoop
	
    	IF(IWAVE==1) THEN
    		CycleLoop=.true.
    		Do While(CycleLoop)
    			IF(TIME<=TWAVH) THEN
    				TRATIO=(TIME-TWAVL)/DTWAV
    				HIN=HINL+(HINH-HINL)*TRATIO
    				T=TL+(TH-TL)*TRATIO
    				CycleLoop=.false.
    			ELSE
    				TWAVL=TWAVH
    				TWAVH=TWAVL+DTWAV
    				HINL=HINH
    				TL=TH
    				READ(fhiWAV,*) HINH,TH
    			End If
    		End Do
    		HIN=HIN*FACTL*facth
    		LO=1.5613*T**2
    	End If
    	IF(IANG==1) THEN
    		CycleLoop=.true.
    		Do While(CycleLoop)
    			IF(TIME<=TANGH) THEN
    				TRATIO=(TIME-TANGL)/DTANG
    				ZIN=ZINL+(ZINH-ZINL)*TRATIO
    				CycleLoop=.false.
    			ELSE
    				TANGL=TANGH
    				TANGH=TANGL+DTANG
    				ZINL=ZINH
    				READ(fhiANG,*) ZINH
    			End If
    		End Do
    		IF(ZIN<-90.OR.ZIN>90) THEN
    			OFFSHORE=.TRUE.
    		ELSE
    			OFFSHORE=.FALSE.
    		End If
    		ZIN=ZIN*DTR
    	End If
    	IF(IELEV==1) THEN
    		CycleLoop=.true.
    		Do While(CycleLoop)
    			IF(TIME<=TELVH) THEN
    				TRATIO=(TIME-TELVL)/DTELV
    				DSURGE=DSURGEL+(DSURGEH-DSURGEL)*TRATIO
    				CycleLoop=.false.
    			ELSE
    				TELVL=TELVH
    				TELVH=TELVL+DTELV
    				DSURGEL=DSURGEH
    				READ(fhiELV,*) DSURGEH
    			End If
    		End Do
    		DSURGE=DSURGE*FACTL
    	End If
    	IF(IWIND==1) THEN
    		CycleLoop=.true.
    		Do While(CycleLoop)
    			IF(TIME<=TWNDH) THEN
    				TRATIO=(TIME-TWNDL)/DTWND
    				W=WL+(WH-WL)*TRATIO
    				ZWIND=ZWINDL+(ZWINDH-ZWINDL)*TRATIO
    				CycleLoop=.false.
    			ELSE
    				TWNDL=TWNDH
    				TWNDH=TWNDL+DTWND
    				WL=WH
    				ZWINDL=ZWINDH
    				READ(fhiWND,*) WH,ZWINDH
    			End If
    		End Do
    		W=W*FACTL
    		ZWIND=ZWIND*atan(1.d0)/4.5d1 !i.e., pi/180 ~ 0.0175
    	End If
    
    	Return
     End Subroutine
! Format constants, QuantityUnits and Timing modules
!***************************************************************************************
Module FormatConst
	Character(Len=*),parameter:: &
		FM7000='(10F8.2)',&
		FM7050='(10F8.1)',&
		FM8000='(10I8)'	
End Module FormatConst
!***************************************************************************************
Module QuantityUnits
	Implicit NONE
	Integer,parameter::StrQuanUnitLen=20
	Integer::UNITS
	Real(8)::FACTL,FACTV
	Character(Len=StrQuanUnitLen)::LUNITS,L2UNITS,VUNITS
	
	Contains
	Subroutine SetQuantityUnits
		use Memoryq_u, ONLY: ERROCC
		Implicit NONE
		IF(UNITS==1) THEN
			FACTL=1.0
			FACTV=1.0
			LUNITS='(m)'
			L2UNITS='(m^2)'
			VUNITS='(m^3/m)'
		ELSEIF(UNITS==2) THEN
			FACTL=0.3048
			FACTV=0.3048**2.0*9.*3.
			LUNITS='(ft)'
			L2UNITS='(ft^2)'
			VUNITS='(yd^3/ft)'
		ELSE
			Call WriUFLOG('ERROR: INVALID UNITS SPECIFICATION')
			ERROCC=.True.
		End If
		
		Return
	End Subroutine SetQuantityUnits

	Subroutine AdjustQuantityUnits(IBFILL,IDX,NGRID,NFILL,FACTH,DTR,DFS,XF,EFILL,DMEAS,&
		RPERC,ZIN,DXC,W,ZWIND,PEFAIL,WEFAIL,HFAIL,EDP,HIN,DTWAV,DTANG,TELEV,DTELV,&
		DTWND,XLAND,DLAND,XLBDUNE,DLBDUNE,XLCDUNE,DLCDUNE,XSCDUNE,DSCDUNE,XBERMS,DBERMS,&
		XBERME,DBERME,XFORS,DFORS,XBFS,XBFE,XSWALL,FACTL,REFELV,ELV,OFFSHORE,DXV)
		Use Memory, ONLY: BAV,BMAX,D50,DT,XSTART
		Use MemoryMain, ONLY: NrDXV
		Implicit NONE
		Integer,intent(inout):: IBFILL,IDX,NGRID,NFILL
		Real(8),intent(inout):: FACTH,DTR,DFS,XF(10),EFILL(10),DMEAS,&
			RPERC,ZIN,DXC,W,ZWIND,PEFAIL,WEFAIL,HFAIL,EDP(3),HIN,DTWAV,DTANG,&
			TELEV,DTELV,DTWND,XLAND,DLAND,XLBDUNE,DLBDUNE,XLCDUNE,DLCDUNE,XSCDUNE,&
			DSCDUNE,XBERMS,DBERMS,XBERME,DBERME,XFORS,DFORS,XBFS,XBFE,XSWALL,FACTL,&
			REFELV,ELV(3)
		Logical,intent(inout)::OFFSHORE
		Real(8),intent(inout):: DXV(NrDXV)
		Integer:: I
		
		XSTART=XSTART*FACTL
		DXC=DXC*FACTL
		IF(IDX==1) THEN
			DO I=1,NGRID
				DXV(I)=DXV(I)*FACTL
			End Do
		End If
		ELV=-ELV*FACTL
		EDP=EDP*FACTL
		REFELV=-REFELV*FACTL
		HIN=HIN*FACTL*facth
		DMEAS=DMEAS*FACTL
		TELEV=TELEV*FACTL
		W=W*FACTL
		XLAND=XLAND*FACTL
		DLAND=-DLAND*FACTL
		XLBDUNE=XLBDUNE*FACTL
		DLBDUNE=-DLBDUNE*FACTL
		XLCDUNE=XLCDUNE*FACTL
		DLCDUNE=-DLCDUNE*FACTL
		XSCDUNE=XSCDUNE*FACTL
		DSCDUNE=-DSCDUNE*FACTL
		XBERMS=XBERMS*FACTL
		DBERMS=-DBERMS*FACTL
		XBERME=XBERME*FACTL
		DBERME=-DBERME*FACTL
		XFORS=XFORS*FACTL
		DFORS=-DFORS*FACTL
		DFS=DFS*FACTL
		XBFS=XBFS*FACTL
		XBFE=XBFE*FACTL
		IF(IBFILL==1) THEN
			DO I=1,NFILL
				XF(I)=XF(I)*FACTL
				EFILL(I)=-EFILL(I)*FACTL
			End Do
		End If
		XSWALL=XSWALL*FACTL
		PEFAIL=PEFAIL*FACTL
		WEFAIL=WEFAIL*FACTL
		HFAIL=HFAIL*FACTL
		!*
		!* CONVERT D50 FROM MM TO METERS
		!*
		D50=D50/1000.
		!*
		!* CONVERT ANGLES FROM DEGREES TO RADIANS
		!*
		DTR=0.01745
		IF(ZIN<-90.OR.ZIN>90) THEN
			OFFSHORE=.TRUE.
		ELSE
			OFFSHORE=.FALSE.
		End If
		ZIN=ZIN*DTR
		ZWIND=ZWIND*DTR
		BMAX=BMAX*DTR
		BAV=BAV*DTR
		!*  
		!* CONVERT TIME STEPS FROM MINUTES TO SECONDS
		!*
		DT=DT*60.
		DTWAV=DTWAV*60.
		DTANG=DTANG*60.
		DTELV=DTELV*60.
		DTWND=DTWND*60.
		!*
		!* CONVERT RPERC FROM PERCENTAGE TO FRACTION
		!*
		RPERC=RPERC/100.
		Return
	End Subroutine
End Module QuantityUnits
!***************************************************************************************
Module Timing
	Use IFPORT, ONLY: RTC
	Implicit NONE
	Real(8):: StartTime

	Contains
	Subroutine InitTiming
		Implicit NONE
		StartTime=rtc()
		return
	End Subroutine InitTiming
	
	Subroutine TimeTake(TimeTaken,TimeUnits)
		Implicit NONE
		Real(8),intent(inout):: TimeTaken
		Character(Len=8),intent(inout):: TimeUnits
		
		TimeTaken=(rtc()-StartTime)/3.6d3
		IF(TimeTaken>=1.) THEN
			TimeUnits='hours.'
		ELSEIF(TimeTaken*6.d1>1.) THEN
			TimeTaken=TimeTaken*6.d1
			TimeUnits='minutes.'
		ELSE
			TimeTaken=TimeTaken*3.6d3
			TimeUnits='seconds.'
		End If
		Return
	End Subroutine TimeTake
End Module Timing
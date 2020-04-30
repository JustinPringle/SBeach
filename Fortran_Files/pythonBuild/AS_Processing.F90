!***************************************************************************************
Subroutine InitialContourLocations(IELV,ELV,XELVI,XELVR,NRMAX,ErrOcc,REFELV,XREFELV)
!* FIND INITIAL LOCATION OF 3 ELEVATIONS TO MONITOR AND REFERENCE ELEVATION
	Use ArrSize,ONLY: NDX,StrLenMax
     !f2py integer, intent(aux):: NDX, StrLenMax
	Use MemoryMain,ONLY: DI,X
     !f2py real(8), intent(aux):: DI,X
!	Use QuantityUnits,ONLY: FACTL,LUNITS
	Implicit NONE
	Integer,intent(out),dimension(0:2):: IELV
	Real(8),intent(in), dimension(0:2):: ELV
      Real(8), intent(in)::REFELV
	Real(8),intent(out), dimension(0:2):: XELVI,XELVR
      Real(8),intent(out)::XREFELV
	Integer,intent(out), dimension(0:2):: NRMAX
	Logical,intent(out):: ErrOcc
	
	Logical:: CycleLoop
	Integer::I,J
	Character(Len=StrLenMax):: TmpStr
	Real(8):: FACTL
      FACTL=1
	IELV(0:2)=1
	IF(ELV(1)==ELV(0)) IELV(1)=0  
	IF(ELV(2)==ELV(1).OR.ELV(2)==ELV(0)) IELV(2)=0
	DO I=0,2
		IF(IELV(I)==0) Cycle
		CycleLoop=.false.
		DO J=NDX-1,1,-1
			IF(DI(J)>=ELV(I).AND.DI(J-1)<=ELV(I)) THEN
                      !write(*,*) DI(0), DI(J), DI(J-1), J
				XELVI(I)=X(J-1)+(X(J)-X(J-1))/(DI(J)-DI(J-1))*(ELV(I)-DI(J-1))
				XELVR(I)=XELVI(I)
				NRMAX(I)=J
				CycleLoop=.TRUE.
				Exit
			End If
		End Do
		if(CycleLoop) Cycle
		Write(*,'(" ERROR: THE ELEVATION: ",F0.2," "," DOES NOT EXIST&
			& ON THE INITIAL PROFILE.")') -ELV(I)/FACTL
!		Call WriUFLOG(TmpStr)
		Return
	End Do
	ErrOcc=.TRUE.
	DO J=NDX-1,1,-1
		IF(DI(J)>=REFELV.AND.DI(J-1)<=REFELV) THEN
			XREFELV=X(J-1)+(X(J)-X(J-1))/(DI(J)-DI(J-1))*(REFELV-DI(J-1))
			ErrOcc=.FALSE.
			Exit
		End If
	End Do
	if(ErrOcc)Then
		Write(*,'(" ERROR: THE ELEVATION: ",F0.2," "," DOES NOT EXIST&
			& ON THE INITIAL PROFILE.")') -REFELV/FACTL
!		Call WriUFLOG(TmpStr)
		Return
	End If
	return
End Subroutine
!***************************************************************************************
!Subroutine PostProcessingResults(NDX,ISTMIN,RunSimCode,D,DX,DI,VOLCH,X,NRSV,IELV,&
!	ICOMP,NEDP1,NEDP2,NEDP3,XREFELV,DRUMAX,XELVI,XELVR,EDP,ELV,SUMQ)
!	Use QuantityUnits
!	Use Memoryq_u,ONLY: ErrOcc 
!	Use MemoryMain, ONLY: DMX,TDMX,DMN,TDMN,HMX,THMX,EMX,TEMX,DPMX,TDPMX,DSV,XSV
!	Use Memory, ONLY: MassError
!	Use IOHandles, ONLY: fhoPRC,fhoXVR,fhoRPT,fhoLOG
!	Use FormatConst, ONLY: FM7000,FM7050,FM8000
!	Implicit NONE
!	Integer,intent(in):: NDX,ISTMIN,NRSV,IELV(3),ICOMP
!	Integer,intent(inout):: RunSimCode,NEDP1,NEDP2,NEDP3
!	Real(8),intent(inout):: D(NDX),DX(NDX),DI(NDX),VOLCH(NDX),X(NDX)
!	Real(8),intent(in):: XREFELV,DRUMAX,XELVI(3),XELVR(3),EDP(3),ELV(3),SUMQ
!	
!	
!	!local variables
!	Integer:: J,IEMXMAX
!	Real(8):: SUMVOL,SUMD,SUMDI,DISTEDP1,DISTEDP2,EMXMAX,DISTEDP3,XXS,XXE
!	!local formats
!	Character(len=*),parameter:: &
!		FM1910='(A7," CALCULATED ELEVATION ",A," AT GIVEN POSITION ON PROFILE")',&
!		FM6000='(" A ",F0.2," ",A," EROSION DEPTH DID NOT OCCUR ANYWHERE ON THE PROFILE.")',&
!		FM6005='(" POSITION OF LANDWARD MOST OCCURRENCE OF A ",F0.2," ",A," EROSION DEPTH:"&
!			&,/,F0.1," ",A)',&
!		FM6006='(" DISTANCE FROM POSITION OF REFERENCE ELEVATION ON INITIAL PROFILE",/,&
!			&" TO POSITION OF LANDWARD MOST OCCURRENCE OF A ",F0.2," ",A," EROSION DEPTH:"&
!			&,/,F0.1," ",A)'
!	Integer:: I
!			
!
!
!	!*
!	!* WRITE MAX AND MIN CALCULATED ELEVATIONS AT EACH POSITION ON PROFILE
!	!*
!	CALL DSHLN(fhoPRC)
!	WRITE(fhoPRC,FM1910)'MAXIMUM',Trim(LUNITS)
!	WRITE(fhoPRC,FM7000)(-DMX(J)/FACTL,J=1,NDX)
!	WRITE(fhoPRC,*) 'CORRESPONDING TIME STEPS'
!	WRITE(fhoPRC,FM8000)(TDMX(J),J=1,NDX)
!	CALL DSHLN(fhoPRC)
!	WRITE(fhoPRC,FM1910)'MINIMUM',Trim(LUNITS)
!	WRITE(fhoPRC,FM7000)(-DMN(J)/FACTL,J=1,NDX)
!	WRITE(fhoPRC,*) 'CORRESPONDING TIME STEPS'
!	WRITE(fhoPRC,FM8000)(TDMN(J),J=1,NDX)
!
!	CALL DSHLN(fhoXVR)
!	WRITE(fhoXVR,'("MAXIMUM CALCULATED WAVE HEIGHT ",A," AT GIVEN POSITION ON PROFILE&
!			&FIRST VALUE CORRESPONDS TO CELL ",I0)') Trim(LUNITS),ISTMIN
!	WRITE(fhoXVR,FM7000)(HMX(J)/FACTL,J=ISTMIN,NDX)
!	WRITE(fhoXVR,*) 'CORRESPONDING TIME STEPS'
!	WRITE(fhoXVR,FM8000)(THMX(J),J=ISTMIN,NDX)
!	CALL DSHLN(fhoXVR)
!	WRITE(fhoXVR,'("MAXIMUM CALCULATED TOTAL WATER ELEVATION + SETUP ",A,&
!			&/," AT GIVEN POSITION ON PROFILE",/,&
!			&"FIRST VALUE CORRESPONDS TO CELL ",I0)') Trim(LUNITS),ISTMIN
!	WRITE(fhoXVR,FM7000)(EMX(J)/FACTL,J=ISTMIN,NDX)
!	WRITE(fhoXVR,*) 'CORRESPONDING TIME STEPS'
!	WRITE(fhoXVR,FM8000)(TEMX(J),J=ISTMIN,NDX)
!	CALL DSHLN(fhoXVR)
!	WRITE(fhoXVR,'("MAXIMUM CALCULATED WATER DEPTH ",A," AT GIVEN POSITION ON PROFILE",/,&
!			&"FIRST VALUE CORRESPONDS TO CELL ",I0)') Trim(LUNITS),ISTMIN
!	WRITE(fhoXVR,FM7000)(DPMX(J)/FACTL,J=ISTMIN,NDX)
!	WRITE(fhoXVR,*) 'CORRESPONDING TIME STEPS'
!	WRITE(fhoXVR,FM8000)(TDPMX(J),J=ISTMIN,NDX)
!
!
!
!	SUMVOL=0.0
!	DO J=1,NDX
!		SUMVOL=SUMVOL+(DI(J)-D(J))*DX(J)
!		VOLCH(J)=SUMVOL
!	End Do
!	CALL DSHLN(fhoXVR)
!	WRITE(fhoXVR,'("CALCULATED CHANGE IN VOLUME LANDWARD OF GIVEN",/,&
!		&"POSITION ON PROFILE ",A)') Trim(VUNITS)
!	WRITE(fhoXVR,FM7050)(VOLCH(J)/FACTV,J=1,NDX)
!
!	!*
!	!* CHECK MASS CONSERVATION BY COMPARING THE AMOUNT OF SAND IN THE
!	!* INITIAL PROFILE AND THE FINAL PROFILE. THE CUBIC SPLINE ROUTINE
!	!* IS USED TO APPROXIMATE THE PROFILE SHAPE
!	!*
!	CALL INTEGRATE(X,DI,NDX,X(1),X(NDX),SUMDI)
!	If(ErrOcc)Then
!		RunSimCode=15
!		return
!	End If
!	CALL INTEGRATE(X,D,NDX,X(1),X(NDX),SUMD)
!	If(ErrOcc)Then
!		RunSimCode=16
!		return
!	End If
!	MassError=(SUMDI-SUMD)/FACTV
!	WRITE(fhoRPT,*)
!	WRITE(fhoRPT,'(" DIFFERENCE IN TOTAL VOLUME BETWEEN FINAL AND INITIAL&
!		& PROFILES:",/,F0.1," ",A)') MassError,Trim(VUNITS)
!	MassError=0.5*(D(1)-DI(1)+D(NDX)-DI(NDX))
!	DO I=2,NDX-1
!		MassError=MassError+D(I)-DI(I)
!	End Do
!	MassError=DX(1)*MassError
!	
!	!*
!	!* CHECK MASS CONSERVATION FOR MEASUREMENTS
!	!*
!	IF(ICOMP==1)THEN
!		XXS=MAX(X(1),XSV(1))
!		XXE=MIN(X(NDX),XSV(NRSV))
!		CALL INTEGRATE(X,DI,NDX,XXS,XXE,SUMDI)
!		If(ErrOcc)Then
!			RunSimCode=15
!			return
!		End If
!		CALL INTEGRATE(XSV,DSV,NRSV,XXS,XXE,SUMD)
!		If(ErrOcc)Then
!			RunSimCode=16
!			return
!		End If
!		WRITE(fhoRPT,*)
!		WRITE(fhoRPT,'(" DIFFERENCE IN TOTAL VOLUME BETWEEN MEASURED AND &
!			&INITIAL PROFILES:",/,F0.1," ",A)') (SUMDI-SUMD)/FACTV,Trim(VUNITS)
!
!		!*
!		!* REPORT SUM OF SQUARES BETWEEN MEASURED AND CALCULATED PROFILES
!		!*
!		WRITE(fhoRPT,*)
!		WRITE(fhoRPT,'(" SUM OF SQUARES OF DIFFERENCES BETWEEN MEASURED AND FINAL&
!			& PROFILES:",/,F0.1," ",A)') SUMQ/FACTL**2.,Trim(L2UNITS)
!	End If
!	!*
!	!* FIND MAXIMUM WATER ELEVATION AND SETUP FOR SIMULATION
!	!*
!	EMXMAX=EMX(ISTMIN)
!	IEMXMAX=ISTMIN
!	DO J=ISTMIN+1,NDX
!		IF(EMX(J)>EMXMAX) THEN
!			EMXMAX=EMX(J)
!			IEMXMAX=J
!		End If
!	End Do
!	WRITE(fhoRPT,*)
!	WRITE(fhoRPT,'( " MAXIMUM VALUE OF WATER ELEVATION + SETUP FOR SIMULATION:"&
!		&,/," ",F0.2," ",A,//," TIME STEP AND POSITION ON"," PROFILE AT WHICH &
!		&MAXIMUM VALUE",/," OF WATER ELEVATION +"," SETUP OCCURRED:",/," ",I0,","&
!		&,F0.1," ",A)')EMXMAX/FACTL,Trim(LUNITS),TEMX(IEMXMAX),X(IEMXMAX)/FACTL,&
!		Trim(LUNITS)
!	write(fhoRPT,*)
!	write(fhoRPT,'( " MAXIMUM ESTIMATED RUNUP ELEVATION:",F0.2," ",A,/," (REFERENCED&
!		& TO VERTICAL DATUM)")') drumax/factl,Trim(LUNITS)
!
!
!	DO J=NDX,1,-1
!		IF(DMN(J)-DI(J)>=EDP(1)) NEDP1=J
!		IF(DMN(J)-DI(J)>=EDP(2)) NEDP2=J
!		IF(DMN(J)-DI(J)>=EDP(3)) NEDP3=J
!	End Do
!	IF(EDP(1)>0.) THEN
!		WRITE(fhoRPT,*)
!		IF(NEDP1==0) THEN
!			WRITE(fhoRPT,FM6000) EDP(1)/FACTL,Trim(LUNITS)
!		ELSE
!			DISTEDP1=XREFELV-X(NEDP1)    
!			WRITE(fhoRPT,FM6005) EDP(1)/FACTL,Trim(LUNITS),X(NEDP1)/FACTL,Trim(LUNITS)
!			WRITE(fhoRPT,*)
!			WRITE(fhoRPT,FM6006) EDP(1)/FACTL,Trim(LUNITS),DISTEDP1/FACTL,Trim(LUNITS)
!		End If
!	End If
!	IF(EDP(2)>0.AND.EDP(2)/=EDP(1)) THEN
!		WRITE(fhoRPT,*)
!		IF(NEDP2==0) THEN
!			WRITE(fhoRPT,FM6000) EDP(2)/FACTL,Trim(LUNITS)
!		ELSE
!			DISTEDP2=XREFELV-X(NEDP2)    
!			WRITE(fhoRPT,FM6005) EDP(2)/FACTL,Trim(LUNITS),X(NEDP2)/FACTL,Trim(LUNITS)
!			WRITE(fhoRPT,*)
!			WRITE(fhoRPT,FM6006) EDP(2)/FACTL,Trim(LUNITS),DISTEDP2/FACTL,Trim(LUNITS)
!		End If
!	End If
!	IF(EDP(3)>0.AND.EDP(3)/=EDP(2).AND.EDP(3)/=EDP(1)) THEN
!		WRITE(fhoRPT,*)
!		IF(NEDP3==0) THEN
!			WRITE(fhoRPT,FM6000) EDP(3)/FACTL,Trim(LUNITS)
!		ELSE
!			DISTEDP3=XREFELV-X(NEDP3)    
!			WRITE(fhoRPT,FM6005) EDP(3)/FACTL,Trim(LUNITS),X(NEDP3)/FACTL,Trim(LUNITS)
!			WRITE(fhoRPT,*)
!			WRITE(fhoRPT,FM6006) EDP(3)/FACTL,Trim(LUNITS),DISTEDP3/FACTL,Trim(LUNITS)
!		End If
!	End If
!	DO J=1,3
!		IF(IELV(J)==0) Cycle
!		WRITE(fhoRPT,*)
!		IF(IELV(J)==1) THEN
!			IF(XELVR(J)<XELVI(J)) THEN
!				WRITE(fhoRPT,'(" MAXIMUM RECESSION OF THE ",F0.2," ",A," ELEVATION &
!					&CONTOUR:",/,F0.1," ",A)') -ELV(J)/FACTL,Trim(LUNITS),&
!					(XELVI(J)-XELVR(J))/FACTL,Trim(LUNITS)
!			ELSE 
!				WRITE(fhoRPT,'(" THE ",F0.2," ",A," ELEVATION CONTOUR DID NOT RECEDE.")')&
!					-ELV(J)/FACTL,Trim(LUNITS)
!			End If
!		End If
!		IF(IELV(J)==-1) THEN 
!			WRITE(fhoRPT,'(" THE ENTIRE PROFILE WAS BELOW THE ",F0.2," ",A," ELEVATION &
!				&CONTOUR",/," AT SOME POINT DURING THE SIMULATION.")') &
!				-ELV(J)/FACTL,Trim(LUNITS)  
!		End If
!	End Do
!	CALL DSHLN(fhoRPT)
!	CALL DSHLN(fhoLOG)
!	WRITE(fhoRPT,*)
!	WRITE(fhoLOG,*)
!
!
!	Return
!End Subroutine PostProcessingResults
!!***************************************************************************************
!Subroutine CompareWithMeasurements(NDX,NRSV,TCOMP,NCOMP,IPP,ICOMP,NRSVM,GOOD,&
!	RunSimCode,R2,FACTL,SUMQ,I)
!	Use Memory, ONLY: XSTART
!	Use MemoryMain, ONLY: D,DSV,X,XSV
!	Use Memoryq_u, ONLY: ErrOcc
!	Use IOHandles, ONLY: fhiPRM
!	Implicit NONE
!	Integer,intent(inout):: NDX,NRSV,TCOMP(10),NCOMP,IPP,ICOMP,NRSVM,I
!	Logical,intent(inout):: GOOD
!	Integer,intent(out):: RunSimCode
!	Real(8),intent(inout):: R2,FACTL,SUMQ
!	
!	Integer:: J,IStat
!	IF(ICOMP==1) THEN
!		DO J=1,NCOMP
!			IF(TCOMP(J)==I) GOOD=.TRUE.
!		End Do
!		IF(GOOD)THEN
!			IPP=IPP+1
!			!*
!			!* READ MEASURED PROFILES TO COMPARE WITH
!			!*
!			CALL HEADER(fhiPRM)
!			If(ErrOcc)Then
!				RunSimCode=12
!				return
!			End IF
!			READ(fhiPRM,*,iostat=istat) NRSV
!			If(istat == 0)Then
!				DO J=1,NRSV
!					READ(fhiPRM,*,iostat=istat) XSV(J),DSV(J)
!					if(istat /= 0)Exit
!					XSV(J)=XSV(J)*FACTL
!					DSV(J)=-DSV(J)*FACTL
!				End Do
!			End If
!			If(istat /= 0)Then
!				Call WriUFLOG('ERROR: INSUFFICIENT NUMBER OF DATA IN MEASURED PROFILE FILE')
!				RunSimCode=13
!				return
!			End if
!			CALL CHKORDER(fhiPRM,XSV,NRSV)
!			If(ErrOcc)Then
!				RunSimCode=14
!				return
!			End If
!			!*
!			!* WRITE MEASURED FINAL PROFILE (ONLY PLOT THOSE POINTS WHICH ARE INSIDE THE GRID)
!			!*
!			IF(IPP==NCOMP)THEN
!				NRSVM=0
!				DO J=1,NRSV
!					IF(XSV(J)>=XSTART.AND.XSV(J)<=X(NDX))NRSVM=NRSVM+1
!				End Do
!			End If
!			!*
!			!* DETERMINE GOODNESS OF FIT FOR THE SIMULATION
!			!*
!			CALL SUMSQ(NDX,X,D,NRSV,XSV,DSV,R2,SUMQ)
!		End If
!	End If
!	Return
!End Subroutine CompareWithMeasurements



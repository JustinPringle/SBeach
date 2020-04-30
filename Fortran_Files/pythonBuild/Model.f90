!***************************************************************************************
Subroutine HNSRATIO(y)
! Initial array for determining Hs of unbroken waves
	Implicit NONE
	Real(8),intent(inout):: y(12:248)
	Integer:: I,J,N
	Real(8):: X1,X2,Y1,Y2,XX,YY

	Real(8):: xt(29),yt(29),x(12:248)
	data (xt(j),j=1,29)/.12,.14,.34,.46,.54,.62,.70,.76,.82,.88,&
		.94,1.0,1.05,1.1,1.16,1.21,1.26,1.32,1.38,1.44,1.5,1.57,1.65,&
		1.74,1.84,2.0,2.16,2.38,2.48/
	data (yt(j),j=1,29)/1.29,1.29,1.295,1.30,1.305,1.31,1.315,&
		1.32,1.325,1.33,1.335,1.34,1.345,1.35,1.355,1.36,1.365,1.37,&
		1.375,1.38,1.385,1.39,1.395,1.4,1.405,1.41,1.413,1.415,1.415/
	x1=xt(1)
	x2=xt(2)
	y1=yt(1)
	y2=yt(2)
	n=3
	do 10 i=12,248
		xx=Real(i)/100.
5		if(xx>=x1.and.xx<=x2) then
			yy=y1+(xx-x1)/(x2-x1)*(y2-y1)
		else
			x1=x2
			y1=y2
			x2=xt(n)
			y2=yt(n)
			n=n+1
			goto 5
		End If
		x(i)=xx
		y(i)=yy
10	continue
	return
end	Subroutine HNSRATIO
!***************************************************************************************
     Subroutine AVHBOT(TIME,IHB,D,NRAV,CRANG,NDX)
! ROUTINE TO DETERMINE AVALANCHING WITH RESPECT TO EXPOSED HARD BOTTOM
!     Use ArrSize
     ! integer, intent(aux):: NDX
     Use Memoryq_u,ONLY: CHECKONLY,ERROCC
     !f2py logical, paramter :: CHECKONLY, ERROCC
     Implicit NONE
     Integer,intent(in)::NDX
     Integer,intent(inout):: NRAV,IHB(0:NDX-1)
     Real(8),intent(inout):: D(0:NDX-1),CRANG,TIME

     Integer:: N1(50),N2(50)
     Integer:: K1,K2,I

	!
	! CHECK FOR EXPOSED HARD BOTTOM
	!
     K1=1
     N1(1)=NRAV
     K2=0
     DO I=NRAV,NDX-1
     IF((IHB(I-1)==1.AND.IHB(I)==0).AND.I/=NRAV)THEN
         K1=K1+1
         N1(K1)=I+1
     End If
     IF(IHB(I-1)==0.AND.IHB(I)==1)THEN
         K2=K2+1
         N2(K2)=I-1
     End If
     End Do
     N2(K2+1)=NDX
	!
	! CHECK FOR AVALANCHING IN EACH SUBAREA SEPARATED BY EXPOSED
	! HARD BOTTOM
	!
!write(*,*) NRAV,N1(K1),N2(K2)
!     DO I=0,K2-1
!!write(*,*) '74'
!         CALL MAXANG(TIME,D,N2(I),CRANG,NDX)
!!         CALL MAXANG(TIME,D,N1(I),CRANG,NDX) !,N2(I)
!         if(ErrOcc) return
!     End Do
!     DO I=K2,K1-1
!!write(*,*) '74'
!!         CALL MAXANG(TIME,D,N2(I),CRANG,NDX)
!         CALL MAXANG(TIME,D,N1(I),CRANG,NDX) !,N2(I)
!         if(ErrOcc) return
!     End Do
     DO I=0,K1-1
!write(*,*) '74'
!         CALL MAXANG(TIME,D,N2(I),CRANG,NDX)
         CALL MAXANG(TIME,D,N1(I),CRANG,NDX) !,N2(I)
         if(ErrOcc) return
     End Do
     RETURN
     End Subroutine AVHBOT

!***************************************************************************************
!     Subroutine HRDBOT(NRU,DP,D,QP,Q,IHB,DBOT,SCF)
!! THIS ROUTINE COMPUTES DEPTH CHANGES IN THE PRESENCE OF
!! HARD BOTTOMS
!    	Use ArrSize, ONLY: NDX
!     !f2py integer,intent(aux)::NDX
!    	Use Memory, ONLY: DT,XSTART !DX,
!     !f2py real(8), intent(aux):: DX,DT,XSTART
!    	Implicit NONE
!      Real(8),intent(in),dimension(0:NDX-1)::DX
!    	Integer,intent(inout):: IHB(0:NDX-1)
!    	Real(8),intent(inout):: DP(0:NDX-1),D(0:NDX-1),QP(0:NDX),Q(0:NDX),DBOT(0:NDX-1),SCF
!    	Integer,intent(inout):: NRU
!    
!    	Integer::ICROSS(50),NCROSS,IDIR(50),IS,ISGN,II,IL1,IL2,ILEX
!    	Integer:: K,NStart,NStop,NS,I,IHBB,IHBBT,IEXP
!    
!    	Real(8),allocatable:: QBOT(:),QAKT(:)
!    	Real(8):: DTemp,DIR,SUMQ
!    
!    	Allocate(QBOT(0:NDX),QAKT(0:NDX))
!    	!
!    	! POTENTIAL TRANSPORT (IGNORING EXPOSED HARD BOTTOMS)
!    	!
!    	DO I=NRU,NDX
!    		QAKT(I)=0.5*(QP(I)+Q(I))
!    		QBOT(I)=QAKT(I)
!    	End Do
!    	!
!    	! LOCATE OFFSHORE (+) AND ONSHORE (-) TRANSPORT AREAS
!    	!
!    	NCROSS=1
!    	ICROSS(1)=NRU
!    	!
!    	! FIND TRANSPORT DIRECTION IN FIRST AREA
!    	!
!    	DO I=NRU+1,NDX-1
!    		IF(QAKT(I)/=0)THEN
!    			ISGN=NINT(ABS(QAKT(I))/QAKT(I))
!    			Exit
!    		End If
!    	End Do
!    	IDIR(1)=ISGN
!    	DO I=NRU+1,NDX-1
!    		IF(QAKT(I)>0) IS=1
!    		IF(QAKT(I)<0) IS=-1
!    		IF(QAKT(I)==0) IS=ISGN
!    		IF(IS/=ISGN)THEN
!    			ISGN=IS
!    			NCROSS=NCROSS+1
!    			ICROSS(NCROSS)=I-1
!    			IDIR(NCROSS)=ISGN
!    		End If
!    	End Do
!    	ICROSS(NCROSS+1)=NDX-1
!    	!
!    	! APPLY SAND CONSERVATION EQUATION AND CHECK FOR EXPOSED HARD
!    	! BOTTOMS. CORRECT TRANSPORT RATES IF NECESSARY; DO THE
!    	! CORRECTIONS INDEPENDENTLY WITHIN EACH TRANSPORT AREA
!    	!
!    	IHBBT=0
!    	DO K=1,NCROSS
!    		IF(IDIR(K)==1)THEN
!    			NSTART=ICROSS(K)
!    			NSTOP=ICROSS(K+1)-1
!    			NS=1
!    			IF(K==NCROSS) NSTOP=NSTOP+1
!    		ELSE
!    			NSTART=ICROSS(K+1)+1
!    			NSTOP=ICROSS(K)+2
!    			NS=-1
!    			IF(K==1) NSTOP=NSTOP-1
!    		End If
!    		DIR=Real(IDIR(K))
!    		IHBB=0
!    		DO I=NSTART,NSTOP,NS
!    			IF(NS==1)THEN
!    				II=I
!    			ELSE
!    				II=I-1
!    			End If
!    			IHB(II)=0
!    			D(II)=DP(II)+DT/DX(II)*(QAKT(I+NS)-QAKT(I))*DIR
!    			IF(D(II)>DBOT(II))THEN
!    				IHBB=1
!    				IHBBT=1
!    			End If
!    			IF(D(II)>=DBOT(II)) IHB(II)=1
!    		End Do
!    		!
!    		! EXPOSED HARD BOTTOM ENCOUNTERED WITHIN THE REGION
!    		!
!    		IF(IHBB==1)THEN
!    			IEXP=0
!    			SUMQ=0
!    			DO I=NSTART,NSTOP,NS
!    				IF(NS==1)THEN
!    					II=I
!    					IL1=I-1
!    					IL2=I
!    				ELSE
!    					II=I-1
!    					IL1=I
!    					IL2=I-1
!    				End If
!    				IF(D(II)<DBOT(II))THEN
!    					IF(IEXP==0)THEN
!    						QBOT(I+NS)=QAKT(I+NS)
!    					ELSE
!    						SUMQ=SUMQ+DX(II)
!    						QBOT(I+NS)=QAKT(I+NS)+(QBOT(ILEX)-QAKT(I+NS))*EXP(-SCF*SUMQ)
!    						IF(ABS(QBOT(I+NS))>ABS(QAKT(I+NS)))QBOT(I+NS)=QAKT(I+NS)
!    					End If
!    					!
!    					! CHECK IF ADDITIONAL HARD BOTTOMS ARE EXPOSED AFTER MODIFYING
!    					! TRANSPORT RATES
!    					!
!    					IF(I>1)THEN
!    						IF(IHB(IL1)==1.AND.IHB(IL2)==0)THEN
!    							IEXP=1
!    							ILEX=IL1
!    							SUMQ=SUMQ+DX(II)
!    							QBOT(I+NS)=QAKT(I+NS)+(QBOT(ILEX)-QAKT(I+NS))*EXP(-SCF*SUMQ)
!    							IF(ABS(QBOT(I+NS))>ABS(QAKT(I+NS)))QBOT(I+NS)=QAKT(I+NS)
!    							DTEMP=DP(II)+DT/DX(II)*(QBOT(I+NS)-QBOT(I))*DIR
!    							IF(DTEMP>=DBOT(II))THEN
!    								QBOT(I+NS)=QBOT(I)+DX(II)/DT*(DBOT(II)-DP(II))*DIR
!    								IHB(II)=1
!    								IEXP=0
!    								SUMQ=0.0
!    							End If
!    						End If
!    					End If
!    				ELSE
!    					!
!    					! CHECK IF TRANSPORT CORRECTIONS START IN THE MIDDLE OF THE GRID
!    					!
!    					IF(I==NSTART)THEN
!    						IF(I==NRU.OR.I==NDX+1)THEN
!    							QBOT(I+NS)=QBOT(I)+DX(II)/DT*(DBOT(II)-DP(II))*DIR
!    						ELSE
!    							!
!    							! FOR STARTING CELL IN THE MIDDLE OF THE GRID, DIVIDE TRANSPORTED
!    							! AMOUNT OUT FROM THE CELL PROPORTIONALLY
!    							!
!    							QBOT(I+NS)=QAKT(I+NS)*(DBOT(II)-DP(II))/(D(II)-DP(II))
!    						End If
!    					ELSE
!    						QBOT(I+NS)=QBOT(I)+DX(II)/DT*(DBOT(II)-DP(II))*DIR
!    					End If
!    				End If
!    			End Do
!    		End If
!    	End Do
!    	IF(IHBBT==0) RETURN
!    	!
!    	! NEW DEPTHS
!    	!
!    	DO I=NRU,NDX-1
!    		D(I)=DP(I)+DT/DX(I)*(QBOT(I+1)-QBOT(I))
!    	End Do
!    	DeAllocate(QBOT,QAKT)
!    
!    	RETURN
!     End Subroutine HRDBOT

!***************************************************************************************
    Subroutine InitialProfile(TPIN,DTI,XT,DX,NDX,NRV)!,DTI,XT,NRV,NDX,DX,XSTART)!XLAND,DLAND,XLBDUNE,DLBDUNE,XLCDUNE,DLCDUNE,XSCDUNE,&
!    	DSCDUNE,XBERMS,DBERMS,XBERME,DBERME,XFORS,DFORS,FACTL,DX,XSTART)
    	Use Memory, ONLY: DT,XSTART
!     use MemoryConst
     !f2py real(8), intent(aux)::DT,XSTART
     Use MemoryMain, ONLY: DI,X
    !f2py real(8), intent(aux):: DI,X
     
     ! real(8), intent(aux),dimension(0:NDX-1)::DX
!     DX(0:NDX-1),
!    	Use Memoryq_u, ONLY: CHECKONLY,ERROCC,PriX,PriD,PriNRV
!    	Use IOHandles, ONLY: fhiPRI
    	Implicit NONE
    !* THIS ROUTINE GENERATES THE INITIAL PROFILE. AN ARBITRARY OR A SCHEMATIZED PROFILE
    !* SHAPE MAY BE CHOSEN.
    	Integer:: NDX,NRV
    	Integer::TPIN
      Real(8),intent(in),dimension(0:NDX-1)::DX
!    	Real(8),intent(out),dimension(0:NDX-1):: DI,X!,DX
!      real(8),intent(inout):: DX(NDX)!,XLAND,DLAND,XLBDUNE,DLBDUNE,XLCDUNE,DLCDUNE,&
!    		XSCDUNE,DSCDUNE,XBERMS,DBERMS,XBERME,FACTL,DBERME,XFORS,DFORS, DX(NDX)
    	Real(8):: SUMX,CDEAN,XDEAN!, XSTART
    	Integer:: ISTAT,I,J
     Real(8), intent(inout):: DTI(NRV),XT(NRV)

    	!*
    	!* DETERMINE THE COEFFICIENT IN DEAN'S EQUILIBRIUM EQUATION AND
    	!* THE CORRESPONDING ENERGY DISSIPATION
    	!*
    	CALL EQCOF(CDEAN)
!    	!*
!    	!* GENERATE X-VALUES
!    	!*
    	SUMX=XSTART
    	DO I=0,NDX-1
    		X(I)=SUMX
!           if(I==0) write(*,*) DX(I)
    		IF(I<NDX-1) SUMX=SUMX+(DX(I)+DX(I+1))/2
    	End Do
    	!*
    	!* READ AN ARBITRARY OR A SCHEMATIZED PROFILE
    	!*
    	IF(TPIN==1)THEN 
    		!*
    		!* ARBITRARY PROFILE
    		!*
!    		READ(fhiPRI,*,iostat=istat) NRV
!    		If(istat == 0)Then
!    			DO I=1,NRV
!    				READ(fhiPRI,*,iostat=istat) XT(I),DTI(I)
!    				if(istat /= 0)Exit
!    				XT(I)=XT(I)*FACTL
!    				DTI(I)=-DTI(I)*FACTL
!    			End Do
!    			If(istat == 0)CALL CHKORDER(fhiPRI,XT,NRV)
!    			If(ErrOcc) Return
!    		End If
!    		If(istat/=0)Then
!    			Call WriUFLOG('ERROR: INSUFFICIENT NUMBER OF DATA IN INITIAL PROFILE FILE.')
!    			ErrOcc=.true.
!    			Return
!    		End If
!    		IF(XT(1)>XSTART.OR.XT(NRV)<X(NDX)) THEN
!    			Call WriUFLOG('ERROR: INITIAL PROFILE DATA DOES NOT COVER RANGE OF SIMULATION')
!    			ErrOcc=.true.
!    			Return
!    		End If
    		DO I=0,NDX-1
    			DO J=0,NRV-2
    				IF(X(I)>=XT(J).AND.X(I)<XT(J+1)) Exit
    			End Do
    			DI(I)=DTI(J)+(DTI(J+1)-DTI(J))*(X(I)-XT(J))/(XT(J+1)-XT(J))
    		End Do
    	!*
    	!* SCHEMATIZED PROFILE
    	!*
!    	ELSE
!    		XDEAN=(DFORS/CDEAN)**1.5
!    		DO I=1,NDX
!    			IF(X(I)>=XLAND.AND.X(I)<=XLBDUNE)THEN
!    				DI(I)=DLAND+(DLBDUNE-DLAND)*(X(I)-XLAND)/(XLBDUNE-XLAND)
!    			ELSEIF(X(I)>XLBDUNE.AND.X(I)<=XLCDUNE)THEN
!    				DI(I)=DLBDUNE+(DLCDUNE-DLBDUNE)*(X(I)-XLBDUNE)/(XLCDUNE-XLBDUNE)
!    			ELSEIF(X(I)>XLCDUNE.AND.X(I)<=XSCDUNE)THEN
!    				DI(I)=DLCDUNE+(DSCDUNE-DLCDUNE)*(X(I)-XLCDUNE)/(XSCDUNE-XLCDUNE)
!    			ELSEIF(X(I)>XSCDUNE.AND.X(I)<=XBERMS)THEN
!    				DI(I)=DSCDUNE+(DBERMS-DSCDUNE)*(X(I)-XSCDUNE)/(XBERMS-XSCDUNE)
!    			ELSEIF(X(I)>XBERMS.AND.X(I)<=XBERME)THEN
!    				DI(I)=DBERMS+(DBERME-DBERMS)*(X(I)-XBERMS)/(XBERME-XBERMS)
!    			ELSEIF(X(I)>XBERME.AND.X(I)<=XFORS)THEN
!    				DI(I)=DBERME+(DFORS-DBERME)*(X(I)-XBERME)/(XFORS-XBERME)
!    			ELSE
!    				DI(I)=CDEAN*(XDEAN+X(I)-XFORS)**0.667
!    			End If
!    		End Do
    	End If
    	RETURN
    End Subroutine InitialProfile

    !***************************************************************************************
     Subroutine VARGRID(DXV,NDXV,DX,NDX)
    !* THIS ROUTINE GENERATES A GRID WITH A VARIABLE STEP LENGTH.
    	Use MemoryMain, ONLY: NrNDXV,NrDXV
     !f2py integer, intent(aux)::NrNDXV,NrDXV
    	Implicit NONE
    	Integer,intent(in):: NDX
    	Integer,intent(inout),dimension(0:NrNDXV-1):: NDXV
    	Real(8),intent(inout),dimension(0:NrDXV-1):: DXV
      Real(8), intent(out),dimension(0:NDX-1):: DX
    	Integer::I,J,NSUM
    
    	!*
    	!* GENERATE GRID
    	!*
      !write(*,*) NrNDXV
    	J=0
    	NSUM=NDXV(0)
!      write(*,*) NSUM
    	DO I=0,NDX-1
    		DX(I)=DXV(J)
!           if(I==0) write(*,*) DX(I)
    		IF(I==NSUM-1.AND.I<NDX-1) THEN
    			J=J+1
    			NSUM=NSUM+NDXV(J)
    		End If
    	End Do
    	RETURN
     End Subroutine VARGRID

!*******************************************************************************
     Subroutine PreProcessWaves(IWAVE,IRAND,OFFSHORE,HBEG,T,DMEAS,RPERC,&
    	HO,HIN,DNDX_Value,ZIN,ZBEG)
    	Use Memory, ONLY: CC,DSURGE,LO
     !f2py real(8), intent(aux)::CC,DSURGE,LO
    	Implicit NONE
    	Integer,intent(inout):: IWAVE,IRAND
    	Logical,intent(inout):: OFFSHORE
     real(8),intent(out):: HBEG,HO,ZBEG
    	Real(8),intent(inout):: T,DMEAS,RPERC,HIN,DNDX_Value,ZIN
     real(8):: ARG,DTOT,FTEMP,CGT,CNT,L1,LT,CG1,CN1,HOUT,hinc,ZO
    	ZO=ZIN
    	IF(IRAND==1)THEN
    		IF(IWAVE==0) HIN=HINC
    		CALL RANDOM_I(RPERC,HIN,HOUT)
    		HIN=HOUT
!           write(*,*) HOUT
    	End If
    	IF(OFFSHORE) HIN=0.0
    	IF(HIN==0.0) THEN
    		HBEG=0.0
    		HO=0.0
    		ZO=0.0
    	Else
		!*
		!* DETERMINE DEEPWATER WAVE HEIGHT AND INPUT WAVE HEIGHT AT END OF GRID
		!*
    		IF(abs(DMEAS)<=0.0001)THEN
    			LT=LO
    			FTEMP=CC*HIN**2*LT/T/2
    			HO=HIN
    		ELSE
    			DTOT=DMEAS+DSURGE
    			CALL LINWAV(T,DTOT,LT,CGT,CNT)
    			FTEMP=CC*HIN**2*CGT
    		End If
		!*
		!* WAVE HEIGHT AN ANGLE AT THE BEGINNING OF THE GRID
		!*
    		CALL LINWAV(T,DNDX_Value,L1,CG1,CN1)
    		ARG=ABS(L1/LT*SIN(ZIN))
    		IF (ARG>1.) THEN
    			HBEG=HIN
    			ZBEG=ZIN
    		ELSE
    			ZBEG=ASIN(L1/LT*SIN(ZIN))
    			HBEG=SQRT(FTEMP*COS(ZIN)/COS(ZBEG)/CC/CG1)
    		End If
		!*
		!* DEEPWATER WAVE HEIGHT
		!*
    		IF(abs(DMEAS)>0.0)THEN
    			ARG=ABS(LO/LT*SIN(ZIN))
    			IF(ARG>1.) THEN
    				HO=HIN
    				ZO=ZIN
    			ELSE
    				ZO=ASIN(LO/LT*SIN(ZIN))
    				HO=SQRT(FTEMP*COS(ZIN)/COS(ZO)/CC/LO*T*2)
    			End If
    		End If
    	End If
    	return
     End Subroutine PreProcessWaves

!***************************************************************************************

!***************************************************************************************
     Subroutine MAXANG(TIME,D,NFS,bmax,NDX)
!* THIS ROUTINE LIMITS THE STEEPNESS OF THE PROFILE BY ALLOWING FOR AVALANCHING IF THE
!* SLOPE LOCALLY EXCEEDS A PREDEFINED MAXIMUM VALUE. THE AVALANCHING IS CONSIDERED TO
!* OCCUR INFINITELY FAST IN COMPARISON WITH THE TIME STEP OF THE MODEL. SAND IS REDISTRIBUTED
!* IN THE NEIGHBORING CELLS AT A SLOPE CORRESPONDING TO THE STABLE SLOPE AFTER AVALANCHING.
!     use ArrSize
     ! integer, intent(aux) :: NDX
     Use Memory,ONLY: DX,DT,XSTART
     !f2py real(8), intent(aux) :: DX,DT,XSTART
     Use Memoryq_u,ONLY: CHECKONLY,ERROCC
     Implicit NONE
     integer, intent(in) :: NDX
     Integer,intent(inout):: NFS
     Integer:: I,NCORR,NDIR,NVAL,J,L,IMAX,NTIME
     Real(8),intent(inout):: D(0:NDX-1),BMAX,TIME
     Real(8):: DSUM,CSUM1,CSUM2,BAV,BPMAX,DCORR,&
     DCR,DDCRIT
     Real(8),allocatable:: DD(:)
     LOGICAL CHECK

     Real(8):: BTEMP
     Integer::II

     Allocate(DD(0:NDX-1))
!write(*,*) NFS,bmax
     NTIME=NINT(TIME/DT)
     L=0
     bav=bmax-0.175
	!*
	!* DETERMINE THE MAXIMUM SLOPE ALONG THE PROFILE
	!*
88   CHECK=.FALSE.
     L=L+1
     BPMAX=0
!     write(*,*) nfs,errocc
     DO I=NFS,NDX-1
         BTEMP=2*(D(I)-D(I-1))/(DX(I)+DX(I-1))
         IF(ABS(BTEMP)>ABS(BPMAX))THEN
             BPMAX=BTEMP
             IMAX=I
         End If
     End Do
	!*
	!* CHECK IF THE MAXIMUM SLOPE EXCEEDS THE ANGLE OF REPOSE
	!*
     IF(ATAN(ABS(BPMAX))>BMAX) CHECK=.TRUE.
     IF(CHECK)THEN
!write(*,*) IMAX
		!*
		!* REDISTRIBUTE THE SAND IN THE NEIGHBORING CELLS ACCORDING TO THE
		!* SLOPE AFTER AVALANCHING HAS OCCURRED (BAV)
		!*
         IF(BPMAX>0)THEN
             NDIR=-1
         ELSE
         NDIR=1
         IMAX=IMAX-1
         End If
         NCORR=IMAX
77       NCORR=NCORR+NDIR
		!*
		!* CONTINUE TO REDISTRIBUTE SAND UNTIL THE FIRST CELL OUTSIDE THE AREA
		!* WHERE SAND IS MOVED HAS A SLOPE LOWER THAN BAV
		!*
         NVAL=1
         J=1
         DSUM=0.
         CSUM1=0.0
         CSUM2=0.0
		!         IF(NTIME<=300) WRITE(21,*) IMAX+NDIR,NCORR,NDIR
!write(*,*) IMAX+NDIR,NCORR,NDIR
         DO 20 I=IMAX+NDIR,NCORR,NDIR
             J=J+1
             IF(I>NDX-1.OR.I<NFS-1)THEN
                if(i>NDX-1) then
                 write(*,*)'"Avalanching beyond seaward boundary at time step:"&
                         &I0" --extend boundary"',NTIME 
!        Call WriFILOG('"AVALANCHING BEYOND SEAWARD BOUNDARY AT TIME STEP: "&
!						&I0" --EXTEND BOUNDARY."',NTIME)
!                else
!        Call WriFILOG('"AVALANCHING BEYOND LANDWARD BOUNDARY AT TIME STEP: "&
!        &I0" --EXTEND BOUNDARY."',NTIME)
                End If
                ErrOcc=.true.
                    return
                End If
                DSUM=DSUM+D(I)*DX(I)
                DCR=(DX(I)+DX(I-NDIR))/2*TAN(BAV)
                DO II=I,NCORR,NDIR
                    CSUM1=CSUM1+DCR*DX(II)
                End Do
                CSUM2=CSUM2+DX(I)
                NVAL=NVAL+1
20          CONTINUE
         DCORR=(-D(IMAX)*CSUM2+DSUM+CSUM1)/(CSUM2+DX(IMAX))
         DD(IMAX)=D(IMAX)+DCORR
!		         IF(NTIME<=300) WRITE(*,*) IMAX,D(IMAX),DD(IMAX),DCORR
		!*
		!* CALCULATE THE NEW DEPTHS
		!*
         J=1
         CSUM1=0.0
            DO I=IMAX+NDIR,NCORR,NDIR
                J=J+1
                DCR=(DX(I)+DX(I-NDIR))/2*TAN(BAV)
                CSUM1=CSUM1+DCR
                DCORR=DD(IMAX)-D(I)-CSUM1
                DD(I)=D(I)+DCORR
            End Do
         IF(NCORR+NDIR>=NFS-1)THEN
             IF(NDIR==-1)THEN
                DDCRIT=2*(DD(NCORR)-D(NCORR+NDIR))/(DX(NCORR)+DX(NCORR+NDIR))
                    IF(DDCRIT>TAN(BAV)) GOTO 77
                    ELSE
                    DDCRIT=2*(D(NCORR+NDIR)-DD(NCORR))/(DX(NCORR+NDIR)+DX(NCORR))
                IF(DDCRIT<-TAN(BAV)) GOTO 77
                End If
            End If
            DO I=IMAX,NCORR,NDIR
                D(I)=DD(I)
            End Do
        End If
	!*
	!* GO BACK AND CHECK IF THE ANGLE OF REPOSE IS EXCEEDED IN OTHER
	!* PARTS OF THE PROFILE
	!*
     IF(CHECK.AND.L<100) GOTO 88
     IF(L>=100)THEN
     write(*,*)'"AVALANCHING ROUTINE NOT CONVERGING AT TIME STEP: "I0&
     &" -- SIMULATION TERMINATED."',NTIME
     ErrOcc=.true.
     return
     End If
     DeAllocate(DD)
     RETURN
     End Subroutine MAXANG

!***************************************************************************************
!!
!! HERE ARE THE SUBROUTINES FOR RANDOM WAVES
!! -----------------------------------------
!!
!!                     |
!!                     |
!!                    \ /
!!                     !
!!
!!
!***************************************************************************************
     Subroutine WAVRAN(DX,NSWALL,HBEG,T,ZBEG,D,H,E,ISTOP,ZW,DISS,BROK,X,ISWLIM,hnsr,w,zwind,da,NDX)
     !* THIS ROUTINE DETERMINES THE TRANSFORMATION OF THE RMS WAVE HEIGHT (BROKEN AND
     !* NON-BROKEN WAVES) ACROSS AN ARBITRARY BEACH PROFILE. DALLY'S BREAKER DECAY MODEL
     !* IS GENERALIZED TO Include IRREGULAR WAVES.
     !* FOR REFERENCE LOOK AT:
     !*
     !* LARSON, M. "MODEL FOR DECAY OF RANDOM WAVES IN THE SURF ZONE," JOURNAL OF WATERWAY,
     !* PORT, COASTAL, AND OCEAN ENGINEERING, VOL 121, NO. 1, PP 1-12.
     !*
     !* IN
     !* --
     !* NDX = NUMBER OF CELLS IN THE GRID
     !* DX = GRID SIZE
     !* HBEG = RMS WAVE HEIGHT AT THE OFFSHORE END OF THE GRID
     !* T = WAVE PERIOD
     !* ZBEG = WAVE ANGLE AT THE OFFSHORE END OF THE GRID
     !* D = PROFILE ELEVATINOS (POSITIVE BELOW STILL WATER LEVEL)
     !*
     !* OUT
     !* ---
     !* H = RMS WAVE HEIGHT, VARIATION WITH x
     !* E = MEAN WAVE SETUP, VARIATION WITH x
     !* ISTOP = MOST SHOREWARD CELL WHERE THE CALCULATION STOPPED
     !* ZW = WAVE ANGLE, VARIATION WITH x
     !* DISS = WAVE BREAKING ENERGY DISSIPATION, VARIATION WITH x
     !* BROK = MEAN RATIO OF BROKEN WAVES, VARIATION WITH x
     
!    	Use ArrSize, ONLY: NDX
!     ! integer, intent(aux):: NDX
    	Use Memory,ONLY: CC,PI,LO
     !f2py real(8), intent(aux):: CC,PI,LO
    	Use MemoryConst,ONLY: KAPPA,GAMSTB,KREF,IRMAX
     !f2py real(8), intent(aux):: KAPPA,GAMSTB,KREF,IRMAX
    	Implicit NONE

    	Integer,intent(inout):: ISTOP,NSWALL,ISWLIM, NDX
      real(8), intent(inout), dimension(0:NDX)::H,ZW,E
     
    	Real(8),intent(inout):: D(0:NDX-1),&
    			HBEG,ZBEG,T,DISS(0:NDX-1),BROK(0:NDX-1),W,ZWIND,&
    			DX(0:NDX-1),X(0:NDX-1),hnsr(12:248),da(0:NDX-1)

    	Real(8):: L1,CG1,CN1,GAMBR,L2,CG2,CN2,REF,HR2N,HR2O,AC,FSTAB,B,BMIN,REFO,&
    			REFN,DMIN,DAKT,CGMIN,HRMSO,HRMSM,SUMPDF,&
    			SUMHR2,HREF,HMAX,HR2,HN2,DIFFP,CTRAN1,CTRAN2,COSMIN, error

    	Real(8),allocatable:: F(:),SXX(:),DH(:),PDF(:),HRPDF(:)
    	LOGICAL START,fndlim
    
    	Integer:: J,JTRY,JINC,JFIX,JJFIX,I,ILOOP,IPDF,IPWAVE,IPMAX
    	Real(8)::CCH2,CGBEG,ETRY,EMEAN,AO,FO,FSO,HH1,HH2,HBRATIO,HNSFACT,HNS2
    	Real(8):: CWIND,FORC
    
    	FSTAB=0.0
    
    	Allocate(F(0:NDX),SXX(0:NDX),DH(0:NDX),PDF(0:IRMAX-1),HRPDF(0:IRMAX-1))
    	GAMBR=0.78
! write(*,*) H
     do i = 0,NDX-1
        	H(I)=0.
           ZW(I) =0.
           E(I) =0.
           BROK(I)=0.
           DISS(I)=0.
           F(I) =0.
           SXX(I)=0.
           DH(I)=0.
           
     end do
    	!*
    	!*
    	!* INITIALIZE ARRAY FOR REFORMED WAVES
    	!*
    	!* IRMAX is the number of 
    
    	DO I=0,IRMAX-1
    		HRPDF(I)=4*HBEG*Real(I+1)/Real(IRMAX)
!           write(*,*) HRPDF(I)
    		PDF(I)=0.0
    	End Do

    	!*
    	!*
    	!*  DETERMINE AVERAGE DEPTH BETWEEN NEIGHBORING CELLS
    	!*
    	!*
!write(*,*) '512',nswall 
    	DH(nswall)=D(nswall)
    	DO I=nswall+1,NDX-1
!write(*,*) I
    		DH(I)=D(I-1)+(D(I)-D(I-1))*DX(I-1)/(DX(I)+DX(I-1))
    	End Do
    	DH(NDX)=1.01*D(NDX-1)
      
    	!*
    	!*
    	!*  INITIALIZE VARIABLES
    	!*
    	!*

    	START=.TRUE.
    	ISWLIM=ndx-1
    	fndlim=.false.
    	BMIN=1.0
    	REFO=0.0
    	HR2O=0.0
    	IPMAX=-1
    	DMIN=1E6

    	!*
    	!*
    	!*  DETERMINE WAVE PROPERTIES IN THE MOST SEAWARD CELL
    	!*
    	!*
!write(*,*) HBEG
    	H(NDX)=HBEG*1.416
!write(*,*) HBEG
    	ZW(NDX)=ZBEG
    	CCH2=CC*HBEG**2*1.416**2.
!write(*,*) '543'
    	CALL LINWAV(T,DH(NDX),L1,CG1,CN1)
!write(*,*) '544'
    	CGBEG=CG1
    	E(NDX)=-PI*1.416**2.*HBEG**2/4./L1/SINH(4*PI*DH(NDX)/L1)
    	F(NDX)=CCH2*CG1
    	SXX(NDX)=CCH2*(CN1*((COS(ZBEG))**2+1)-0.5)


    	ISTOP=0 !THIS IS THE VALUE WHEN THERE IS NO OTHER STOPPING POINT FOUND DURING CALCULATIONS

    	!*
    	!* PROCEED FROM THE END OF THE GRID IN THE SHOREWARD DIRECTION
    	!* WHEN CALCULATING CROSS-SHORE WAVE PROPERTIES
    	!*
    	!* MAIN LOOP IN SPACE
    	!* ------------------
    	!*

    	MSpace: DO I=NDX-1,0,-1
            
		!*
		!* DETERMINE WAVE PROPERTIES USING LINEAR WAVE THEORY
		!*
		!* EXTRAPOLATE SETUP CLOSE TO THE SHORELINE
		!*
!           write(*,*) i,H(i+1)
    		ETRY=E(I+1)
    		IF(DH(I)<0) ETRY=E(I+1)+DX(I)/DX(I+1)*(E(I+1)-E(I+2))
    		IF(DH(I)+ETRY<=0)THEN
    			ISTOP=I+1
!                 write(*,*) I,'Stopping at 561'
    			Exit MSpace !GOTO 77
    		End If
    
    		IF(I<NSWALL)THEN
    			ISTOP=I+1
!                 write(*,*) I,'Stopping at 567'
    			Exit MSpace !GOTO 77
    		End If
    
    		CALL LINWAV(T,DH(I)+ETRY,L2,CG2,CN2)
    
    		!*
    		!*
    		!* CALCULATE WAVE REFRACTION FROM SNELLS LAW
    		!*
    		!*
    		ZW(I)=ASIN(L2/L1*SIN(ZW(I+1)))
    
    		!*
    		!*
    		!* DEFINE CONSTANTS
    		!*
    		!*
    		CTRAN1=CG2*COS(ZW(I))!0.25*(CG1+CG2)*(COS(ZW(I))+COS(ZW(I+1)))!(CG2*COS(ZW(I)))!0.25*(CG1+CG2)*(COS(ZW(I))+COS(ZW(I+1)))
    		CTRAN2=SQRT(CG1*COS(ZW(I+1))/CG2/COS(ZW(I)))
    
    		!*
    		!* APPLY DALLY'S BREAKER DECAY MODEL
    		!* ---------------------------------
    		!*
    		!*
    		EMEAN=0.5*(ETRY+E(I+1)) !average setup

    		!*
    		!*
    		!* RMS HEIGHT NEGLECTING WAVE BREAKING, Eq 27
    		!*
    		!*
!write(*,*) CGBEG
    		HRMSO=HBEG*SQRT(CGBEG*COS(ZBEG)/CTRAN1)!H(I+1)*SQRT(CG1*COS(ZW(i+1))/CTRAN1)!HBEG*CGBEG*COS(ZBEG)/CTRAN1)!CG1*COS(ZW(I)))!CTRAN1)
    
    		!*
    		!*
    		!* PERCENT OF UNBROKEN WAVES, beta=1-alpha, alpha by Eq 28
    		!*
    		!*
    		B=1.-EXP(-(GAMBR*(D(I)+EMEAN)/HRMSO)**2)

    		!*
    		!*
    		!* NEVER LET B INCREASE, BOOK KEEPING
    		!*
    		!*
!           write(*,*) i,B
    		IF(B<BMIN)THEN
    			BMIN=B
    			DMIN=D(I)+EMEAN
    			HRMSM=HRMSO
    			CGMIN=0.5*(CG1+CG2)
    			COSMIN=0.5*(COS(ZW(I))+COS(ZW(I+1)))
    		End If

    		!*
    		!* CHECK IF REFORMATION MAY OCCUR
    		!* ------------------------------
    		!*
    		!*
    		IF(DH(I)+ETRY>DH(I+1)+E(I+1))THEN !depth increase in the shoreward direction
    			IF(START)THEN
    				START=.FALSE.
    				IF(I+2<=NDX+1)THEN
    					FO=F(I+2)
    				ELSE
    					FO=F(I+1)
    				End If
    				FSO=FSTAB
    				AO=1-BMIN-REFO
    			End If
    
    			!*
    			!*
    			!* PREDICT PERCENTAGE OF WAVES REFORMED
    			!*
    			!*
    			IF(F(I+1)>=1.001*FSTAB.and.fo>=1.001*fso)THEN
    				REFN=1-BMIN-AO*((F(I+1)-FSTAB)/(FO-FSO))**KREF
    			ELSE
    				REFN=1-BMIN
    			End If

    			!*
    			!*
    			!* NEVER LET PERCENTAGE OF REFORMED WAVES DECREASE
    			!*
    			!*
    			IF(REFN-REFO<0) REFN=REFO
    			DIFFP=REFN-REFO
    
    			!*
    			!*
    			!* REFORMED WAVE HEIGHT
    			!*
    			!*
    			HREF=GAMSTB*(D(I)+EMEAN)
    
    			!*
    			!*
    			!* SHOAL DISTRIBUTION OF REFORMED WAVES
    			!*
    			!*
    			IF(IPMAX>0)THEN
    				DO J=0,IRMAX-1
    					HRPDF(J)=HRPDF(J)*CTRAN2
    				End Do
    			End If

    			!*
    			!*
    			!* STORE PDF FOR REFORMED WAVES
    			!*
    			!*

    			DO J=1,IRMAX-2
    				HH1=0.5*(HRPDF(J-1)+HRPDF(J))
    				HH2=0.5*(HRPDF(J)+HRPDF(J+1))
    				IF(HREF>=HH1.AND.HREF<HH2)THEN
    					IPDF=J
    					EXIT
    				End If
    				IF(J==1.AND.HREF<HH1)THEN
    					IPDF=1
    					EXIT
    				End If
    				IF(J==IRMAX-2.AND.HREF>=HH2)THEN
    					IPDF=IRMAX
    					EXIT
    				End If
    			End Do
    			IF(IPDF>IPMAX) IPMAX=IPDF
    			PDF(IPDF)=PDF(IPDF)+DIFFP

    			!*
    			!*
    			!* RMS HEIGHT FOR UNBROKEN WAVES
    			!*
    			!*
    			HN2=(HRMSO**2-EXP(-(GAMBR*DMIN/HRMSM)**2)*(HRMSO**2+&
    				(GAMBR*DMIN)**2*CGMIN*COSMIN/CTRAN1))/BMIN
    
    			!*
    			!*
    			!* RMS HEIGHT FOR REFORMED WAVES
    			!*
    			!*
    			SUMPDF=0.0
    			SUMHR2=0.0
    			DO J=0,IPMAX-1
    				SUMPDF=SUMPDF+PDF(J)
    				SUMHR2=SUMHR2+PDF(J)*HRPDF(J)**2
    			End Do
    			IF(SUMPDF>0)THEN
    				HR2N=SUMHR2/SUMPDF
    			ELSE
    				HR2N=0.0
    			End If
    		ELSE
    
    			!*
    			!* NO REFORMATION
    			!* --------------
    			!*
    			!* RMS HEIGHT FOR UNBROKEN WAVES
    			!*
    			!*
    			DAKT=D(I)+EMEAN
    			IF(DAKT<=DMIN.OR.D(I)<0)THEN
!write(*,*) I
    				HN2=(HRMSO**2-EXP(-(GAMBR*DAKT/HRMSO)**2)*(HRMSO**2+&
    					(GAMBR*DAKT)**2))/BMIN
    
    			ELSE
    				HN2=(HRMSO**2-EXP(-(GAMBR*DMIN/HRMSM)**2)*(HRMSO**2+&
    					(GAMBR*DMIN)**2*CGMIN*COSMIN/CTRAN1))/BMIN
    			END IF
    			START=.TRUE.
    
    			!*
    			!*
    			!* SHOAL DISTRIBUTION OF REFORMED WAVES
    			!*
    			!*
    			IF(IPMAX>0)THEN
    				DO J=0,IRMAX-1
    					HRPDF(J)=HRPDF(J)*CTRAN2
    				End Do
    			End If
    
    			!*
    			!*
    			!* RMS HEIGHT FOR REFORMED WAVES (HR2N)
    			!*
    			!*
    			SUMPDF=0.0
    			SUMHR2=0.0
    			HMAX=GAMBR*DAKT
    			jtry=int(hmax/hbeg/ctran2*irmax/4.)
!                 write(*,*) HMAX,jtry
    			if(jtry<1) jtry=1
    			IPWAVE=IRMAX
    			DO J=jtry-1,IRMAX-1
    				IF(HRPDF(J)>HMAX)THEN
    					IPWAVE=J-1
    					EXIT
    				End If
    			End Do
    			ILOOP=MIN(IPWAVE,IPMAX)
    			DO J=0,ILOOP-1
    				SUMPDF=SUMPDF+PDF(J)
    				SUMHR2=SUMHR2+PDF(J)*HRPDF(J)**2
    			End Do
    			IF(SUMPDF>0)THEN
    				HR2N=SUMHR2/SUMPDF
    			ELSE
    				HR2N=0.0
    			End If
    			REFN=SUMPDF
    			DO J=ILOOP,IPMAX-1
    				PDF(J)=0.0
    			End Do
    		End If
    
    		!*
    		!*
    		!* STORE WAVE PROPERTIES
    		!*
    		!*
    		REF=0.5*(REFN+REFO)
    		HR2=0.5*(HR2N+HR2O)
!           write(*,*) HR2, REF
    		BROK(I)=1.-BMIN-REF
    		if(BROK(I)<0.0) BROK(I)=0.0
    		if(.not.FNDLIM.and.BMIN<0.99) then
    			ISWLIM=i
    			fndlim=.true.
    		End If
    
    		if(HN2<=0) HN2=0.0
    		AC=KAPPA*DX(I)/(D(I)+EMEAN)
!           if(H(i+1)/(D(I)+EMEAN)<GAMSTB) AC=0

    		!*
    		!*
    		!* "EQUIVALENT" STABLE ENERGY FLUX
    		!*
    		!*
    		hbratio=dmin*gambr/hrmsm
    		if(hbratio<0.12) then
    			hnsfact=1.289
    		elseif(hbratio>2.48) then
    			hnsfact=1.416
    		else
    			hnsfact=hnsr(nint(hbratio*100.))
    		End If
    		hns2=HN2*hnsfact**2.
!           write(*,*) i,BMIN,Hns2, REF
    		FSTAB=CC*(CG1+CG2)/2*(BMIN*Hns2+REF*HR2+(1-BMIN-REF)*&
    			(D(I)+EMEAN)**2*GAMSTB**2)
    
    		!*
    		!*
    		!* CALCULATE THE ENERGY FLUX
    		!*
    		!*
!           write(*,*) AC,I, DX(I), D(I),EMEAN
    		F(I)=(F(I+1)*(COS(ZW(I+1))-0.5*AC)+AC*FSTAB)/(COS(ZW(I))+0.5*AC)
!           error = ABS(F(I)-F(I+1))/F(I+1)*100
           
    		IF(F(I)<0)THEN
    			ISTOP=I
!                 write(*,*) 'Stopping at 839'
    			Exit MSpace
           
    		End If
!          if(error<0.003) F(I)=F(I+1)
           

    		!*
    		!*
    		!* CALCULATE THE ENERGY DISSIPATION
    		!*
    		!*
    		DISS(I)=KAPPA/(D(I)+EMEAN)**2*(0.5*(F(I)+F(I+1))-FSTAB)
!           write(*,*) DISS(I)
    		IF(DISS(I)<0) DISS(I)=0.0
    
    		!*
    		!* AVOID NEGATIVE ENERGY DISSIPATION IN THE TRANSITION WHEN ALL WAVES CEASE TO BREAK
    		!*
    		IF(F(I)*COS(ZW(I))>F(I+1)*COS(ZW(I+1)))THEN
    			F(I)=F(I+1)*COS(ZW(I+1))/COS(ZW(I))
    		End If
    
    		!*
    		!*
    		!* CALCULATE THE WAVE HEIGHT FROM THE ENERGY FLUX
    		!*
    		!*
    		H(I)=SQRT(F(I)/CC/CG2)
!           write(*,*) I,H(I)
    
    		!*
    		!*
    		!* CALCULATE THE RADIATION STRESS SXX
    		!*
    		!*
    		SXX(I)=CC*H(I)**2.*(CN2*((COS(ZW(I)))**2+1)-0.5)
    
    		!*
    		!*
    		!* CALCULATE THE WIND STRESS FACTOR
    		!*
    		!*
    		CWIND=1.2875*1E-3
    		IF(W>=7.5) CWIND=(0.8+0.065*W)*1E-3
    
    		!*
    		!*
    		!* CALCULATE THE WAVE SET-DOWN OR SET-UP
    		!*
    		!*
    		FORC=(SXX(I+1)-SXX(I))/4./CC+E(I+1)**2+2*D(I)*E(I+1)&
    			+CWIND*1.25*W**2*COS(ZWIND)*DX(I)/4./CC
    		IF(D(I)**2+FORC<0)THEN
    			ISTOP=I
!write(*,*) 'Stopping at 894'
    			Exit MSpace
    		End If
    		E(I)=-D(I)+SQRT(D(I)**2+FORC)
    		IF(DH(I)+E(I)<0.or.DA(I)+E(I)<0)THEN
    			ISTOP=I
!                 write(*,*) I,E(I),da(I),'breaking at 873'
    			Exit MSpace
    		End If
    		if(E(i)-E(i+1)>H(i+1)*1./2.) then
			istop=i
!write(*,*) 'Stopping at 905'
    			Exit MSpace
        End If

            if(I>2) then
    			if(da(i)<0.0.and.(da(i)-da(i-1))<=0) then
    				do jinc=i-2,0,-1
    					if(da(jinc)-da(i)>=da(i)+e(i)) Exit
    					if(da(jinc)<=da(i)) Exit
    				End Do
    				IF(da(jinc)-da(i)>=da(i)+e(i)) THEN
    					do jfix=i,ndx-1
    						if(e(jfix)+da(i)<=0) then
							do jjfix=i-1,jfix-2
								e(jjfix)=-da(i)-(-da(i)-e(jfix))*((x(jjfix)-dx(jjfix)/2)&
									-(x(i)-dx(i)/2))/((x(jfix)-dx(jfix)/2)-(x(i)-dx(i)/2))
							End Do
							istop=i
!write(*,*) 'Stopping at 923'
							Exit MSpace
						End If
					End Do
    				END IF
    			End If
    		End If
		!*
		!*
		!*  UPDATE WAVE PARAMETERS
		!*
		!*
		!*
    		CN1=CN2
    		CG1=CG2
    		L1=L2
    		REFO=REFN
    		HR2O=HR2N

    	End Do MSpace
    	!*
    	!* MAIN LOOP IN SPACE ENDS
    	!* -----------------------
    	!*
    
    	!*
    	!*
    	!*  CALCULATE SETUP AND WAVE HEIGHT IN THE MIDDLE OF THE CELL
    	!*
    	!*
    	DO I=0,NDX-1
    		IF(I<ISTOP) THEN
    			E(I)=0.
    		ELSE
    			E(I)=(E(I)+E(I+1))/2.
    		End If
    	End Do
!write(*,*) 'Deallocating'
    	DeAllocate(F,SXX,DH,PDF,HRPDF)
    	RETURN
     End Subroutine WAVRAN

!***************************************************************************************
     Subroutine TRRAN(HO,ZO,T,VF,X,D,ISTOP,ZW,DISS,BROK,EDHDX,Q,ISWLIM,NDX)
!* ROUTINE FOR CALCULATING THE CROSS-SHORE TRANSPORT RATE UNDER IRREGULAR WAVES
!*
    	Use Memory, ONLY: K,RO,ROS,G,CEQ,CCORR,CRED,D50,SRATIO,BI,LAMM,CC2
      !f2py real(8), intent(aux):: K,RO,ROS,G,CEQ,CCORR,CRED,D50,SRATIO,BI,LAMM,CC2
    	Use MemoryConst, ONLY: GAMBR
     !f2py real(8), intent(aux):: GAMBR
    	Implicit NONE
    	Integer,intent(in):: NDX
    	Integer,intent(inout):: ISTOP,ISWLIM
    	Real(8),intent(inout):: D(0:NDX-1),X(0:NDX-1),DISS(0:NDX-1),BROK(0:NDX-1),HO,&
    		T,Q(0:NDX),VF,ZW(0:NDX),ZO,EDHDX(0:NDX-1)
    	Real(8),allocatable:: SIGTR(:),QB(:),QO(:),CL(:),DUB(:)
    	Real(8):: CEQEFF
    	Real(8):: TRBRK,TRDEC
    	Integer:: J,I,DStat


    	Allocate(SIGTR(0:NDX-1),QB(0:NDX-1),QO(0:NDX-1),CL(0:NDX-1),DUB(0:NDX-1),STAT=DStat)
	!
	! TRANSPORT DIRECTION AND WEIGHTING FACTORS
	!
     do I=0,NDX-1
        	SIGTR(I)=0.0
        	QO(I)=0.0
     end do
    	CALL RANDIR(HO,ZO,T,VF,D,ZW,ISTOP,SIGTR,NDX)
!      ISWLIM =ndx-1
	!
	! TRANSPORT CONTRIBUTION FROM BROKEN WAVES
	!
    	DO I=ISTOP,ISWLIM-1
    		CEQEFF=CEQ-EDHDX(I)
    		IF(CEQEFF<0) CEQEFF=0.0
    		IF(DISS(I)>BROK(I)**CC2*CEQEFF)THEN
    			QB(I)=k*(DISS(I)-BROK(I)**CC2*CEQEFF)
    		ELSE
    			QB(I)=0.0
    		End If
    		!
    		! CALCULATE VALUES AT GRID POINTS ACROSS SHORE BEFORE SUMMING
    		! UP CONTRIBUTIONS FROM OFFSHORE ZONES
    		!
    		!         CALL LINWAV(T,D(I),L,CGB(I),CN)
    		CL(I)=LAMM*(1E3*D50/GAMBR/D(I))**0.47
    
    		IF(I>ISTOP.AND.I<NDX)THEN
    			DUB(I)=(BROK(I-1)-BROK(I+1))/2.
    		ELSEIF(I==ISTOP)THEN
    			DUB(I)=BROK(I)-BROK(I+1)
    		ELSEIF(I==NDX-1)THEN
    			DUB(I)=BROK(I-1)-BROK(I)
    		End If
    
    		!
    		! ONLY INCREASE IN THE RATIO OF BROKEN WAVES PRODUCE TRANSPORT
    		!
!           write(*,*) DUB(I)
    		IF(DUB(I)<0) DUB(I)=0
    
    		qo(i)=0.0
    	End Do
    	do i=ISWLIM,ndx-1
    		qb(i)=0.0
    		qo(i)=0.0
    	End Do
    
    
    	!
    	! TRANSPORT CONTRIBUTION FROM PRE-BREAKING ZONES
    	!
    	DO I=ISTOP,ISWLIM-1
    		trbrk=qb(i)*dub(i)
    		if(trbrk==0.0) cycle
    		DO J=I,ndx-1
    			trdec=trbrk*exp(-cl(i)*(x(j)-x(i)))
    			qo(j)=qo(j)+sigtr(i)*trdec
    			if(trdec<=0.01*trbrk) Exit
    		End Do
    	End Do
    	!
    	! TOTAL TRANSPORT (AT Q-POINTS)
    	!
    	!* Q BOUNDARIES
    	Q(ISTOP)=SIGTR(ISTOP)*QB(ISTOP)+QO(ISTOP)
    	if(q(istop)<0) q(istop)=0.0
    	do I=ISTOP+1,NDX-1
    		Q(I)=0.5*(SIGTR(I-1)*QB(I-1)+QO(I-1)+SIGTR(I)*QB(I)+QO(I))
           
    		if(Q(i)<0) Q(i)=0.0
    	End Do
    	Q(NDX)=Q(NDX-1)
    
    	DeAllocate(SIGTR,QB,QO,CL,DUB)
    	RETURN
     End Subroutine TRRAN

!***************************************************************************************
     Subroutine SWSHDOM(X,D,E,HO,NSWALL,DFS,NRU,NFS,IMIN,BLR,ISWALL,IZERO,LRU,ZO,DRUMAX,&
    	ITIME,IDRUMAX,BROK,dru,NDX)
! CROSS-SHORE TRANSPORT RATES IN THE SWASH ZONE
    	USE MEMORY,ONLY:K,RO,ROS,G,CEQ,CCORR,CRED,DX,DT,XSTART, &
    		D50,SRATIO,BI,CC,PI,LO,CEXO,CEXP,AEXP,CEXR,IELEV,DSURGE,PSURGE,LAMM,CC2
      !f2py real(8), intent(aux):: K,RO,ROS,G,CEQ,CCORR,CRED,DX,DT,XSTART,D50,SRATIO,BI,CC,PI,LO,CEXO,CEXP,AEXP,CEXR,DSURGE,PSURGE,LAMM,CC2
    	!f2py integer, intent(aux):: IELEV
     Implicit NONE
    	Integer,INTENT(IN):: NDX,NSWALL,ITIME,ISWALL
    	Integer,INTENT(INOUT):: NFS,NRU,IDRUMAX,IZERO,IMIN
    	Real(8),INTENT(IN):: X(0:NDX-1),D(0:NDX-1),E(0:NDX),DFS,HO,BROK(0:NDX-1),ZO
    	Real(8),INTENT(INOUT):: LRU,DRUMAX,DRU
    	LOGICAL,INTENT(INOUT):: BLR
    	Integer:: I,ISTART,IHIGH,ILOW,II
    	Real(8):: DMIN,HFACT,DCHECK,SUMVOL,TVOL,ADDVOL,DBROK,BI1
    	!*
    	!* DETERMINE THE RUNUP HEIGHT FROM EMPIRICAL RELATIONSHIP
    	!*
    	HFACT=1.00D0*(COS(ZO))**0.5D0
    	DRU=-1.47D0*HFACT*HO*(TAN(BI)/SQRT(HFACT*HO/LO))**0.79D0  !<-RUNUP HEIGHT *HFACT*HO
    	LRU=-DRU/TAN(BI) !<-RUNUP LENGHT FOR A FORESHORE SLOPE OF BI
    	!*
    	!* COMPUTE "SATURATION VOLUME" FOR RUNUP LIMIT DETERMINATION
    	!*
    	TVOL=0.5*DRU**2./TAN(BI)
    
    	IF(-DRU+DSURGE>DRUMAX) THEN ! MAXIMUM RUNUP HEIGHT THROUGHOUT THE SIMULATION
    		DRUMAX=-DRU+DSURGE
           IDRUMAX=ITIME
    	END IF
!if(itime ==724) write(*,*) DRUMAX

    	!*
    	!* DETERMINE X-DISTANCES TO THE TRANSPORT REGION BOUNDARIES.
    	!* IF A BOUNDARY DEPTH OCCUR IN THE MIDDLE OF A GRID CELL
    	!* THE SEAWARD BOUNDARY NUMBER OF THE CELL IS USED
    	!*
    	NFS=NSWALL
    	DO I=NDX-2,NSWALL,-1
    		!* END OF SURF ZONE
    		IF((D(I)+E(I))<DFS.AND.(D(I+1)+E(I+1))>=DFS)THEN
    			NFS=I+1
    			IF(MINVAL(D(NSWALL:NFS))>0) NFS=NSWALL
    			EXIT
    		END IF
    	End Do
    
    	!*
    	!* LIMIT OF RUNUP
    	!*
    	SUMVOL=0
    	IZERO=0
    	BLR=.FALSE.
    	NRU=-999
!     write(*,*) NFS
    	DO I=NFS,NDX-1
    		IF(D(I)>=0) THEN
    			ISTART=I!-1	!<-- FIRST LOCATION WITH POSITIVE DEPTH (EXCLUDING SETUP)
    			EXIT
    		END IF
    	End Do
    	DO I=ISTART,NSWALL,-1
    		ADDVOL=-D(I)*DX(I)
    		IF(ADDVOL<0.0) THEN
    			IF(SUMVOL==0.0) THEN
    				ADDVOL=0.0
    			ELSE
    				ADDVOL=-ADDVOL
    			END IF
    		END IF
    		IF(SUMVOL==0.0.AND.ADDVOL>0.0) IZERO=I
    		SUMVOL=SUMVOL+ADDVOL
    		IF(SUMVOL>TVOL.OR.D(I)<DRU) THEN
    			NRU=I!+1	!<--	RUN-UP LIMIT, BEING EITHER WHEN THE WATER VOL ~ IDEALISED SWASH VOL OR
    					!			WHEN ELEVATION EXCEEDS THE RUN-UP HEIGHT (NOT SIGN AS LOOKING AT
    					!			DEPTHS HERE, NOT ELEVATIONS, THUS NEG DEPTH IS LESS THEN NEG RUN-UP
    					!			HEIGHT
    					!		RUN-UP LIMIT STARTS FROM WHEN DEPTH (EXCLUDING SETUP) IS POSITIVE
    					!			THAT IS, RUN-DOWN LIMIT IS APPROXIMATED BY MSL (TIDE+SURGE)
    			EXIT
    		END IF
    	End Do
!     write(*,*) 'NRU', NRU, TVOL, DRU
    	IF(NRU==-999) THEN !* RUNUP IS LIMITED BY LANDWARD BOUNDARY IF NO SEAWALL
    		NRU=NSWALL
    		IF(ISWALL==0) BLR=.TRUE.
    	END IF
    
    	IF(NRU>NFS)NRU=NFS
    	DMIN=MINVAL(D(NRU:NFS-1))
    	IMIN=MINLOC(D(NRU:NFS-1),DIM=1)+NRU-1
    	!*
    	!* CALCULATE NEW FORESHORE SLOPE
    	!*
    	IHIGH=NRU!MAX(IMIN,NRU)
    	DO I=NDX-1,NRU,-1
    		DBROK=0.D0
    		IF(BROK(I)>0.5D0) THEN
    			DBROK=D(I)	!<-- DEPTH FIRST POINT COMING IN FROM OFFSHORE THAT HAS 50% WAVES BROKEN
    			EXIT
    		END IF
    	End Do
    	IF(DBROK>0.D0) THEN
    		ILOW=NDX-1
    		DO I=IHIGH,NDX-1
    			IF(D(I)>DBROK) THEN
    				ILOW=I-1	!<--	LOCATION WHEN STARTING FROM THE RUN-UP LIMIT THAT HAS THE SAME DEPTH
    							!			AS THE FIRST OFFSHORE POISTION WITH 50% WAVES BROKEN
    				EXIT
    			END IF
    		End Do
    		IF(ILOW>IHIGH) THEN
    			BI1=ATAN((D(ILOW)-D(IHIGH))/(X(ILOW)-X(IHIGH)))
    			IF(BI1>0.0) BI=BI1
    		END IF
    	END IF
    	RETURN
     End Subroutine SWSHDOM

!***************************************************************************************
     Subroutine EXRUN(NRU,NFS,Q,IMIN,X,D,OVERWASH,IZERO,LRU,BLR,dfs,dru,NDX)
!* THIS ROUTINE DETERMINES THE TRANSPORT RATE ACCORDING TO A LINEAR DECAY FROM THE END
!* OF THE SURF ZONE TO THE RUNUP LIMIT.
!    	use ArrSize
    	Use Memory,ONLY:DX,DT,XSTART,D50,SRATIO,BI,DSURGE
     !f2py real(8), intent(aux):: DX,DT,XSTART,D50,SRATIO,BI,DSURGE
    	Use Memoryq_u,ONLY:rrm
     !f2py real(8), intent(aux):: rrm
    	Implicit NONE
    	Integer,intent(inout):: IZERO,NRU,NFS,IMIN,NDX
    	Real(8),intent(inout):: Q(0:NDX),D(0:NDX-1),X(0:NDX-1),LRU,DFS,dru
    	LOGICAL,intent(inout):: OVERWASH,BLR
    
    	Integer:: I
    	Real(8):: DNM,RRATIO,DBORE,UBORE,SLOPE,ANG,FAC,R,z
    	OVERWASH=.FALSE.
    	!*
    	!* IF Q(NFS)=0 --> NO TRANSPORT IN SWASH ZONE
    	!*
!write(*,*) Q(NFS)
    	IF(Q(NFS)==0.0) THEN
    		BLR=.FALSE.
    		DO I=NRU,NFS-1
    			Q(I)=0.0
    		End Do
    		RETURN
    	End If
!     write(*,*) BLR
    	!*
    	!* CHECK IF OVERWASH OCCURS
    	!*
    	IF(NRU<NFS.AND.NRU<IMIN.AND.Q(NFS)>0.0&
    			& .AND.IZERO>0) THEN
    		OVERWASH=.TRUE.
    	End If
    	!*
    	!* OVERWASH ALGORITHM
    	!*
    	IF(OVERWASH) THEN

    		IF(IZERO>IMIN) THEN
    			DNM=X(IZERO)-X(IMIN)
    		ELSE
    			DNM=DX(IZERO)
    		End If
!           RRATIO=2.*(LRU-DNM)/LRU
    		rratio=(lru-dnm)/(lru-(x(izero)-x(nfs)))
          
!write(*,*) rratio, lru, dnm
           R = -dru+dsurge
           z = -D(IMIN)+dsurge
    		IF(RRATIO<0.0)THEN!.OR.lru < dnm) THEN
!write(*,*) rratio, lru, dnm
    			OVERWASH=.FALSE.
    			BLR=.FALSE.
    			NRU=IMIN
    			DNM=(X(NFS)-DX(NFS)/2)-(X(NRU)-DX(NRU)/2)
    			DO I=NRU+1,NFS-1
    				Q(I)=Q(NFS)*(((X(I)-DX(I)/2)-(X(NRU)-DX(NRU)/2))/DNM)**2
    			End Do
    		Else
!write(*,*) nru
    			dbore=dfs*rratio
    			ubore=2*(9.81*dbore)**0.5
    			Q(IMIN)=-0.005*2*sqrt(2*9.81/R)*(R-z)**2!-(0.005/(dbore**0.5))*ubore**3/9.81!dbore**2!-rrm*ubore**3.
!write(*,*) Q(IMIN), Q(nfs)
    			DNM=(X(NFS)-DX(NFS)/2)-(X(IMIN)-DX(IMIN)/2)
    			DO I=IMIN+1,NFS-1
    				Q(I)=Q(IMIN)+(Q(NFS)-Q(IMIN))*&
    				&         (((X(I)-DX(I)/2)-(X(IMIN)-DX(IMIN)/2))/DNM)**2    
    			End Do
    			DNM=(X(IMIN)-DX(IMIN)/2)-(X(NRU)-DX(NRU)/2)
    			DO I=NRU,IMIN-1
    				Q(I)=Q(IMIN)/(1+0.12*((X(I)-DX(I)/2)-(X(IMIN)-DX(IMIN)/2)))!*((X(I)-DX(I)/2)-(X(NRU)-DX(NRU)/2))/DNM
    			End Do
!write(*,*) Q(IMIN), Q(nru),nru,x(nru)
    		End If
    	ELSE
    		DNM=(X(NFS)-DX(NFS)/2)-(X(NRU)-DX(NRU)/2)
    		DO I=NRU,NFS-1
    			Q(I)=Q(NFS)*(((X(I)-DX(I)/2)-(X(NRU)-DX(NRU)/2))/DNM)**2
!                write(*,*) Q(I)
    		End Do
    	End If
    	!*
    	!*    SLOPE DEPENDENCY ADDED AS FIRST APPROXIMATION
    	!*
    	DO I=NRU+1,NFS-1
    		SLOPE=(D(I-1)-D(I))/((DX(I-1)+DX(I))/2.)
    
    		ANG=-ATAN(SLOPE)
    		FAC=tan(ANG)/tan(BI) !BI - FORESHORE SLOPE IN THE RUNUP ZONE
    		if(FAC<0.) FAC=0.
    		!      FAC=((1-SIN(ANG))/(1+SIN(ANG)))**1.0
    		IF(Q(I)>0) Q(I)=Q(I)*FAC
    		if(Q(I)<0.and.ANG>BI) Q(I)=Q(I)/FAC
    	End Do
    	RETURN
     End Subroutine EXRUN

!***************************************************************************************
     Subroutine CONSAN(NRU,D,Q,DP,QP,NSWALL,PEFAIL,ISWFAIL,SWFAIL,time,istop,iswall,&
    	factl,lunits,IHBOT,IHB,DBOT,SCF,NDX)
!* THIS ROUTINE DETERMINES THE NEW DEPTHS FROM THE SAND CONSERVATION EQUATION
    	Use ArrSize, ONLY: StrLenMax
     !f2py character, intent(aux):: StrLenMax
    	Use Memory, ONLY: DX,DT,XSTART,EPS,BMAX,BAV,IELEV,DSURGE,PSURGE
     !f2py real(8), intent(aux):: DX,DT,XSTART,EPS,BMAX,BAV,DSURGE,PSURGE
     !f2py integer, intent(aux):: IELEV
    	Use Memoryq_u, ONLY: CHECKONLY,ERROCC
     !f2py logical, intent(aux):: CHECKONLY, ERROCC
    	Implicit NONE
    	Integer,intent(in):: NDX
    	Integer,intent(inout):: IHBOT,IHB(0:NDX-1),NSWALL,ISWFAIL,iswall,istop
    	Character,intent(inout):: lunits*4
    	Real(8),intent(inout):: DP(0:NDX-1),D(0:NDX-1),QP(0:NDX),Q(0:NDX),PEFAIL,DBOT(0:NDX-1),SCF,FACTL,time
    	Real(8):: CRANG,DX1,DX2
    	Real(8),allocatable:: QSMO(:)
    	Integer:: I,NRU,NBEG,NEND,ilim
    	LOGICAL:: SWFAIL
    	Character(Len=StrLenMax):: TmpStr

    	Allocate(QSMO(0:NDX))
    	!*
    	!* SMOOTH TRANSPORT RATES
    	!*
    	NBEG=NRU
!     write(*,*) NBEG
    	NEND=NDX-1
    	DO I=NBEG,NEND
    		DX1=DX(I)/(DX(I)+DX(I-1))
    		DX2=DX(I-1)/(DX(I)+DX(I-1))
    		QSMO(I)=(Q(I-1)*DX1+Q(I)+Q(I+1)*DX2)/2
    	End Do
    
    	DO I=NBEG,NEND
    		Q(I)=QSMO(I)
    	End Do
    	!*
    	!* SUBTRACT SURGE BEFORE NEW DEPTHS ARE CALCULATED
    	!*
    	DO I=0,NDX-1
    		D(I)=D(I)-DSURGE
    		DP(I)=DP(I)-PSURGE
    	End Do
    	!*
    	!* SOLVE FOR THE DEPTHS ALONG THE GRID
    	!*
    	IF(IHBOT/=1)THEN
    		DO I=NRU,NDX-1
    			D(I)=DP(I)+DT/DX(I)*0.5*(QP(I+1)-QP(I)+Q(I+1)-Q(I))
    		End Do
    	ELSE
    		!*
    		!* HARD BOTTOM MAY RESTRICT TRANSPORT
    		!*
!    write(*,*) '1524'
    		CALL HRDBOT(NRU,DP,D,QP,Q,IHB,DBOT,SCF,NDX)
    	End If
    	!*
    	!* CHECK FOR SEAWALL FAILURE
    	!*
    	IF(ISWALL==1.AND.(SWFAIL.OR.(ISWFAIL==1.AND.-D(NSWALL)<=PEFAIL))) THEN
    		IF(.NOT.SWFAIL) THEN
    			SWFAIL=.TRUE.
!			Call WriUFRPT('')
    			WRITE(*,'("SEAWALL FAILED AT TIME STEP: ",I0," DUE TO PROFILE ELEVATION="F0.2" ",A)') &
    				nint(TIME/DT),-D(NSWALL)/FACTL,Trim(LUNITS)
!			Call WriUFRPT(TmpStr)
!			Call WriUFRPT('')
     		End If
    		ISWALL=0
    		ISWFAIL=0
    		NSWALL=1
    	End If
    	!*
    	!* CHECK IF THE ANGLE OF REPOSE IS EXCEEDED ANYWHERE ALONG THE PROFILE
    	!* AND IN THAT CASE ALLOW FOR AVALANCHING TO MODIFY THE SLOPE BY
    	!* REDISTRIBUTING SAND TO NEIGHBORING CELLS
    	!*
    	!!
    	!! VARIABLE "TIME" REMOVED FROM CALL STATEMENT
    	!!
    	IF(IHBOT/=1)THEN
    		ilim=nswall+1
    		crang=bmax
    		call maxang(time,d,ilim,crang,NDX)
    		if(ErrOcc) return
    		ilim=max(nswall+1,istop-1)
    		crang=0.30D0
    		CALL MAXANG(time,D,ilim,crang,NDX)
    		if(ErrOcc) return
    	ELSE
    		ILIM=NSWALL+1
    		CRANG=BMAX

    		CALL AVHBOT(TIME,IHB,D,ILIM,CRANG,NDX)
    		If(ErrOcc) return;
    		ILIM=MAX(NSWALL+1,ISTOP-1)
    		CRANG=0.30
    		CALL AVHBOT(TIME,IHB,D,ILIM,CRANG,NDX)
    		If(ErrOcc) return;
    	End If
    	!*
    	!* UPDATE ARRAYS DP AND QP
    	!*
    	DO I=0,NDX-1
    		IF(I>=NRU) QP(I)=Q(I)
    		DP(I)=D(I)
    	End Do
    	QP(NDX)=Q(NDX)
    	!*
    	!* ASSIGN ZERO TRANSPORT RATES TO NODES OUTSIDE CALCULATION DOMAIN
    	!*
    	DO I=0,NRU-1
    		QP(I)=0
    	End Do
    	DeAllocate(QSMO)
    	RETURN
     End Subroutine CONSAN

!***************************************************************************************
     Subroutine HRDBOT(NRU,DP,D,QP,Q,IHB,DBOT,SCF,NDX)
! THIS ROUTINE COMPUTES DEPTH CHANGES IN THE PRESENCE OF
! HARD BOTTOMS
!	Use ArrSize, ONLY: NDX
      !integer, intent(aux):: NDX
    	Use Memory, ONLY: DX,DT,XSTART
      !f2py real(8), intent(aux):: DX,DT,XSTART
    	Implicit NONE
      Integer, intent(in):: NDX
    	Integer,intent(inout):: IHB(0:NDX-1)
    	Real(8),intent(inout):: DP(0:NDX-1),D(0:NDX-1),QP(0:NDX),Q(0:NDX),DBOT(0:NDX-1),SCF
    	Integer,intent(inout):: NRU
    
    	Integer::ICROSS(50),NCROSS,IDIR(50),IS,ISGN,II,IL1,IL2,ILEX
    	Integer:: K,NStart,NStop,NS,I,IHBB,IHBBT,IEXP
    
    	Real(8),allocatable:: QBOT(:),QAKT(:)
    	Real(8):: DTemp,DIR,SUMQ
    
    	Allocate(QBOT(0:NDX),QAKT(0:NDX))
    	!
    	! POTENTIAL TRANSPORT (IGNORING EXPOSED HARD BOTTOMS)
    	!

    	DO I=NRU,NDX
    		QAKT(I)=0.5*(QP(I)+Q(I))
    		QBOT(I)=QAKT(I)
    	End Do
    	!
    	! LOCATE OFFSHORE (+) AND ONSHORE (-) TRANSPORT AREAS
    	!
    	NCROSS=1
    	ICROSS(1)=NRU
    	!
    	! FIND TRANSPORT DIRECTION IN FIRST AREA
    	!
    	DO I=NRU,NDX
    		IF(QAKT(I)/=0)THEN
    			ISGN=NINT(ABS(QAKT(I))/QAKT(I))
    			Exit
    		End If
    	End Do
    	IDIR(1)=ISGN
    	DO I=NRU,NDX
    		IF(QAKT(I)>0) IS=1
    		IF(QAKT(I)<0) IS=-1
    		IF(QAKT(I)==0) IS=ISGN
!write(*,*) ISGN,IS
    		IF(IS/=ISGN)THEN
    			ISGN=IS
    			NCROSS=NCROSS+1
    			ICROSS(NCROSS)=I-1
    			IDIR(NCROSS)=ISGN
    		End If
    	End Do
    	ICROSS(NCROSS+1)=NDX-1
    	!
    	! APPLY SAND CONSERVATION EQUATION AND CHECK FOR EXPOSED HARD
    	! BOTTOMS. CORRECT TRANSPORT RATES IF NECESSARY; DO THE
    	! CORRECTIONS INDEPENDENTLY WITHIN EACH TRANSPORT AREA
    	!
    	IHBBT=0
!write(*,*) NCROSS
    	DO K=0,NCROSS
!write(*,*) K
    		IF(IDIR(K)==1)THEN
    			NSTART=ICROSS(K)
    			NSTOP=ICROSS(K+1)-1
    			NS=1
    			IF(K==NCROSS) NSTOP=NSTOP+1
    		ELSE
    			NSTART=ICROSS(K+1)+1
    			NSTOP=ICROSS(K)+2
    			NS=-1
    			IF(K==1) NSTOP=NSTOP-1
    		End If
    		DIR=Real(IDIR(K))
    		IHBB=0
    		DO I=NSTART,NSTOP,NS
!write(*,*) I
    			IF(NS==1)THEN
    				II=I
    			ELSE
    				II=I-1
    			End If
    			IHB(II)=0
    			D(II)=DP(II)+DT/DX(II)*(QAKT(I+NS)-QAKT(I))*DIR
!write(*,*) D(II)
    			IF(D(II)>DBOT(II))THEN
    				IHBB=1
    				IHBBT=1
    			End If
    			IF(D(II)>=DBOT(II)) IHB(II)=1
    		End Do
    		!
    		! EXPOSED HARD BOTTOM ENCOUNTERED WITHIN THE REGION
    		!
    		IF(IHBB==1)THEN
    			IEXP=0
    			SUMQ=0
    			DO I=NSTART,NSTOP,NS
    				IF(NS==1)THEN
    					II=I
    					IL1=I-1
    					IL2=I
    				ELSE
    					II=I-1
    					IL1=I
    					IL2=I-1
    				End If
    				IF(D(II)<DBOT(II))THEN
    					IF(IEXP==0)THEN
    						QBOT(I+NS)=QAKT(I+NS)
    					ELSE
    						SUMQ=SUMQ+DX(II)
    						QBOT(I+NS)=QAKT(I+NS)+(QBOT(ILEX)-QAKT(I+NS))*EXP(-SCF*SUMQ)
    						IF(ABS(QBOT(I+NS))>ABS(QAKT(I+NS)))QBOT(I+NS)=QAKT(I+NS)
    					End If
    					!
    					! CHECK IF ADDITIONAL HARD BOTTOMS ARE EXPOSED AFTER MODIFYING
    					! TRANSPORT RATES
    					!
    					IF(I>1)THEN
    						IF(IHB(IL1)==1.AND.IHB(IL2)==0)THEN
    							IEXP=1
    							ILEX=IL1
    							SUMQ=SUMQ+DX(II)
    							QBOT(I+NS)=QAKT(I+NS)+(QBOT(ILEX)-QAKT(I+NS))*EXP(-SCF*SUMQ)
    							IF(ABS(QBOT(I+NS))>ABS(QAKT(I+NS)))QBOT(I+NS)=QAKT(I+NS)
    							DTEMP=DP(II)+DT/DX(II)*(QBOT(I+NS)-QBOT(I))*DIR
    							IF(DTEMP>=DBOT(II))THEN
    								QBOT(I+NS)=QBOT(I)+DX(II)/DT*(DBOT(II)-DP(II))*DIR
    								IHB(II)=1
    								IEXP=0
    								SUMQ=0.0
    							End If
    						End If
    					End If
    				ELSE
    					!
    					! CHECK IF TRANSPORT CORRECTIONS START IN THE MIDDLE OF THE GRID
    					!
    					IF(I==NSTART)THEN
    						IF(I==NRU.OR.I==NDX-1)THEN
    							QBOT(I+NS)=QBOT(I)+DX(II)/DT*(DBOT(II)-DP(II))*DIR
    						ELSE
    							!
    							! FOR STARTING CELL IN THE MIDDLE OF THE GRID, DIVIDE TRANSPORTED
    							! AMOUNT OUT FROM THE CELL PROPORTIONALLY
    							!
    							QBOT(I+NS)=QAKT(I+NS)*(DBOT(II)-DP(II))/(D(II)-DP(II))
    						End If
    					ELSE
    						QBOT(I+NS)=QBOT(I)+DX(II)/DT*(DBOT(II)-DP(II))*DIR
    					End If
    				End If
    			End Do
    		End If
    	End Do
    	IF(IHBBT==0) RETURN !.and.NCROSS>0
    	!
    	! NEW DEPTHS
    	!
!     if(IHBBT==0.and.NCROSS==0) then
    	DO I=NRU,NDX-1
!write(*,*) QBOT(i), QAKT(I)
    		D(I)=DP(I)+DT/DX(I)*(QBOT(I+1)-QBOT(I))
    	End Do
!     endif
    	DeAllocate(QBOT,QAKT)
    
    	RETURN
     End Subroutine HRDBOT

!***************************************************************************************
     Subroutine RANDIR(HO,ZO,T,VF,D,ZW,ISTOP,SIGTR,NDX)
! THIS ROUTINE DETERMINES THE LOCAL DIRECTION OF CROSS-SHORE
! TRANSPORT UNDER IRREGULAR WAVES.
    	Use MemoryConst,ONLY: GAMBR,prixDIR
      !f2py real(8), intent(aux):: GAMBR, prixDIR
    	Implicit NONE
    	Integer,intent(inout):: ISTOP
    	Integer,intent(in):: NDX
    	Real(8),intent(inout):: SIGTR(0:NDX-1)
    	Real(8),intent(in):: D(0:NDX-1),HO,T,VF,ZW(0:NDX),ZO
    
    
    	Real(8):: HWT,L,CN,HL,HRMSO,CGO,DAKT,HBO,CGB
    	Real(8):: SIG,ZAKT,E1,E2
    	Integer:: I
    
    	!
    	! DEEPWATER RMS WAVE HEIGHT
    	!
    	HRMSO=HO
    	CGO=0.78065*T
    	!
    	! DETERMINE DEAN PARAMETER AND DEEPWATER STEEPNESS
    	!
    	HWT=HRMSO/VF/T
    	HL=HRMSO/1.5613/T**2
    	!
    	! COMPUTE TRANSPORT DIRECTION AND WEIGHTS ALONG GRID
    	!
    	DAKT=D(NDX-1)
    	SIG=1.0
    	L=0.;CGB=0.;CN=0.
!write(*,*) ISTOP
    	DO I=NDX-1,ISTOP,-1
    		IF(D(I)<DAKT)THEN
    			DAKT=D(I)
    			ZAKT=0.5*(ZW(I)+ZW(I+1))
    			CALL LINWAV(T,DAKT,L,CGB,CN)
    			HBO=SQRT((GAMBR*DAKT)**2*CGB/CGO*COS(ZAKT)/COS(ZO))
    			!
    			! DIFFERENCE BETWEEN EROSIONAL AND ACCRETIONARY WAVES
    			!
    			E1=EXP(-HL/HWT**3/prixDIR)
    			E2=EXP(-(HBO/HRMSO)**2)
    			if(e2<0.0001) then
    				SIG=0.
    			Else
    				SIG=min(1.,2*E1/E2-1.)
    			End If
    		End If
    		SIGTR(I)=SIG
    	End Do
    
    	RETURN
     End Subroutine RANDIR

!***************************************************************************************
     Subroutine EQCOF(CDEAN)
    !* THIS ROUTINE DETERMINES THE COEFFICIENT IN DEAN'S BEACH
    !* PROFILE EQUATION AND THE EQUILIBRIUM ENERGY DISSIPATION 
    ! used in schematisized profile
    	Use Memory,ONLY:K,RO,ROS,G,CEQ,CCORR,CRED,D50,SRATIO,BI
     !f2py real(8), intent(aux):: K,RO,ROS,G,CEQ,CCORR,CRED,D50,SRATIO,BI
    	Implicit NONE
    	Real(8),intent(inout):: CDEAN
    	Real(8):: DMM
    
    	DMM=D50*1.0D3
    	!*
    	!* CALCULATE THE EQUILIBRIUM PROFILE SHAPE FACTOR
    	!*
    	IF(DMM<0.4) CDEAN=0.41*DMM**0.94
    	IF(DMM>=0.4.AND.DMM<10) CDEAN=0.23*DMM**0.32
    	IF(DMM>=10.AND.DMM<40) CDEAN=0.23*DMM**0.28
    	IF(DMM>=40) CDEAN=0.46*DMM**0.11
    	!*
    	!* CALCULATE THE EQUILIBRIUM ENERGY DISSIPATION
    	!*
!write(*,*) CDEAN
!     cdean=0.1
    	CEQ=CCORR*CDEAN**1.5*RO*G**1.5*0.127
    	RETURN
     End Subroutine EQCOF

!***************************************************************************************
     Subroutine LINWAV(T,DD,L,CG,CN)
!* THIS ROUTINE SOLVES THE DISPERSION EQUATION USING A PADE APPROXIMATION
    	Use Memory,ONLY:K,RO,ROS,G,CEQ,CCORR,CRED,CC,PI,LO
     !f2py real(8), intent(aux):: K,RO,ROS,G,CEQ,CCORR,CRED,CC,PI,LO
    	Implicit NONE
    	Real(8),intent(in):: T,DD
    	Real(8),intent(out):: L,CG,CN
    	Real(8):: Y,F,TWOPI,C
	!*
	!* USE A PADE APPROXIMATION TO CALCULATE THE WAVE LENGTH
	!*
    	TWOPI=2*PI
!write(*,*) T
    	Y=(TWOPI/T)**2*DD/G
!write(*,*) F
    	F=Y+1./(1.+Y*(0.66667+Y*(0.3555+Y*(0.16084+Y*(0.06320+ &
    		Y*(0.02174+Y*(0.00654+Y*(0.00171+Y*(0.00039+Y*0.000111 &
    		)))))))))
!write(*,*) F
    	L=TWOPI/SQRT(Y*F/DD**2)
    	C=4*PI*DD/L
    	IF(C<15) THEN
    		CN=0.5*(1+C/SINH(C))
    	ELSE
    		CN=0.5
    	End If
    	CG=CN*L/T
    	RETURN
     End Subroutine LINWAV

!***************************************************************************************
     Subroutine KVISC(TETA,ATETA)
    !* THIS ROUTINE DETERMINES THE CUBIC SPLINE COEFFICIENTS TO BE USED FOR CALCULATING THE
    !* KINEMATIC VISCOSITY. THE COEFFICIENTS ARE STORED IN THE ARRAY ATETA AND THE
    !* CORRESPONDING TEMPERATURES IN THE ARRAY TETA.
    	Use ArrSize
     !f2py integer, intent(aux)::NDX
    	Implicit NONE
    	Real(8),allocatable:: NYVAL(:)
    	Real(8),intent(inout):: TETA(0:NDX-1),ATETA(0:NDX-1,0:3)
    	Integer:: I
    
    	Allocate(NYVAL(0:NDX-1))
    	NYVAL(0)=1.79
    	NYVAL(1)=1.51
    	NYVAL(2)=1.31
    	NYVAL(3)=1.14
    	NYVAL(4)=1.00
    	NYVAL(5)=0.894
    	NYVAL(6)=0.799
    	NYVAL(7)=0.724
    	NYVAL(8)=0.658
    	DO I=0,8
    		TETA(I)=Real(I-1)*5.0
    		NYVAL(I)=NYVAL(I)*1E-6
    	End Do
    	CALL CUBSPL(TETA,NYVAL,8,0,ATETA)
    	DeAllocate(NYVAL)
    	RETURN
     End Subroutine KVISC

!***************************************************************************************
     Subroutine CUBSPL(X,Y,N,IEND,A)
    !* THIS ROUTINE CALCULATES THE COEFFICIENTS FOR THE CUBIC SPLINE POLYNOMIALS
    	use ArrSize
     !f2py integer, intent(aux):: NDX
    	Implicit NONE
    	Integer,intent(inout):: N,IEND
    	Real(8),intent(inout):: X(0:NDX-1),Y(0:NDX-1),A(0:NDX-1,0:3)
    
    	Integer:: NM1,NM2,J,I
    	Real(8),allocatable:: S(:)
    	Real(8):: DX1,DY1,DX2,DY2,DXN2,DXN1
    	Allocate(S(NDX))
    	NM2=N-2
    	NM1=N-1
    	DX1=X(1)-X(0)
    	DY1=(Y(1)-Y(0))/DX1*6.
    	!*
    	!* DETERMINE COEFFICIENTS IN MATRIX FOR SOLUTION OF THE SYSTEM OF EQUATIONS
    	!*
    	DO I=0,NM2
    		DX2=X(I+2)-X(I+1)
    		DY2=(Y(I+2)-Y(I+1))/DX2*6.
    		A(I,1)=DX1
    		A(I,2)=2.*(DX1+DX2)
    		A(I,3)=DX2
    		A(I,4)=DY2-DY1
    		DX1=DX2
    		DY1=DY2
    	End Do
    	!*
    	!* ADJUST MATRIX ACCORDING TO END CONDITIONS
    	!*
    	IF(IEND==1)THEN
    		A(0,1)=A(0,1)+X(1)-X(0)
    		A(NM2,0)=A(NM2,0)+X(N)-X(NM1)
    	ELSEIF(IEND==2)THEN
    		DX1=X(1)-X(0)
    		DX2=X(2)-X(1)
    		A(0,1)=(DX1+DX2)*(DX1+2.*DX2)/DX2
    		A(0,2)=(DX2*DX2-DX1*DX1)/DX2
    		DXN2=X(NM1)-X(NM2)
    		DXN1=X(N)-X(NM1)
    		A(NM2,0)=(DXN2*DXN2-DXN1*DXN1)/DXN2
    		A(NM2,1)=(DXN1+DXN2)*(DXN1+2.*DXN2)/DXN2
    	End If
    	!*
    	!* SOLVE SYSTEM OF EQUATIONS
    	!*
    	DO I=1,NM2
    		A(I,1)=A(I,2)-A(I,0)/A(I-1,1)*A(I-1,2)
    		A(I,3)=A(I,3)-A(I,0)/A(I-1,1)*A(I-1,3)
    	End Do
    	A(NM2,3)=A(NM2,3)/A(NM2,1)
    	DO I=1,NM2
    		J=NM1-I
    		A(J,3)=(A(J,3)-A(J,2)*A(J+1,3))/A(J,1)
    	End Do
    	DO I=0,NM2
    		S(I+1)=A(I,3)
    	End Do
    	!*
    	!* DETERMINE THE SECOND DERIVATIVE AT THE END POINTS
    	!*
    	IF(IEND==0)THEN
    		S(0)=0
    		S(N)=0
    	ELSEIF(IEND==1)THEN
    		S(0)=S(1)
    		S(N)=S(N-1)
    	ELSE
    		S(0)=((DX1+DX2)*S(1)+DX1*S(2))/DX2
    		S(N)=((DXN2+DXN1)*S(NM1)-DXN1*S(NM2))/DXN2
    	End If
    	!*
    	!* DETERMINE THE COEFFICIENTS FOR THE SPLINE POLYNOMIALS
    	!*
    	DO I=0,N-1
    		A(I,0)=(S(I+1)-S(I))/6./(X(I+1)-X(I))
    		A(I,1)=S(I)/2
    		A(I,2)=(Y(I+1)-Y(I))/(X(I+1)-X(I))-(X(I+1)-X(I))*(2*S(I)+S(I+1))/6
    		A(I,3)=Y(I)
    	End Do
    	DeAllocate(S)
    	RETURN
     End Subroutine CUBSPL

!***************************************************************************************
    Subroutine FALVEL(TEMPC,TETA,ATETA,VF)
    !* THIS ROUTINE CALCULATES THE FALL VELOCITY FROM HALLERMEIER'S FORMULA IN THE SPM. THE
    !* FALL VELOCITY IS DETERMINED FOR DIFFERENT GRAIN SIZES AND TEMPERATURES.
    	use ArrSize
     !f2py integer, intent(aux):: NDX
    	Use Memory, ONLY: D50,SRATIO,BI
     !f2py real(8), intent(aux):: D50, SRATIO,BI
    	Implicit NONE
    	Real(8),intent(inout),dimension(0:NDX-1):: TETA
      Real(8), intent(inout)::TEMPC
      Real(8), intent(inout), dimension(0:NDX-1,0:3)::ATETA
      real(8), intent(out)::VF
    	Real(8):: B,NY
    	Integer:: I
    
    	!*
    	!* DETERMINE KINEMATIC VISCOSITY
    	!*
    	DO I=0,7
    		IF(TETA(I)<=TEMPC.AND.TETA(I+1)>=TEMPC)THEN
    			NY=ATETA(I,0)*(TEMPC-TETA(I))**3+ATETA(I,1)*(TEMPC-&
    			&      TETA(I))**2+ATETA(I,2)*(TEMPC-TETA(I))+ATETA(I,3)
    			Exit
    		End If
    	End Do
    	B=(SRATIO-1)*9.81*D50**3/NY**2
    	!*
    	!* CALCULATE FALL VELOCITY
    	!*
    	IF(B<=39)THEN
    		VF=(SRATIO-1)*9.81*D50**2/18/NY
    	ELSEIF(B>39.AND.B<=1E4)THEN
    		VF=((SRATIO-1)*9.81)**0.7*D50**1.1/6/NY**0.4
    	ELSE
    		VF=((SRATIO-1)*9.81*D50/0.91)**0.5
    	End If
    	RETURN
     End Subroutine FALVEL

!***************************************************************************************
     Subroutine SlopeTransPort(istop,NSWALL,EDHDX,D,DX,NDX)
    	Use Memory, ONLY: EPS,K
     !f2py real(8), intent(aux):: EPS,K
    	Implicit NONE
    	Integer,intent(in):: NDX,istop,NSWALL
    	Real(8),intent(inout):: EDHDX(0:NDX-1)
    	Real(8),intent(in):: D(0:NDX-1),DX(0:NDX-1)
    	Integer:: ISTOPI,J

    	if(istop<=0)then
    		istopi=1
    	else
    		istopi=istop
    	End If
    
    	DO J=ISTOPI,NDX-1
    		!*
    		!* DEPTH-DEPENDENT SLOPE COEFFICIENT
    		!*
    		EDHDX(J)=EPS/k*2*(D(J)-D(J-1))/(DX(J)+DX(J-1))
    		IF(NSWALL>=ISTOPI-1) EDHDX(J)=0.0
    	End Do

    	Return
     End Subroutine SlopeTransPort

!***************************************************************************************
     Subroutine RANDOM_I(RPERC,HIN,HOUT)
    !* THIS ROUTINE GENERATES A RANDOM WAVE HEIGHT HOUT WHICH IS UNIFORMLY DISTRIBUTED AROUND
    !* HIN WITH THE OUTER LIMITS HIN +/- RPERC*HIN (0 < RPERC < 1.0)
     Implicit NONE
    	Real(8), intent(in):: RPERC,HIN
     Real(8), intent(out)::HOUT
    	Real(4)::RST
    	CALL RANDOM_NUMBER(RST)
    	HOUT=HIN*(1-RPERC)+dble(RST)*2*RPERC*HIN
    	RETURN
     End Subroutine RANDOM_I

!***************************************************************************************
     Subroutine AddSetUpToDepth(ISTOP,DSMO,E,NDX)
    	Implicit NONE
    	Integer,intent(in):: NDX
    	Integer,intent(inout):: ISTOP
    	Real(8),intent(inout):: DSMO(0:NDX-1)
    	Real(8),intent(in):: E(0:NDX)

    	Integer:: J
    	DO J=NDX-1,ISTOP,-1
    		DSMO(J)=DSMO(J)+E(J)
    		if(DSMO(J)<0) then
    			ISTOP=J+1
    			Exit
    		End If
    	End Do
    	Return
     End Subroutine AddSetUpToDepth

!***************************************************************************************
     Subroutine DOMSMO(DX,D,HO,DSMO,nswall,NDX)
    ! ROUTINE THAT SMOOTHS PROFILE
    	use Memory,ONLY: CC,PI,LO
      !f2py real(8), intent(aux):: CC,PI,LO
    	Implicit NONE
    	Integer,intent(inout):: NDX,NSWALL
    	Real(8),intent(inout), dimension(0:NDX-1):: D,DX
     real(8), intent(out), dimension(0:NDX-1):: DSMO
     real(8), intent(inout):: HO
    	Real(8):: SUM,HBEST
    	Real(8):: DXPREV,SUM2
    	Integer:: MAV,I,J

	!*
	!* FIRST ESTIMATE OF BREAKING WAVE USED FOR SMOOTHING
	!*
    	IF(HO==0.0) THEN
    		HBEST=0.0
    	ELSE
    		HBEST=HO*0.525*(HO/LO)**(-0.25)
    	End If
	!*
	!* DETERMINE NUMBER OF CELLS FOR MAV
	!*
    	DXPREV=0.0
	!*
	!* DETERMINE MAV PROFILE
	!*
    	DO I=NDX-1,nswall,-1
    		IF(DX(I)/=DXPREV) THEN
    			MAV=MAX(INT(3*HBEST/DX(I))+1,5)
    			IF(MOD(MAV,2)==0) MAV=MAV+1
    			DXPREV=DX(I)
!write(*,*) DXPREV, MAV
    		End If
    		IF(I>NDX-MAV/2)THEN
    			DSMO(I)=D(I)
    		ELSEIF(I-MAV/2<nswall)THEN
    			DO J=0,I-1
    				DSMO(J)=D(J)
    			End Do
    		ELSEIF(D(I-MAV/2)>0)THEN
    			SUM=0
    			SUM2=0
    			DO J=I-MAV/2,I+MAV/2
    				SUM=SUM+DX(J)
    				SUM2=SUM2+D(J)*DX(J)
    			End Do
    			DSMO(I)=SUM2/SUM
    		ELSE
    			DO J=0,I
    				DSMO(J)=D(J)
    			End Do
    			Exit
    		End If
    	End Do
    	RETURN
     End Subroutine DOMSMO

!****************************************************************************************
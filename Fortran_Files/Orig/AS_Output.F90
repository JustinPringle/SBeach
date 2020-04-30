!***************************************************************************************
Subroutine WriScr(MSG)
	Use IOHandles, ONLY: fhoSCR
	Implicit NONE
	Character(Len=*),intent(in)::MSG
	Write(fhoSCR,'("'//Trim(MSG)//'")')
	return
End Subroutine WriScr
!***************************************************************************************
Subroutine WriUFLOGRPT(MSG)
	Use IOHandles, ONLY: fhoRPT,fhoLOG
	Implicit NONE
	Character(Len=*),intent(in)::MSG
	Write(fhoRPT,*) Trim(MSG)
	Write(fhoLOG,*) Trim(MSG)
	return
End Subroutine WriUFLOGRPT
!***************************************************************************************
Subroutine WriFILOG(MSG,I)
	Use IOHandles, ONLY: fhoSCR,fhoLOG
	Implicit NONE
	Character(Len=*),intent(in)::MSG
	Integer,intent(in):: I
	Write(fhoSCR,'(X,'//Trim(MSG)//')') I
	Write(fhoLOG,'('//Trim(MSG)//')') I
	return
End Subroutine WriFILOG
!***************************************************************************************
Subroutine WriUFLOG(MSG)
	Use IOHandles, ONLY: fhoSCR,fhoLOG
	Implicit NONE
	Character(Len=*),intent(in)::MSG
	Write(fhoSCR,*) Trim(MSG)
	Write(fhoLOG,*) Trim(MSG)
	return
End Subroutine WriUFLOG
!***************************************************************************************
Subroutine WriUFSCROnly(MSG)
	Use IOHandles, ONLY: fhoSCR
	Implicit NONE
	Character(Len=*),intent(in)::MSG
	Write(fhoSCR,*) Trim(MSG)
	return
End Subroutine WriUFSCROnly
!***************************************************************************************
Subroutine WriUFRPT(MSG)
	Use IOHandles, ONLY: fhoSCR,fhoRPT
	Implicit NONE
	Character(Len=*),intent(in)::MSG
	Write(fhoSCR,*) Trim(MSG)
	Write(fhoRPT,*) Trim(MSG)
	return
End Subroutine WriUFRPT
!***************************************************************************************
Subroutine OpenOutputFiles(NAME,TITLE)
	Use IOHandles
	Implicit NONE
	Character(Len=*),intent(inout)::NAME,TITLE

	open(fhoEXCEL,File=Trim(Name)//Trim(feoEXCEL),status='unknown')

	Open(fhoLOG,FILE=Trim(NAME)//Trim(feoLOG))
	Open(fhoXVR,FILE=Trim(NAME)//Trim(feoXVR))
	Open(fhoPRC,FILE=Trim(NAME)//Trim(feoPRC))
	Open(fhoRPT,FILE=Trim(NAME)//Trim(feoRPT))
	Call TITLEBOX(fhoLOG,NAME,feoLOG,TITLE)
	Call TITLEBOX(fhoXVR,NAME,feoXVR,TITLE)
	Call TITLEBOX(fhoPRC,NAME,feoPRC,TITLE)
	Call TITLEBOX(fhoRPT,NAME,feoRPT,TITLE)
	Call DSHLN(fhoRPT)
	Call DSHLN(fhoLOG)

	return
End Subroutine OpenOutputFiles
!***************************************************************************************
Subroutine CloseAllFiles
	Use IOHandles
	Implicit NONE
	Integer:: IORes

	Close(fhiCFG,IOStat=IORes)
	Close(fhiPRI,IOStat=IORes)
	Close(fhiPRM,IOStat=IORes)
	Close(fhiHDB,IOStat=IORes)
	Close(fhiWAV,IOStat=IORes)
	Close(fhiANG,IOStat=IORes)
	Close(fhiELV,IOStat=IORes)
	Close(fhiWND,IOStat=IORes)
	Close(fhoLOG,IOStat=IORes)
	Close(fhoXVR,IOStat=IORes)
	Close(fhoPRC,IOStat=IORes)
	Close(fhoRPT,IOStat=IORes)
	Close(fhoEXCEL,IOStat=IORes)

	return
End Subroutine CloseAllFiles
!***************************************************************************************
Subroutine WRTX(UNIT,NDX,LUNITS,X,FACTL)
	Use FormatConst, ONLY: FM7050
	Implicit NONE
	Integer,intent(inout):: NDX,UNIT
	Real(8),intent(inout):: X(NDX),FACTL
	Character,intent(inout):: LUNITS*4
	Integer::I
	WRITE(UNIT,5) NDX
	CALL DSHLN(UNIT)
	WRITE(UNIT,10) Trim(LUNITS)
	WRITE(UNIT,FM7050) (X(I)/FACTL,I=1,NDX)
5	FORMAT('NUMBER OF GRID CELLS: ',I0)
10	FORMAT('POSITION OF GRID CELLS RELATIVE TO INITIAL PROFILE ',A4)
	RETURN
End Subroutine WRTX
!***************************************************************************************
Subroutine DSHLN(UNIT)
	Use IOHandles, ONLY: fhoSCR
	Implicit NONE
	Integer,intent(inout):: UNIT
	Character LINE*80
	LINE=Repeat('_',len(LINE))
	IF(UNIT==fhoSCR) THEN
		Call WriUFSCROnly(trim(LINE))
	ELSE
		WRITE(UNIT,'(A)') trim(LINE)
	End If
	RETURN
End Subroutine DSHLN
!***************************************************************************************
Subroutine TITLEBOX(UNIT,NAME,EXT,TITLE)
	Implicit NONE
	Integer,intent(inout):: UNIT
	Character,intent(in)::NAME*(*)
	Character,intent(inout):: EXT*(*),TITLE*(*)
	Character:: STARLN*38
	STARLN='**************************************'
	WRITE(UNIT,5) Trim(STARLN)
	WRITE(UNIT,10) trim(NAME)//Trim(EXT)
	WRITE(UNIT,5) Trim(STARLN)
	WRITE(UNIT,*)
	WRITE(UNIT,15) Trim(TITLE)
	WRITE(UNIT,*)
5	FORMAT(22X,A)
10	FORMAT(22X,'*  Output file: ',A,'  *')
15	FORMAT('RUN: ',A)
	RETURN
End Subroutine TITLEBOX
!***************************************************************************************
Subroutine ConfigReport(NDX,NDT,NWR,IANG,IWAVE,IRAND,ISEED,IWIND,ISWFAIL,ICOMP,IBFILL,&
		IDX,UNITS,NGRID,TPIN,ISWALL,IHBOT,T,DFS,TEMPC,DMEAS,RPERC,ZIN,DXC,W,ZWIND,&
		PEFAIL,WEFAIL,HFAIL,EDP,HIN,DTWAV,DTANG,TELEV,DTELV,DTWND,XSWALL,REFELV,ELV,&
		NDXV,WRI,DXV)
	Use Memory, ONLY: BMAX,D50,DT,EPS,IELEV,K,LAMM
	Use IOHandles, ONLY: fhoRPT
	Use MemoryMain, ONLY: NrNDXV,NrWRI,NrDXV
	Implicit NONE
	Integer,intent(inout):: NDX,NDT,NWR,IANG,IWAVE,IRAND,ISEED,IWIND, &
		ISWFAIL,ICOMP,IBFILL,IDX,UNITS,NGRID,TPIN,ISWALL,IHBOT
	Real(8),intent(inout):: T,DFS,TEMPC,DMEAS,RPERC,ZIN,DXC,W,ZWIND,PEFAIL, &
		WEFAIL,HFAIL,EDP(3),HIN,DTWAV,DTANG,TELEV,DTELV,DTWND,XSWALL,REFELV,ELV(3)
	Integer,intent(inout):: NDXV(NrNDXV),WRI(NrWRI)
	Real(8),intent(inout):: DXV(NrDXV)
	
	Integer:: I

	!* WRITE CONFIGURATION TO REPORT FILE
	Write(fhoRPT,'(T30,"<MODEL CONFIGURATION>")')
	Write(fhoRPT,*)
	Write(fhoRPT,'("INPUT UNITS (SI=1, AMERICAN CUST.=2): "I0)') UNITS
	Write(fhoRPT,'("NUMBER OF CALCULATION CELLS: ",I0)') NDX
	Write(fhoRPT,'("GRID TYPE (CONSTANT=0, VARIABLE=1): ",I0)') IDX
	IF(IDX==0) THEN
		Write(fhoRPT,'("CONSTANT CELL WIDTH: ",F0.1)') DXC
	ELSE
		Write(fhoRPT,'("NUMBER OF GRID CELL REGIONS: ",I0)') NGRID
		DO I=1,NGRID
			Write(fhoRPT,'("NUMBER CELLS AND CELL WIDTH IN REGION ",I0,": "&
				&,I0,",",F0.1)') I,NDXV(I),DXV(I)
		End Do
	End If
	Write(fhoRPT,'("NUMBER OF TIME STEPS AND VALUE OF TIME STEP IN MINUTES:"&
			&," ",I0,",",F0.1)')NDT,DT
	IF(NWR>0) THEN
		WRITE(fhoRPT,'("TIME STEP(S) OF INTERMEDIATE OUTPUT: ",5(I0:", "),/,&
			&T38,4(I0:", "),I0)')(WRI(I),I=1,NWR)
	End If
	IF(NWR == -1 .AND. WRI(1) > 0)THEN
		Write(fhoRPT, '("TIME STEP INTERVAL OF INTERMEDIATE OUTPUT: ", I0)') WRI(1)
	End If
	IF(ICOMP==0) WRITE(fhoRPT,*)'NO COMPARSION WITH MEASURED PROFILE.'
	DO I=1,3
		WRITE(fhoRPT,'("PROFILE ELEVATION CONTOUR ",I0,": ",F0.2)') I,ELV(I)
	End Do
	DO I=1,3
		WRITE(fhoRPT,'("PROFILE EROSION DEPTH ",I0,": ",F0.2)') I,EDP(I)
	End Do
	WRITE(fhoRPT,2295) REFELV
	WRITE(fhoRPT,2210) K*10**6
	WRITE(fhoRPT,2211) EPS
	WRITE(fhoRPT,2212) LAMM
	WRITE(fhoRPT,2213) TEMPC
	WRITE(fhoRPT,*)
	WRITE(fhoRPT,2215) IWAVE
	IF(IWAVE==0) THEN
		WRITE(fhoRPT,2216) HIN,T
	ELSE
		WRITE(fhoRPT,2217) DTWAV
	End If
	WRITE(fhoRPT,2218) IANG
	IF(IANG==0) THEN
		WRITE(fhoRPT,2219) ZIN
	ELSE
		WRITE(fhoRPT,2220) DTANG
	End If
	WRITE(fhoRPT,2221) DMEAS
	IF(IRAND==0) THEN
		WRITE(fhoRPT,*) 'NO RANDOMIZATION OF INPUT WAVE HEIGHT.'
	ELSE
		WRITE(fhoRPT,2222) ISEED,RPERC
	End If
	WRITE(fhoRPT,2223) IELEV
	IF(IELEV==0) THEN
		WRITE(fhoRPT,2224) TELEV
	ELSE
		WRITE(fhoRPT,2225) DTELV
	End If
	WRITE(fhoRPT,2226) IWIND
	IF(IWIND==0) THEN
		WRITE(fhoRPT,2227) W,ZWIND
	ELSE
		WRITE(fhoRPT,2228) DTWND
	End If
	WRITE(fhoRPT,*)
	WRITE(fhoRPT,2229) TPIN
	WRITE(fhoRPT,2231) DFS
	WRITE(fhoRPT,2232) D50
	WRITE(fhoRPT,2233) BMAX
	WRITE(fhoRPT,*)
	IF(IBFILL==0) THEN
		WRITE(fhoRPT,'(X,A)') 'NO BEACH FILL IS PRESENT.'
	ELSE
		WRITE(fhoRPT,'(X,A)') 'A BEACH FILL IS PRESENT.'
	End If
	WRITE(fhoRPT,*)
	IF(ISWALL==0) THEN
		WRITE(fhoRPT,'(X,A)') 'NO SEAWALL IS PRESENT.'
	ELSE
		WRITE(fhoRPT,2234) XSWALL
		IF(ISWFAIL==0) THEN
			WRITE(fhoRPT,'(X,A)') 'SEAWALL FAILURE IS NOT ALLOWED.'
		ELSE
			WRITE(fhoRPT,2235) PEFAIL
			WRITE(fhoRPT,2236) WEFAIL
			WRITE(fhoRPT,2237) HFAIL
		End If
	End If

	WRITE(fhoRPT,*)
	IF(IHBOT==0) THEN
		WRITE(fhoRPT,'(X,A)') 'NO HARD BOTTOM IS PRESENT.'
	ELSE
		WRITE(fhoRPT,'(X,A)') 'HARD BOTTOM IS PRESENT.'
	End If

	CALL DSHLN(fhoRPT)
	WRITE(fhoRPT,*)
	WRITE(fhoRPT,2238)
	2201  FORMAT('INPUT UNITS (SI=1, AMERICAN CUST.=2): 'I0)
	2202  FORMAT('NUMBER OF CALCULATION CELLS: ',I0)
	2203  FORMAT('GRID TYPE (CONSTANT=0, VARIABLE=1): ',I0)
	2204  FORMAT('CONSTANT CELL WIDTH: ',F0.1)
	2205  FORMAT('NUMBER OF GRID CELL REGIONS: ',I0)
	2206  FORMAT('NUMBER CELLS AND CELL WIDTH IN REGION ',I0,': ',I0,',',F0.1)
	2207  FORMAT('NUMBER OF TIME STEPS AND VALUE OF TIME STEP IN MINUTES: ',I0,',',F0.1)
	2275  FORMAT('TIME STEP(S) OF INTERMEDIATE OUTPUT: ',5(I0:', '),/,T38,4(I0:', '),I0)
	2208  FORMAT('PROFILE ELEVATION CONTOUR ',I0,': ',F0.2)
	2209  FORMAT('PROFILE EROSION DEPTH ',I0,': ',F0.2)
	2295  FORMAT('REFERENCE ELEVATION: ',F0.2)
	2210  FORMAT('TRANSPORT RATE COEFFICIENT (m^4/N): ',F0.2,'E-6')
	2211  FORMAT('COEFFICIENT FOR SLOPE DEPENDENT TERM (m^2/s): ',F0.4)
	2212  FORMAT('TRANSPORT RATE DECAY COEFFICIENT MULTIPLIER: ',F0.2)
	2213  FORMAT('WATER TEMPERATURE IN DEGREES C : ',F0.1)
	2214  FORMAT('WAVE TYPE (MONOCHROMATIC=1, IRREGULAR=2): ',I0)
	2215  FORMAT('WAVE HEIGHT AND PERIOD INPUT (CONSTANT=0, VARIABLE=1): ',I0)
	2216  FORMAT('CONSTANT WAVE HEIGHT AND PERIOD: ',F0.2,',',F0.1)
	2217  FORMAT('TIME STEP OF VARIABLE WAVE HEIGHT AND PERIOD INPUT IN MINUTES: ',F0.1)
	2218  FORMAT('WAVE ANGLE INPUT (CONSTANT=0, VARIABLE=1): ',I0)
	2219  FORMAT('CONSTANT WAVE ANGLE: ',F0.1)
	2220  FORMAT('TIME STEP OF VARIABLE WAVE ANGLE INPUT IN MINUTES: ',F0.1)
	2221  FORMAT('WATER DEPTH OF INPUT WAVES (DEEP WATER=0.0): ',F0.1)
	2222  FORMAT('SEED VALUE FOR WAVE HEIGHT RANDOMIZER AND % VARIABILITY: ',I0,',',F0.1)
	2223  FORMAT('TOTAL WATER ELEVATION INPUT (CONSTANT=0, VARIABLE=1): ',I0)
	2224  FORMAT('CONSTANT TOTAL WATER ELEVATION: ',F0.2)
	2225  FORMAT('TIME STEP OF VARIABLE TOTAL WATER ELEVATION INPUT IN MINUTES: ',F0.1)
	2226  FORMAT('WIND SPEED AND ANGLE INPUT (CONSTANT=0, VARIABLE=1): ',I0)
	2227  FORMAT('CONSTANT WIND SPEED AND ANGLE: ',F0.1,',',F0.1)
	2228  FORMAT('TIME STEP OF VARIABLE WIND SPEED AND ANGLE INPUT IN MINUTES: ',F0.1)
	2229  FORMAT('TYPE OF INPUT PROFILE (ARBITRARY=1, SCHEMATIZED=2): 'I0)
	2231  FORMAT('DEPTH CORRESPONDING TO LANDWARD END OF SURF ZONE: ',F0.2)
	2232  FORMAT('EFFECTIVE GRAIN SIZE DIAMETER IN MILLIMETERS: ',F0.2)
	2233  FORMAT('MAXIMUM PROFILE SLOPE PRIOR TO AVALANCHING IN DEGREES: ',F0.1)
	2234  FORMAT('POSITION OF SEAWALL RELATIVE TO INITIAL PROFILE: ',F0.1)
	2235  FORMAT('PROFILE ELEVATION AT SEAWALL WHICH CAUSES FAILURE: ',F0.2)
	2236  FORMAT('WATER ELEVATION AT SEAWALL WHICH CAUSES FAILURE: ',F0.2)
	2237  FORMAT('WAVE HEIGHT AT SEAWALL WHICH CAUSES FAILURE: ',F0.2)
	2238  FORMAT(T32,'<COMPUTED RESULTS>')
	return
End Subroutine ConfigReport
!***************************************************************************************
Subroutine WriteResults(NWR,I,NDX,WRI,D,DI,X,H,E,DSURGE,NRMAX,ELV,XELVR,&
	ISTOP,NDT,ISTMIN,IELV)
	Use IOHandles, ONLY: fhoPRC,fhoEXCEL,fhoXVR
	Use QuantityUnits, ONLY: LUNITS,FACTL
	Use MemoryMain, ONLY: DMX,TDMX,DMN,TDMN,HMX,THMX,EMX,TEMX,DPMX,TDPMX,NrWRI
	Use FormatConst, ONLY: FM7000
	Implicit NONE
	Integer,intent(inout):: WRI(NrWRI)
	Integer,intent(in):: NWR,I,NDX,ISTOP,NDT
	Integer,intent(inout):: ISTMIN,IELV(3)
	Real(8),intent(in):: D(NDX),DI(NDX),X(NDX),H(NDX+1),E(NDX+1),DSURGE,&
		ELV(3)
	Real(8),intent(inout):: XELVR(3)
	Integer,intent(inout):: NRMAX(3)
		
	Logical:: WRT,CycleLoop
	Integer:: J,JX
	Real(8),allocatable:: HMID(:)
	Real(8):: XTEMP
	
	Allocate(HMID(NDX))
	WRT=.FALSE.
	IF(NWR==-1)Then !i.e., regularly space outputs
		If(WRI(1) > 0)Then
			If(MOD(I,WRI(1))==0)Then
				WRT=.TRUE.
			End If
		End If
	Else !i.e., as specified by user
		DO J=1,NWR
			IF(WRI(J)==I)Then
				WRT=.TRUE.
				Exit
			End If
		End Do
	End If
	
	IF(WRT)THEN !Intermediate profile results 
		CALL DSHLN(fhoPRC)
		WRITE(fhoPRC,1902) Trim(LUNITS),I
		WRITE(fhoPRC,FM7000) (-D(J)/FACTL,J=1,NDX)
	End If
	IF(I==NDT)THEN !Final profile result
		CALL DSHLN(fhoPRC)
		WRITE(fhoPRC,1903) Trim(LUNITS)
		WRITE(fhoPRC,FM7000) (-D(J)/FACTL,J=1,NDX)
		do J=1,NDX;write(fhoEXCEL,'(F0.3,2("	"F0.3))') X(j),-DI(j),-D(j);End Do
	End If
1902	FORMAT('PROFILE ELEVATION ',A4,' AFTER TIME STEP: ',I0)
1903	FORMAT('FINAL PROFILE ELEVATION ',A4)
	DO J=ISTOP,NDX
		HMID(J)=(H(J)+H(J+1))/2.
	End Do
	IF(WRT) THEN
		CALL DSHLN(fhoXVR)
		WRITE(fhoXVR,1802) Trim(LUNITS),I,ISTOP
1802	FORMAT('WAVE HEIGHT ',A,' AT TIME STEP: ',I0,/,'FIRST VALUE CORRESPONDS TO CELL',I0)
		WRITE(fhoXVR,FM7000) (HMID(J)/FACTL,J=ISTOP,NDX)
		CALL DSHLN(fhoXVR)
		WRITE(fhoXVR,1803) Trim(LUNITS),I,ISTOP
1803	FORMAT('TOTAL WATER ELEVATION + SETUP ',A,' AT TIME STEP: ',I0,/,'FIRST VALUE CORRESPONDS TO',&
			&' CELL',I0)
		WRITE(fhoXVR,FM7000) ((DSURGE+E(J))/FACTL,J=ISTOP,NDX)
	End If            
	!*
	!* CHECK FOR MAX AND MIN ELEVATION AT EACH POSITION ON PROFILE
	!*
	
	DO J=1,NDX
		IF(D(J)<DMX(J)) THEN
			DMX(J)=D(J)
			TDMX(J)=I
		End If
		IF(D(J)>DMN(J)) THEN
			DMN(J)=D(J)
			TDMN(J)=I
		End If
		IF(J>=ISTOP) THEN
			IF(HMID(J)>=HMX(J)) THEN
				HMX(J)=HMID(J)
				THMX(J)=I
			End If
			IF((DSURGE+E(J))>=EMX(J)) THEN
				EMX(J)=DSURGE+E(J)
				TEMX(J)=I
			End If
			IF((DSURGE+D(J)+E(J))>=DPMX(J)) THEN
				DPMX(J)=DSURGE+D(J)+E(J)
				TDPMX(J)=I
			End If
		End If
	End Do
	IF(ISTOP<ISTMIN) ISTMIN=ISTOP
	!*
	!* CHECK FOR ELEVATION RECESSION
	!*
	DO J=1,3
		CycleLoop=.false.
		IF(IELV(J)/=1) Cycle
		DO JX=NRMAX(J),2,-1
			IF(D(JX)<ELV(J))Then
				CycleLoop=.TRUE.
				Exit
			End If
			IF(D(JX)>=ELV(J).AND.D(JX-1)<=ELV(J)) THEN
				XTEMP=X(JX-1)+(X(JX)-X(JX-1))/(D(JX)-D(JX-1))*(ELV(J)-D(JX-1))
				IF(XTEMP<XELVR(J)) THEN
					XELVR(J)=XTEMP
					NRMAX(J)=JX
				End If
				CycleLoop=.TRUE.
				Exit
			End If
		End Do
		If(.NOT.CycleLoop) IELV(J)=-1
	End Do
	Deallocate(HMID)
	return

End Subroutine	
!***************************************************************************************
Subroutine LoopStatus(MSGT,FLOOD,OVERWASH,BLR,I,NDT)
	Use ArrSize, ONLY: StatusUpdate
	Use IOHandles, ONLY: fhoLOG,fhoSCR
	Implicit NONE
	Character(Len=*),intent(inout):: MSGT
	Logical,intent(in):: FLOOD,OVERWASH,BLR
	Integer,intent(in):: I,NDT
	Character:: MSGO*10,MSGB*24,MSGTO*11,MSGTL*11
	Logical:: SIGNCH
	
	MSGTO=''
	if(mod(I,StatusUpdate)==0) then
		WRITE(fhoSCR,'("+Processing...",I0,"% ",A)') &
			NINT(Real(I)/Real(NDT)*100.),MSGT//MSGO//MSGB
	end if
	IF(MSGTO/=MSGT) THEN
		SIGNCH=.TRUE.
		MSGTL=MSGT
	ELSE
		SIGNCH=.FALSE.
		MSGTL='           '
	End If

	MSGTO=MSGT
	IF(FLOOD) THEN
		WRITE(fhoSCR,915)
		WRITE(fhoLOG,916) I,MSGTL
	ELSE
		IF(OVERWASH) THEN
			MSGO='[OVERWASH]'
		ELSE
			MSGO=Repeat(' ',Len(MSGO))
		End If
		IF(BLR) THEN
			MSGB='[BOUNDARY LIMITED RUNUP]'
		ELSE
			MSGB=Repeat(' ',Len(MSGB))
		End If
		IF(OVERWASH.OR.BLR.OR.SIGNCH) THEN
			WRITE(fhoLOG,901) I,Trim(MSGTL),Trim(MSGO),Trim(MSGB)
		END IF

	End If

900		FORMAT(\\,X,A,A,\)
901		FORMAT('TIME STEP: ',I0,2X,A,2X,A,2X,A)
915		FORMAT(' [ENTIRE PROFILE INUNDATED]',6X,' ',\)
916		FORMAT('TIME STEP: ',I0,2X,A11,2X,'[ENTIRE PROFILE INUNDATED]')

	Return
End Subroutine 
!***************************************************************************************

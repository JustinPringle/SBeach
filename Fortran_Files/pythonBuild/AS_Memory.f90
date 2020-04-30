!
     Module ArrSize
    	Implicit NONE
    	Integer,parameter:: NDXMMax=1000000,StrLenMax=1000
    	Integer:: NDX,StatusUpdate
    	Contains
    	Subroutine InitArrSize
    		Implicit NONE
    		NDX=-999
    		return
    	End Subroutine InitArrSize
     end module ArrSize
!!***************************************************************************************
     Module Memory
    	Implicit NONE
    	Real(8),parameter:: G=9.81d0,AEXP=0.5d0,CEXR=0.1d0,RO=1.d3,ROS=2.6d3,&
    		PI=3.14159d0,CCORR=0.75d0,CRED=0.0d0
    	Real(8):: eps
    	Integer:: IELEV
    	Real(8):: K,CEQ,DT,XSTART,D50,SRATIO,BI,CC,LO,CEXO(15),CEXP,BMAX,BAV,&
    		DSURGE,PSURGE,LAMM,FACTH,MassError
    	Real(8):: CC2
    	Real(8),Allocatable::DX(:)
    	Contains
    	Subroutine InitMemory
    		Implicit NONE
    		eps=0.002d0
    		CC2=0.5
    		CC=G*RO/8.
    		SRATIO=ROS/RO
    		PSURGE=0
    		DSURGE=0
    		FACTH=0.706
    		return
    	End Subroutine InitMemory
	Subroutine DoneMemory
		Implicit NONE
		Integer:: DStat
		if(allocated(DX)) Deallocate(DX,STAT=DStat);Call DeallError(DStat,'DX')

		return
	End Subroutine DoneMemory
     END Module Memory
!!***************************************************************************************
     Module Memoryq_u
    	Implicit NONE
    	Real(8)::rrm
    	Logical::SREMGRID,CHECKONLY,ERROCC
    	Integer:: PriNRV,HBNRV
    	Real(8),allocatable::PriX(:),PriD(:),HBX(:),HBD(:)
    	Contains
    	Subroutine InitMemoryq_u
    		Implicit NONE
    		CHECKONLY=.FALSE.
    		ERROCC=.FALSE.
    		return
    	End Subroutine InitMemoryq_u

	Subroutine DoneMemoryq_u
		Implicit NONE
		Integer:: DStat
		if(allocated(PriX)) Deallocate(PriX,STAT=DStat);Call DeallError(DStat,'PriX')
		if(allocated(PriD)) Deallocate(PriD,STAT=DStat);Call DeallError(DStat,'PriD')
		if(allocated(HBX)) Deallocate(HBX,STAT=DStat);Call DeallError(DStat,'HBX')
		if(allocated(HBD)) Deallocate(HBD,STAT=DStat);Call DeallError(DStat,'HBD')

		return
	End Subroutine DoneMemoryq_u
     END Module Memoryq_u
!!***************************************************************************************
     Module MemoryConst
    	Implicit NONE
    	Real(8),parameter:: KAPPA=0.15d0,GAMSTB=0.4d0,KREF=0.5d0,GAMBR=0.78d0,prixDIR=0.00070
    	Integer,parameter:: IRMAX=100
     END Module MemoryConst
!!***************************************************************************************
     Module MemoryMain
    	Implicit NONE
    	Integer,allocatable:: THMX(:),TEMX(:),TDPMX(:),TDMX(:),TDMN(:)
    	Real(8),allocatable:: DI(:),DSV(:),DMX(:),DPMX(:),HMX(:),VOLCH(:)
    	Integer,allocatable:: IHB(:)
    	Real(8),allocatable:: TETA(:),ATETA(:,:),brok(:),zw(:),DSMO(:),EDHDX(:),DBOT(:),&
    		DMN(:),EMX(:),DP(:),D(:),X(:),Q(:),QP(:),H(:),E(:),DISS(:),XSV(:),DA(:)
    	Integer,allocatable:: NDXV(:),WRI(:)
    	Real(8),allocatable:: DXV(:)
    	Integer:: NrNDXV,NrWRI,NrDXV
!	
    	Contains
    	Subroutine allocateArr
    	 !       Implicit None
    	        Use ArrSize, ONLY: NDX
    	        !f2py integer, intent(aux)::NDX
    	        Use Memory, ONLY: DX
    	        !f2py real(8), intent(aux)::DX
    	        Use Memoryq_u, ONLY: Prix,PriD,HBX,HBD
    	        !f2py real(8), intent(aux)::PriX,PriD,HBX,HBD

    	        Allocate(DX(0:NDX-1),PriX(0:NDX-1),PriD(0:NDX-1),HBX(0:NDX-1),&
                         HBD(0:NDX-1),THMX(0:NDX-1),TEMX(0:NDX-1),&
                         TDPMX(0:NDX-1),TDMX(0:NDX-1),TDMN(0:NDX-1),&
                         IHB(0:NDX-1),DI(0:NDX-1),DP(0:NDX-1),D(0:NDX-1),&
                         X(0:NDX-1),Q(0:NDX),QP(0:NDX),H(0:NDX),&
                         E(0:NDX),DISS(0:NDX-1),XSV(0:NDX-1),&
                         DSV(0:NDX-1),TETA(0:NDX-1),ATETA(0:NDX-1,4),&
                         DMX(0:NDX-1),DMN(0:NDX-1),HMX(0:NDX-1),&
                         EMX(0:NDX-1),DPMX(0:NDX-1),ZW(0:NDX),&
                         VOLCH(0:NDX-1),brok(0:NDX-1),DSMO(0:NDX-1),&
                         EDHDX(0:NDX-1),DBOT(0:NDX-1),DA(0:NDX-1))
                         

    	        return
    	End Subroutine allocateArr
    	
	Subroutine DoneMemoryMain
		Implicit NONE
		Integer:: DStat

		if(allocated(THMX)) Deallocate(THMX,STAT=DStat);Call DeallError(DStat,'THMX')
		if(allocated(TEMX)) Deallocate(TEMX,STAT=DStat);Call DeallError(DStat,'TEMX')
		if(allocated(TDPMX)) Deallocate(TDPMX,STAT=DStat);Call DeallError(DStat,'TDPMX')
		if(allocated(TDMX)) Deallocate(TDMX,STAT=DStat);Call DeallError(DStat,'TDMX')
		if(allocated(TDMN)) Deallocate(TDMN,STAT=DStat);Call DeallError(DStat,'TDMN')
		if(allocated(DI)) Deallocate(DI,STAT=DStat);Call DeallError(DStat,'DI')
		if(allocated(DSV)) Deallocate(DSV,STAT=DStat);Call DeallError(DStat,'DSV')
		if(allocated(DMX)) Deallocate(DMX,STAT=DStat);Call DeallError(DStat,'DMX')
		if(allocated(DPMX)) Deallocate(DPMX,STAT=DStat);Call DeallError(DStat,'DPMX')
		if(allocated(HMX)) Deallocate(HMX,STAT=DStat);Call DeallError(DStat,'HMX')
		if(allocated(VOLCH)) Deallocate(VOLCH,STAT=DStat);Call DeallError(DStat,'VOLCH')
		if(allocated(IHB)) Deallocate(IHB,STAT=DStat);Call DeallError(DStat,'IHB')
		if(allocated(TETA)) Deallocate(TETA,STAT=DStat);Call DeallError(DStat,'TETA')
		if(allocated(ATETA)) Deallocate(ATETA,STAT=DStat);Call DeallError(DStat,'ATETA')
		if(allocated(brok)) Deallocate(brok,STAT=DStat);Call DeallError(DStat,'brok')
		if(allocated(zw)) Deallocate(zw,STAT=DStat);Call DeallError(DStat,'zw')
		if(allocated(DSMO)) Deallocate(DSMO,STAT=DStat);Call DeallError(DStat,'DSMO')
		if(allocated(EDHDX)) Deallocate(EDHDX,STAT=DStat);Call DeallError(DStat,'EDHDX')
		if(allocated(DBOT)) Deallocate(DBOT,STAT=DStat);Call DeallError(DStat,'DBOT')
		if(allocated(DMN)) Deallocate(DMN,STAT=DStat);Call DeallError(DStat,'DMN')
		if(allocated(EMX)) Deallocate(EMX,STAT=DStat);Call DeallError(DStat,'EMX')
		if(allocated(DP)) Deallocate(DP,STAT=DStat);Call DeallError(DStat,'DP')
		if(allocated(D)) Deallocate(D,STAT=DStat);Call DeallError(DStat,'D')
		if(allocated(X)) Deallocate(X,STAT=DStat);Call DeallError(DStat,'X')
		if(allocated(Q)) Deallocate(Q,STAT=DStat);Call DeallError(DStat,'Q')
		if(allocated(QP)) Deallocate(QP,STAT=DStat);Call DeallError(DStat,'QP')
		if(allocated(H)) Deallocate(H,STAT=DStat);Call DeallError(DStat,'H')
		if(allocated(E)) Deallocate(E,STAT=DStat);Call DeallError(DStat,'E')
		if(allocated(DISS)) Deallocate(DISS,STAT=DStat);Call DeallError(DStat,'DISS')
		if(allocated(XSV)) Deallocate(XSV,STAT=DStat);Call DeallError(DStat,'XSV')
		if(allocated(NDXV)) Deallocate(NDXV,STAT=DStat);Call DeallError(DStat,'NDXV')
		if(allocated(WRI)) Deallocate(WRI,STAT=DStat);Call DeallError(DStat,'WRI')
		if(allocated(DXV)) Deallocate(DXV,STAT=DStat);Call DeallError(DStat,'DXV')
		if(allocated(DA)) Deallocate(DA,STAT=DStat);Call DeallError(DStat,'DA')
		
		return
	End Subroutine DoneMemoryMain
     END Module MemoryMain
!!***************************************************************************************
     Module IOHandles
    	Integer,parameter:: &
    		fhiCFG=1010,&
    		fhiPRI=109, &
    		fhiPRM=1016,&
    		fhiHDB=1017,&
    		fhiWAV=1011,&
    		fhiANG=1012,&
    		fhiELV=1013,&
    		fhiWND=1033,&
    		fhoSCR=6, &
    		fhoLOG=1020,&
    		fhoXVR=1018,&
    		fhoPRC=1019,&
    		fhoRPT=1022,&
    		fhoEXCEL=1500
!    	Character(Len=*),parameter:: &
!    		feiCFG='.CFG',&
!    		feiPRI='.PRI', &
!    		feiPRM='.PRM',&
!    		feiHDB='.HDB',&
!    		feiWAV='.WAV',&
!    		feiANG='.ANG',&
!    		feiELV='.ELV',&
!    		feiWND='.WND',&
!    		feoLOG='.LOG',&
!    		feoXVR='.XVR',&
!    		feoPRC='.PRC',&
!    		feoRPT='.RPTNew',&
!    		feoEXCEL='.Excel.AsSupplied'
!		
     End Module IOHandles
!
!
!	
!	

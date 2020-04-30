Module AppVersion
	Use ArrSize, ONLY: StrLenMax
	Implicit NONE
	Character(Len=StrLenMax):: UserID,UserNr,VerNr,VerDate

	Contains
		Subroutine GetUserInfo
			Implicit NONE
			UserID=''
			UserNr=''
			VerNr='Work in progress'
			VerDate='August 2006'
			return
		End Subroutine GetUserInfo
		Subroutine LOGO
			Use IOHandles, ONLY: fhoSCR
			Implicit NONE

			WRITE(fhoSCR,'(X,"SBEACH")')
			WRITE(fhoSCR,'(X,"VERSION: ",A)') Trim(VerNr)
			WRITE(fhoSCR,'(X,"User Nr: ",A)') Trim(UserNr)
			WRITE(fhoSCR,'(X,"User ID: ",A)') Trim(UserID)
			WRITE(fhoSCR,'(X,A)') trim(VerDate)
			RETURN
		End Subroutine LOGO
End Module
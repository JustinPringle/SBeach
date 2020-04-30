Subroutine DeallError(DStat,VarN)
	Implicit NONE
	Integer,intent(in)::DStat
	Character(Len=*),intent(in)::VarN
	
	If(DStat/=0) then
		Call WriScr(' Failed to clear '//VarN//'.')
		Stop
	End if
	Return
End Subroutine DeallError

Subroutine WriScr(MSG)
	Use IOHandles, ONLY: fhoSCR
	Implicit NONE
	Character(Len=*),intent(in)::MSG
	Write(fhoSCR,'("'//Trim(MSG)//'")')
	return
End Subroutine WriScr

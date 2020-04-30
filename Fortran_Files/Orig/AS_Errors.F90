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
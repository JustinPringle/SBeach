module CompareBedSchemes
	Implicit NONE

	Real(8):: SumDdDt(1000,2,4) = 0.0d0


	contains
	Subroutine calctotal(D,DP,DI,NDX,flag)
		Implicit NONE
		Integer,intent(in):: flag,NDX
		Real(8),intent(inout):: D(NDX),DP(NDX)
		Real(8),intent(in):: DI(NDX)
		Integer:: I
		do I=1,NDX
			SumDdDt(I,flag,1)=SumDdDt(I,flag,1)+D(I)-DI(I)
			SumDdDt(I,flag,2)=max(SumDdDt(I,flag,2),abs(D(I)-DI(I)))
			SumDdDt(I,flag,3)=SumDdDt(I,flag,3)+abs(D(I)-DI(I))
			SumDdDt(I,flag,4)=SumDdDt(I,flag,4)+(D(I)-DI(I))**2
		End Do
		I=73;write(1202,'(F0.8)') abs(D(I)-DI(I))
		D=DI
		DP=DI
	End Subroutine

	Subroutine PrintDdDt(NDX,flag)
		Implicit NONE
		Integer,intent(in):: NDX,flag
		Integer:: I

		if(flag/=1) stop 'error with checking'
		do I=1,NDX
			write(1201,'(F0.8,3("	"F0.8))') SumDdDt(I,1,:)
		End Do
	End Subroutine
	
end module


	
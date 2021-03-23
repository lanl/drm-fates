	SUBROUTINE TerminateProgram()
	
	
	use constants
	
	implicit none	
	
	print*,'Fatal error, the program will be terminated. Press any key.'
	write(msgoutfile,*)'Fatal error, the program will be terminated. Press any key.'
	read(*,*)
	stop
	end
!======================================================================================
!======================================================================================
	subroutine OptionNotSupported(str)
	implicit none
	character(*) :: str
	
	print*,trim(str)//': option not supported. Press any key to terminate'
	!read(*,*)
	stop
	end
!===================================================================================
!===================================================================================

	subroutine GenerateQUWinds()

	use diffusion_module
      use flags_module
      
	implicit none
	
	
	integer :: diffiter
	
	if(flag%isfire == 0) print*,'Solving for Mass Consistency'
   if(diff%flag .gt. 1)then
		do diffiter = 1,diff%step
         call divergence
!  call sor3d routine
         call sor3d
!  call Euler equations to get new updated velocities
         call euler
!   call Diffusion operator
         call diffusion
      enddo
   else
!  call divergence routine to calculate divergence of uo field
      call divergence
!  call sor3d routine
      call sor3d
!  call Euler equations to get new updated velocities
!  note that Euler call outfile.f
      call euler
	endif
	END

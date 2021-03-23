      SUBROUTINE DISTURBEDFLOW()

         
         use bld_module
			use landuse_module
         use grid_module
         use winds_module
         
         implicit none
			
         real mag,mago,costheta
			integer i,j,k
			
         quwinds%undisturbed=1.0
         if(bld%number .gt. 0 .or. (landuse%flag .eq. 1 .and. &
					(landuse%veg_flag .eq. 1 .or. landuse%urb_flag .eq. 1)))then
			   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,mag,mago,costheta)
            do k=2,qugrid%nz-1
               do j=1,qugrid%ny-1
                  do i=1,qugrid%nx-1
                     if(bld%icellflag(i,j,k) .gt. 0)then
                        mag=sqrt(((0.5*(quwinds%u(i,j,k)+quwinds%u(i+1,j,k)))**2) &
                                +((0.5*(quwinds%v(i,j,k)+quwinds%v(i,j+1,k)))**2) &
                                +((0.5*(quwinds%w(i,j,k)+quwinds%w(i,j,k+1)))**2))
                        mago=sqrt(quwinds%uint(i,j,k)**2+quwinds%vint(i,j,k)**2)
                        costheta=(                                &
                           quwinds%u(i,j,k)*quwinds%uint(i,j,k)+  &
                           quwinds%v(i,j,k)*quwinds%vint(i,j,k))/(mag*mago)
                        if(abs(mag-mago) .gt. 0.2*mago .or. costheta .lt. 0.8)then
                           quwinds%undisturbed(i,j,k)=0.0
                        endif
                     endif
                  enddo
               enddo
            enddo
            !$OMP END PARALLEL DO
         endif
      END SUBROUTINE
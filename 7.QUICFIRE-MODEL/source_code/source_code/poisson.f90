      subroutine poisson
         
			use sor_module
         use bld_module
         use grid_module
         use winds_module
         
			implicit none
         
			integer iter,i,j,k
         
         do iter=1,10
            !$omp parallel do private(i,j)
            do k= 2,qugrid%nz-1
               do j=2,qugrid%ny-1
                  do i= 2,qugrid%nx-1
                     if(bld%icellflag(i,j,k) .eq. 9 .and. bld%icellflag(i-1,j,k) .eq. 9)then
                        quwinds%uo(i,j,k)=sor%invomegarelax*sor%denom(i,j,k)*                        &
									((sor%bc%e(i,j,k)*quwinds%uo(i+1,j,k)+sor%bc%f(i,j,k)*quwinds%uo(i-1,j,k))        &
                                  +(sor%bc%g(i,j,k)*quwinds%uo(i,j+1,k)+sor%bc%h(i,j,k)*quwinds%uo(i,j-1,k)) &
                                  +(sor%bc%m(i,j,k)*quwinds%uo(i,j,k+1)+sor%bc%n(i,j,k)*quwinds%uo(i,j,k-1)))
                     endif
                     if(bld%icellflag(i,j,k) .eq. 9 .and. bld%icellflag(i,j-1,k) .eq. 9)then
                        quwinds%vo(i,j,k)=sor%invomegarelax*sor%denom(i,j,k)*(					         &
									(sor%bc%e(i,j,k)*quwinds%vo(i+1,j,k)+sor%bc%f(i,j,k)*quwinds%vo(i-1,j,k))         &
                           +(sor%bc%g(i,j,k)*quwinds%vo(i,j+1,k)+sor%bc%h(i,j,k)*quwinds%vo(i,j-1,k))        &
                           +(sor%bc%m(i,j,k)*quwinds%vo(i,j,k+1)+sor%bc%n(i,j,k)*quwinds%vo(i,j,k-1)))
                     endif
                     if(bld%icellflag(i,j,k) .eq. 9 .and. bld%icellflag(i,j,k-1) .eq. 9)then
                        quwinds%wo(i,j,k)=sor%invomegarelax*sor%denom(i,j,k)*(                       &
									(sor%bc%e(i,j,k)*quwinds%wo(i+1,j,k)+sor%bc%f(i,j,k)*quwinds%wo(i-1,j,k))         &
                           + (sor%bc%g(i,j,k)*quwinds%wo(i,j+1,k)+sor%bc%h(i,j,k)*quwinds%wo(i,j-1,k))       &
                           + (sor%bc%m(i,j,k)*quwinds%wo(i,j,k+1)+sor%bc%n(i,j,k)*quwinds%wo(i,j,k-1)))
                     endif
                  enddo
               enddo
            enddo
            !$omp end parallel do
         enddo
         return
      end

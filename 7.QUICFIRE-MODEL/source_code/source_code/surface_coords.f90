!                                 Notice
!  This program was prepared by the University of California (University)
!  under Contract W-7405-ENG-36 with the U.S. Department of Energy (DOE).
!  All rights in the program are reserved by DOE on behalf of the Government
!  and the University pursuant to the contract. You are authorized to use
!  this program for Government purposes but it is not to be released or
!  distributed to the public.
!  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,
!  NOR THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES,
!  MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY
!  OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS, OF
!  ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
!  THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
      subroutine sufacecoords
          
			use file_handling_module
         use bld_module
         use grid_module
         
         implicit none
         
			integer :: i,j,k
         integer idxtotal,idx,v_idx,ncells
         real dxdy
         real, allocatable :: x(:),y(:),dxdz(:),dydz(:)
         real, allocatable :: vertData(:,:),areaData(:)
         integer, allocatable :: faceData(:,:),cellData(:,:)
         
         ncells=int(0.5*real(qugrid%nx*qugrid%ny*qugrid%nz)) ! make a conservative guess as to the necessary size of arrays
         dxdy=qugrid%dx*qugrid%dy
         allocate(x(qugrid%nx),y(qugrid%ny),dxdz(qugrid%nz-1),dydz(qugrid%nz-1))
         allocate(cellData(ncells,4),faceData(ncells,4),vertData(4*ncells,3),areaData(ncells))
         !$omp parallel do default(shared) private(i)
         do i=1,qugrid%nx
            x(i)=(real(i)-1.)*qugrid%dx
         enddo
         !$omp end parallel do
         !$omp parallel do default(shared) private(j)
         do j=1,qugrid%ny
            y(j)=(real(j)-1.)*qugrid%dy
         enddo
         !$omp end parallel do
         !$omp parallel do default(shared) private(k)
         do k=2,qugrid%nz-1
            dxdz(k)=qugrid%dz_array(k)*qugrid%dx
            dydz(k)=qugrid%dz_array(k)*qugrid%dy
         enddo
         !$omp end parallel do
         k=2
         idxtotal=0
         !!$omp parallel do default(shared) private(i,j,idx,v_idx)
         do j=1,qugrid%ny-1
            do i=1,qugrid%nx-1
               if(bld%icellflag(i,j,k) .gt. 0)then
                  !!$omp critical
                  idxtotal=idxtotal+1
                  idx=idxtotal
                  !!$omp end critical
                  v_idx=4*(idx-1)
                  cellData(idx,1)=i
                  cellData(idx,2)=j
                  cellData(idx,3)=k
                  areaData(idx)=dxdy
                  cellData(idx,4)=3
                  faceData(idx,1)=v_idx+1
                  faceData(idx,2)=v_idx+2
                  faceData(idx,3)=v_idx+3
                  faceData(idx,4)=v_idx+4
                  vertData(v_idx+1,1)=x(i)
                  vertData(v_idx+1,2)=y(j)
                  vertData(v_idx+1,3)=qugrid%z(k-1)
                  vertData(v_idx+2,1)=x(i+1)
                  vertData(v_idx+2,2)=y(j)
                  vertData(v_idx+2,3)=qugrid%z(k-1)
                  vertData(v_idx+3,1)=x(i+1)
                  vertData(v_idx+3,2)=y(j+1)
                  vertData(v_idx+3,3)=qugrid%z(k-1)
                  vertData(v_idx+4,1)=x(i)
                  vertData(v_idx+4,2)=y(j+1)
                  vertData(v_idx+4,3)=qugrid%z(k-1)
               endif
            enddo
         enddo
         !!$omp end parallel do
         !!$omp parallel do default(shared) private(i,j,k,idx,v_idx)
         do j=2,qugrid%ny-2
            do k=2,qugrid%nz-1
               do i=2,qugrid%nx-2
                  if(bld%icellflag(i,j,k) .eq. 0)then
                     if(bld%icellflag(i,j,k+1) .gt. 0)then
                         !!$omp critical
                         idxtotal=idxtotal+1
                         idx=idxtotal
                         !!$omp end critical
                         v_idx=4*(idx-1)
                         cellData(idx,1)=i
                         cellData(idx,2)=j
                         cellData(idx,3)=k+1
                         areaData(idx)=dxdy
                         cellData(idx,4)=3
                         faceData(idx,1)=v_idx+1
                         faceData(idx,2)=v_idx+2
                         faceData(idx,3)=v_idx+3
                         faceData(idx,4)=v_idx+4
                         vertData(v_idx+1,1)=x(i)
                         vertData(v_idx+1,2)=y(j)
                         vertData(v_idx+1,3)=qugrid%z(k)
                         vertData(v_idx+2,1)=x(i+1)
                         vertData(v_idx+2,2)=y(j)
                         vertData(v_idx+2,3)=qugrid%z(k)
                         vertData(v_idx+3,1)=x(i+1)
                         vertData(v_idx+3,2)=y(j+1)
                         vertData(v_idx+3,3)=qugrid%z(k)
                         vertData(v_idx+4,1)=x(i)
                         vertData(v_idx+4,2)=y(j+1)
                         vertData(v_idx+4,3)=qugrid%z(k)
                     endif
                     if(bld%icellflag(i,j,k-1) .gt. 0)then
                         !!$omp critical
                         idxtotal=idxtotal+1
                         idx=idxtotal
                         !!$omp end critical
                         v_idx=4*(idx-1)
                         cellData(idx,1)=i
                         cellData(idx,2)=j
                         cellData(idx,3)=k-1
                         areaData(idx)=dxdy
                         cellData(idx,4)=3
                         faceData(idx,1)=v_idx+1
                         faceData(idx,2)=v_idx+2
                         faceData(idx,3)=v_idx+3
                         faceData(idx,4)=v_idx+4
                         vertData(v_idx+1,1)=x(i)
                         vertData(v_idx+1,2)=y(j)
                         vertData(v_idx+1,3)=qugrid%z(k-1)
                         vertData(v_idx+2,1)=x(i+1)
                         vertData(v_idx+2,2)=y(j)
                         vertData(v_idx+2,3)=qugrid%z(k-1)
                         vertData(v_idx+3,1)=x(i+1)
                         vertData(v_idx+3,2)=y(j+1)
                         vertData(v_idx+3,3)=qugrid%z(k-1)
                         vertData(v_idx+4,1)=x(i)
                         vertData(v_idx+4,2)=y(j+1)
                         vertData(v_idx+4,3)=qugrid%z(k-1)
                     endif
                     if(bld%icellflag(i+1,j,k) .gt. 0)then
                         !!$omp critical
                         idxtotal=idxtotal+1
                         idx=idxtotal
                         !!$omp end critical
                         v_idx=4*(idx-1)
                         cellData(idx,1)=i+1
                         cellData(idx,2)=j
                         cellData(idx,3)=k
                         areaData(idx)=dydz(k)
                         cellData(idx,4)=1
                         faceData(idx,1)=v_idx+1
                         faceData(idx,2)=v_idx+2
                         faceData(idx,3)=v_idx+3
                         faceData(idx,4)=v_idx+4
                         vertData(v_idx+1,1)=x(i+1)
                         vertData(v_idx+1,2)=y(j)
                         vertData(v_idx+1,3)=qugrid%z(k-1)
                         vertData(v_idx+2,1)=x(i+1)
                         vertData(v_idx+2,2)=y(j)
                         vertData(v_idx+2,3)=qugrid%z(k)
                         vertData(v_idx+3,1)=x(i+1)
                         vertData(v_idx+3,2)=y(j+1)
                         vertData(v_idx+3,3)=qugrid%z(k)
                         vertData(v_idx+4,1)=x(i+1)
                         vertData(v_idx+4,2)=y(j+1)
                         vertData(v_idx+4,3)=qugrid%z(k-1)
                     endif
                     if(bld%icellflag(i-1,j,k) .gt. 0)then
                         !!$omp critical
                         idxtotal=idxtotal+1
                         idx=idxtotal
                         !!$omp end critical
                         v_idx=4*(idx-1)
                         cellData(idx,1)=i-1
                         cellData(idx,2)=j
                         cellData(idx,3)=k
                         areaData(idx)=dydz(k)
                         cellData(idx,4)=1
                         faceData(idx,1)=v_idx+1
                         faceData(idx,2)=v_idx+2
                         faceData(idx,3)=v_idx+3
                         faceData(idx,4)=v_idx+4
                         vertData(v_idx+1,1)=x(i)
                         vertData(v_idx+1,2)=y(j)
                         vertData(v_idx+1,3)=qugrid%z(k-1)
                         vertData(v_idx+2,1)=x(i)
                         vertData(v_idx+2,2)=y(j)
                         vertData(v_idx+2,3)=qugrid%z(k)
                         vertData(v_idx+3,1)=x(i)
                         vertData(v_idx+3,2)=y(j+1)
                         vertData(v_idx+3,3)=qugrid%z(k)
                         vertData(v_idx+4,1)=x(i)
                         vertData(v_idx+4,2)=y(j+1)
                         vertData(v_idx+4,3)=qugrid%z(k-1)
                     endif
                     if(bld%icellflag(i,j+1,k) .gt. 0)then
                         !!$omp critical
                         idxtotal=idxtotal+1
                         idx=idxtotal
                         !!$omp end critical
                         v_idx=4*(idx-1)
                         cellData(idx,1)=i
                         cellData(idx,2)=j+1
                         cellData(idx,3)=k
                         areaData(idx)=dxdz(k)
                         cellData(idx,4)=2
                         faceData(idx,1)=v_idx+1
                         faceData(idx,2)=v_idx+2
                         faceData(idx,3)=v_idx+3
                         faceData(idx,4)=v_idx+4
                         vertData(v_idx+1,1)=x(i)
                         vertData(v_idx+1,2)=y(j+1)
                         vertData(v_idx+1,3)=qugrid%z(k-1)
                         vertData(v_idx+2,1)=x(i)
                         vertData(v_idx+2,2)=y(j+1)
                         vertData(v_idx+2,3)=qugrid%z(k)
                         vertData(v_idx+3,1)=x(i+1)
                         vertData(v_idx+3,2)=y(j+1)
                         vertData(v_idx+3,3)=qugrid%z(k)
                         vertData(v_idx+4,1)=x(i+1)
                         vertData(v_idx+4,2)=y(j+1)
                         vertData(v_idx+4,3)=qugrid%z(k-1)
                     endif
                     if(bld%icellflag(i,j-1,k) .gt. 0)then
                         !!$omp critical
                         idxtotal=idxtotal+1
                         idx=idxtotal
                         !!$omp end critical
                         v_idx=4*(idx-1)
                         cellData(idx,1)=i
                         cellData(idx,2)=j-1
                         cellData(idx,3)=k
                         areaData(idx)=dxdz(k)
                         cellData(idx,4)=2
                         faceData(idx,1)=v_idx+1
                         faceData(idx,2)=v_idx+2
                         faceData(idx,3)=v_idx+3
                         faceData(idx,4)=v_idx+4
                         vertData(v_idx+1,1)=x(i)
                         vertData(v_idx+1,2)=y(j)
                         vertData(v_idx+1,3)=qugrid%z(k-1)
                         vertData(v_idx+2,1)=x(i)
                         vertData(v_idx+2,2)=y(j)
                         vertData(v_idx+2,3)=qugrid%z(k)
                         vertData(v_idx+3,1)=x(i+1)
                         vertData(v_idx+3,2)=y(j)
                         vertData(v_idx+3,3)=qugrid%z(k)
                         vertData(v_idx+4,1)=x(i+1)
                         vertData(v_idx+4,2)=y(j)
                         vertData(v_idx+4,3)=qugrid%z(k-1)
                     endif
                  endif
               enddo
            enddo
         enddo
         !!$omp end parallel do
         open(unit=ID_FILE_QU_SURFACE,file=TRIM(workingDirectory)// &
				TRIM(fileSeparator)//"QU_surface.bin",form='unformatted',status="unknown")
         write(ID_FILE_QU_SURFACE)idxtotal
         write(ID_FILE_QU_SURFACE)((cellData(i,j),i=1,idxtotal),j=1,4)
         write(ID_FILE_QU_SURFACE)(areaData(i),i=1,idxtotal)
         write(ID_FILE_QU_SURFACE)((faceData(i,j),i=1,idxtotal),j=1,4)
         write(ID_FILE_QU_SURFACE)((vertData(i,j),i=1,4*idxtotal),j=1,3)
         close(ID_FILE_QU_SURFACE)
         deallocate(x,y,dxdz,dydz)
         deallocate(cellData,faceData,vertData)
      end
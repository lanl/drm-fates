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
      subroutine pentagon(ibuild)
!************************************************************************
! pentagon - define pentagon geometry	
! 	  - called by defbuild.f90					
!	  - calls none							
!
!ERP 10/22/2004			
!************************************************************************
         
         use constants
         use bld_module
         use flags_module
         use grid_module
         use winds_module
         
         implicit none
			
			integer,intent(IN) :: ibuild
			integer :: i,j,k
         real x1,y1,x2,y2,x3,y3,x4,y4,x5,y5			
         real x1in,y1in,x2in,y2in,x3in,y3in,x4in,y4in,x5in,y5in
         real xco,yco,uo_h,vo_h,radius_in,radius_out
! local variables added for pentagon flow field man 4/27/05
         integer pentIstart,pentIstop,pentJstart,pentJstop,in_wake,in_court
         integer in_courtwake_u,in_courtwake_v !,dummy_idx
         real norm1,norm2,norm3,norm4,norm5,courtwake_exponent
         real upwind_dir,L_courtwake
         real ret_coeff_court,ret_coeff_courtwake,ret_coeff_roof
         real point_coeff,flat_coeff,expanding_velocity
         real xpent,ypent,xpentu,ypentu,xpentv,ypentv,dNu,dNv
         real zb,delta_y,x_eff,x_nearu,x_nearv,x_u,y_u,x_v,y_v
         real upwind_rel,xw1,yw1,thetaw1,xw2,yw2,thetaw2,ynorm,xnorm
         real xw3,yw3,thetaw3,xw4,yw4,thetaw4,normw1,normw2,normw3,tol
         real L_pentwake,vel_mag,ret_coeff_wake,Ur_u,Vr_u,Ur_v,Vr_v
         real farwake_factor,farwake_exponent,R_u,R_v,theta_u,theta_v
         real zeta_ew,farwake_velocity,u_ref,v_ref
         real ew_start,ew_denom_vert,ew_denom_horz,xnorm2,xnorm3 !,zb_eff
         real Ht_eff,Wp_eff,xf,norm_f,xpf,xf1,xf2,yf1,yf2
         real, allocatable :: uz_ref(:),vz_ref(:)
! new Pentagon Shaped building erp 10/22/04
! bb is the radius of the circle circumscribing the pentagon
! xco and yco are the coordinates of the center of the ellipse
         if(bld%bb(ibuild) .gt. 0.)then
            print*,'Applying Pentagon Building Flow Parameterizations'
            allocate(uz_ref(nz),vz_ref(nz))
! the pentagon was read in with bld%xfo and bld%yfo equalling the center of the
! circle forming the pentagon rename now.
            xco = bld%xfo(ibuild) + bld%bb(ibuild)!CENTER of pentagon
            yco = bld%yfo(ibuild)
! set reference velocities man 1/12/05
            uo_h=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
            vo_h=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
! find upwind direction and determine the type of flow regime
            upwind_dir=atan2(vo_h,uo_h)
            point_coeff=abs(sin(2.5*(upwind_dir-bld%gamma(ibuild)-pi/2.)))
            flat_coeff=abs(cos(2.5*(upwind_dir-bld%gamma(ibuild)-pi/2.)))
            upwind_rel=upwind_dir-bld%gamma(ibuild)
            if(upwind_rel .gt. pi)upwind_rel=upwind_rel-2*pi
            if(upwind_rel .le. -pi)upwind_rel=upwind_rel+2*pi
! define the outer 5 corners of the pentagon 
            x1= bld%bb(ibuild)*cos(18*pi/180)
            y1= bld%bb(ibuild)*sin(18*pi/180)
            x2= bld%bb(ibuild)*cos(-54*pi/180)
            y2= bld%bb(ibuild)*sin(-54*pi/180)
            x3= bld%bb(ibuild)*cos(-126*pi/180)
            y3= bld%bb(ibuild)*sin(-126*pi/180)
            x4= bld%bb(ibuild)*cos(162*pi/180)
            y4= bld%bb(ibuild)*sin(162*pi/180)
            x5= bld%bb(ibuild)*cos(pi/2)
            y5= bld%bb(ibuild)*sin(pi/2)
! Define the 5 corners of the inner pentagon pentagon
! The inner Pentagon Radius is smaller than the outer by a factor
! or 2.5
! Find the corners that will define the limits of the wake
            x1in= (bld%bb(ibuild)/2.5)*cos(18*pi/180)
            y1in= (bld%bb(ibuild)/2.5)*sin(18*pi/180)
            x2in= (bld%bb(ibuild)/2.5)*cos(-54*pi/180)
            y2in= (bld%bb(ibuild)/2.5)*sin(-54*pi/180)
            x3in= (bld%bb(ibuild)/2.5)*cos(-126*pi/180)
            y3in= (bld%bb(ibuild)/2.5)*sin(-126*pi/180)
            x4in= (bld%bb(ibuild)/2.5)*cos(162*pi/180)
            y4in= (bld%bb(ibuild)/2.5)*sin(162*pi/180)
            x5in= (bld%bb(ibuild)/2.5)*cos(pi/2)
            y5in= (bld%bb(ibuild)/2.5)*sin(pi/2)
! Define the angles of the normal vectors of the inner pentagon faces
            norm1=(-126.*pi/180.)+bld%gamma(ibuild)
            norm2=(162.*pi/180.)+bld%gamma(ibuild)
            norm3=(pi/2.)+bld%gamma(ibuild)
            norm4=(18.*pi/180.)+bld%gamma(ibuild)
            norm5=(-54.*pi/180.)+bld%gamma(ibuild)
! Define the extent of the wake according to the wind direction relative to the pentagon
            tol=0.0001
            if(upwind_rel .le. pi+tol .and.   &
                   upwind_rel .ge. 144*pi/180-tol)then
               normw1=norm5
               normw2=0.
               normw3=norm4
               norm_f=norm2
               xw1=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw1=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw2=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw3=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw4=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
               xf1=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yf1=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               xf2=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yf2=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
            elseif(upwind_rel .lt. 144*pi/180-tol .and.   &
                          upwind_rel .gt. 108*pi/180+tol)then
               normw1=norm1
               normw2=norm5
               normw3=norm4
               xw1=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw1=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw2=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw3=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw4=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
            elseif(upwind_rel .le. 108*pi/180+tol .and.   &
                          upwind_rel .ge. 72*pi/180-tol)then
               normw1=norm1
               normw2=0.
               normw3=norm5
               norm_f=norm3
               xw1=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw1=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw2=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw3=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw4=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
               xf1=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yf1=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               xf2=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yf2=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
            elseif(upwind_rel .lt. 72*pi/180-tol .and.   &
                          upwind_rel .gt. 36*pi/180+tol)then
               normw1=norm2
               normw2=norm1
               normw3=norm5
               xw1=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw1=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw2=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw3=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw4=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
            elseif(upwind_rel .le. 36*pi/180+tol .and.   &
                          upwind_rel .ge. 0.-tol)then
               normw1=norm2
               normw2=0.
               normw3=norm1
               norm_f=norm4
               xw1=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw1=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw2=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw3=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw4=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
               xf1=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yf1=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               xf2=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yf2=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
            elseif(upwind_rel .lt. 0.-tol .and.   &
                          upwind_rel .gt. -36*pi/180+tol)then
               normw1=norm3
               normw2=norm2
               normw3=norm1
               xw1=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw1=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw2=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw3=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw4=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
            elseif(upwind_rel .le. -36*pi/180+tol .and.   &
                          upwind_rel .ge. -72*pi/180-tol)then
               normw1=norm3
               normw2=0.
               normw3=norm2
               norm_f=norm5
               xw1=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw1=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw2=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw3=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw4=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
               xf1=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yf1=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               xf2=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yf2=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
            elseif(upwind_rel .lt. -72*pi/180-tol .and.   &
                          upwind_rel .gt. -108*pi/180+tol)then
               normw1=norm4
               normw2=norm3
               normw3=norm2
               xw1=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw1=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw2=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw3=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yw4=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
            elseif(upwind_rel .le. -108*pi/180+tol .and.   &
                          upwind_rel .ge. -144*pi/180-tol)then
               normw1=norm4
               normw2=0.
               normw3=norm3
               norm_f=norm1
               xw1=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw1=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw2=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw3=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw4=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
               xf1=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yf1=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               xf2=x1*cos(upwind_rel)+y1*sin(upwind_rel)
               yf2=-x1*sin(upwind_rel)+y1*cos(upwind_rel)
            else
               normw1=norm5
               normw2=norm4
               normw3=norm3
               xw1=x5*cos(upwind_rel)+y5*sin(upwind_rel)
               yw1=-x5*sin(upwind_rel)+y5*cos(upwind_rel)
               thetaw1=atan2(yw1,xw1)
               xw2=x4*cos(upwind_rel)+y4*sin(upwind_rel)
               yw2=-x4*sin(upwind_rel)+y4*cos(upwind_rel)
               thetaw2=atan2(yw2,xw2)
               xw3=x3*cos(upwind_rel)+y3*sin(upwind_rel)
               yw3=-x3*sin(upwind_rel)+y3*cos(upwind_rel)
               thetaw3=atan2(yw3,xw3)
               xw4=x2*cos(upwind_rel)+y2*sin(upwind_rel)
               yw4=-x2*sin(upwind_rel)+y2*cos(upwind_rel)
               thetaw4=atan2(yw4,xw4)
            endif
!Calculate inner circumscribed circle radius erp 1/05/05
            radius_out= bld%bb(ibuild)
            radius_in = bld%bb(ibuild)*cos(pi/5)
            bld%bb(ibuild)=radius_out
            bld%Wti(ibuild)=2.*radius_out !set minor axis to input Width
            bld%Lti(ibuild)=2.*radius_out !set major axis to input Length
!now update flow features
            bld%Weff(ibuild)=2.*radius_out
            bld%Leff(ibuild)=2.*radius_out
! define region for non-local mixing
	 
! calculate building half widths	
! the width Wti is divided by a factor of 2 based on the Rockle 
! formula for multiple buildings
            bld%Wt(ibuild)=bld%Wti(ibuild)/2.
            bld%Lt(ibuild)=bld%Lti(ibuild)/2.
            bld%xfo(ibuild) = nint(xco - radius_out) !new bld%xfo of pentagon based on inner circle
            bld%istart(ibuild)=nint(bld%xfo(ibuild)/qugrid%dx) !front of the building	!var dz
            bld%iend(ibuild)= bld%istart(ibuild) + nint((bld%Lti(ibuild)-1)/qugrid%dx)	!back of the bld	!var dz
            bld%jend(ibuild)=nint(bld%yfo(ibuild)/qugrid%dy) + int((bld%Wt(ibuild)-1)/qugrid%dy) !far side of bld
            bld%jstart(ibuild)=nint(bld%yfo(ibuild)/qugrid%dy)-int(bld%Wt(ibuild)/qugrid%dy) !close side of bld
! end erp 1/05/05 change
! define length scales and retarding coefficients 1/25/05
            L_courtwake=sqrt(((0.75*point_coeff)**2.)+((0.85*flat_coeff)**2.))*bld%bb(ibuild)/2.5
            ret_coeff_court=sqrt(((0.6*point_coeff)**2.)+((0.35*flat_coeff)**2.))
            ret_coeff_courtwake=sqrt(((0.35*point_coeff)**2.)+((.3*flat_coeff)**2.))
            ret_coeff_wake=sqrt(((0.3*point_coeff)**2.)+((0.4*flat_coeff)**2.))
            ret_coeff_roof=0.3
            courtwake_exponent=sqrt(((0.2*point_coeff)**2.)+((0.3*flat_coeff)**2.))
            ew_denom_vert=sqrt(((0.4*point_coeff)**2.)+((0.5*flat_coeff)**2.))
            ew_denom_horz=sqrt(((0.15*point_coeff)**2.)+((0.1*flat_coeff)**2.))
            ew_start=sqrt(((0.6*radius_out*point_coeff)**2.)+((radius_in*flat_coeff)**2.))
            L_pentwake=bld%Ht(ibuild)*sqrt(((2.5*point_coeff)**2.)+((3*flat_coeff)**2.))
            farwake_exponent=sqrt(((0.70*point_coeff)**2.)+((0.58*flat_coeff)**2.))
            farwake_factor=sqrt(((15.0*point_coeff)**2.)+((26.3*flat_coeff)**2.))
            bld%Lr(ibuild)=L_pentwake
            if(flat_coeff .gt. 0.8)then
               bld%lf(ibuild)=2*bld%Ht(ibuild)
            endif
            if(abs(upwind_dir) .lt. pi/8.)then
               pentIstart=max(2,bld%istart(ibuild)-int(2*bld%Ht(ibuild)/qugrid%dy))
               pentIstop=qugrid%nx-1
               pentJstart=2
               pentJstop=qugrid%ny-1
            endif
            if(abs(upwind_dir) .gt. 7.*pi/8.)then
               pentIstart=2
               pentIstop=min(qugrid%nx-1,bld%iend(ibuild)+int(2*bld%Ht(ibuild)/qugrid%dx))
               pentJstart=2
               pentJstop=qugrid%ny-1
            endif
            if(upwind_dir .ge. pi/8. .and. upwind_dir .le. 3.*pi/8.)then
               pentIstart=max(2,bld%istart(ibuild)-int(2*bld%Ht(ibuild)/qugrid%dy))
               pentIstop=qugrid%nx-1
               pentJstart=max(2,bld%jstart(ibuild)-int(2*bld%Ht(ibuild)/qugrid%dy))
               pentJstop=qugrid%ny-1
            endif
            if(upwind_dir .le. -pi/8. .and. upwind_dir .ge. -3.*pi/8.)then
               pentIstart=max(2,bld%istart(ibuild)-int(2*bld%Ht(ibuild)/qugrid%dy))
               pentIstop=qugrid%nx-1
               pentJstart=2
               pentJstop=min(qugrid%ny-1,bld%jend(ibuild)+int(2*bld%Ht(ibuild)/qugrid%dy))
            endif
            if(upwind_dir .ge. 5.*pi/8. .and. upwind_dir .le. 7.*pi/8.)then
               pentIstart=2
               pentIstop=min(qugrid%nx-1,bld%iend(ibuild)+int(2*bld%Ht(ibuild)/qugrid%dx))
               pentJstart=max(2,bld%jstart(ibuild)-int(2*bld%Ht(ibuild)/qugrid%dy))
               pentJstop=qugrid%ny-1
            endif
            if(upwind_dir .le. -5.*pi/8. .and. upwind_dir .ge. -7.*pi/8.)then
               pentIstart=2
               pentIstop=min(qugrid%nx-1,bld%iend(ibuild)+int(2*bld%Ht(ibuild)/qugrid%dx))
               pentJstart=2
               pentJstop=min(qugrid%ny-1,bld%jend(ibuild)+int(2*bld%Ht(ibuild)/qugrid%dy))
            endif
            if(upwind_dir .gt. 3.*pi/8. .and. upwind_dir .lt. 5.*pi/8.)then
               pentIstart=2
               pentIstop=qugrid%nx-1
               pentJstart=max(2,bld%jstart(ibuild)-int(2*bld%Ht(ibuild)/qugrid%dy))
               pentJstop=qugrid%ny-1
            endif
            if(upwind_dir .lt. -3.*pi/8. .and. upwind_dir .gt. -5.*pi/8.)then
               pentIstart=2
               pentIstop=qugrid%nx-1
               pentJstart=2
               pentJstop=min(qugrid%ny-1,bld%jend(ibuild)+int(2*bld%Ht(ibuild)/qugrid%dy))
            endif
            do k=1,nz
               uz_ref(k)=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),k)
               vz_ref(k)=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),k)
            enddo
lp003:      do k=bld%kstart(ibuild),nz-1  !erp 7/23/03
               zb=zm(k)-bld%zfo(ibuild)		!var dz
lp002:         do j=pentJstart,pentJstop
lp001:            do i=pentIstart,pentIstop
                     x_u=((real(i)-1)*qugrid%dx-xco)*cos(upwind_dir)+		&	!var dz   &
                                  ((real(j)-0.5)*qugrid%dy-yco)*sin(upwind_dir)		!var dz
                     y_u=-((real(i)-1)*qugrid%dx-xco)*sin(upwind_dir)+		&	!var dz   &
                                  ((real(j)-0.5)*qugrid%dy-yco)*cos(upwind_dir)		!var dz
                     x_v=((real(i)-0.5)*qugrid%dx-xco)*cos(upwind_dir)+	&	!var dz   &
                                  ((real(j)-1)*qugrid%dy-yco)*sin(upwind_dir)			!var dz
                     y_v=-((real(i)-0.5)*qugrid%dx-xco)*sin(upwind_dir)+	&	!var dz   &
                                  ((real(j)-1)*qugrid%dy-yco)*cos(upwind_dir)			!var dz
!R_u and theta_u are polar coordinates of the cell in the rotated x_u, y_u coordinate system
                     R_u=sqrt(x_u**2+y_u**2)
                     R_v=sqrt(x_v**2+y_v**2)
                     theta_u=atan2(y_u,x_u)
                     theta_v=atan2(y_v,x_v)
! MAN 7/29/2005 add front recirculation for Pentagon building
                     if(bld%lf(ibuild) .gt. 0)then
! u values
                        if(yf2 .eq. yf1)then
                           xpf=xf1
                        else
                           xpf=((xf2-xf1)/(yf2-yf1))*(y_u-yf1)+xf1
                        endif
                        if(x_u .lt. xpf)then
                           if(y_u .ge. 0)then
                              Wp_eff=yf2
                           else
                              Wp_eff=yf1
                           endif
                           if(abs(y_u/Wp_eff) .le. 1.0 .and. zb .lt. 0.6*bld%Ht(ibuild))then
                              xf=-bld%lf(ibuild)*sqrt((1-(y_u/Wp_eff)**2)*(1-(zb/(0.6*bld%Ht(ibuild)))**2))
                              if(x_u .gt. xpf+xf)then
                                 quwinds%uo(i,j,k)=0.4*quwinds%uo(i,j,k)
                                 quwinds%wo(i,j,k)=0.
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=2
                                    endif
                                 endif
                                 if(zb .le. 0.5*bld%Ht(ibuild))then
                                    xf=0.4*xf
                                    if(x_u .gt. xpf+xf)then
                                       vel_mag=sqrt(uo_h**2.+vo_h**2.)
                                       quwinds%uo(i,j,k)=-uo_h*(-0.6*cos(((pi*zb)/(0.5*bld%Ht(ibuild))))+0.05)   &
                                                              *(-0.6*sin(((pi*abs(x_u-xpf))/(0.4*bld%lf(ibuild)))+0))
                                       quwinds%wo(i,j,k)=-vel_mag*(0.1*cos(((pi*abs(x_u-xpf))/(0.4*bld%lf(ibuild))))-0.05)
                                       if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                          if(icellflag(i,j,k) .ne. 0)then
                                             icellflag(i,j,k)=2
                                          endif
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                        endif
! v values
                        xpf=((xf2-xf1)/(yf2-yf1))*(y_v-yf1)+xf1
                        if(x_v .lt. xpf)then
                           if(y_v .ge. 0)then
                              Wp_eff=yf2
                           else
                              Wp_eff=yf1
                           endif
                           if(abs(y_v/Wp_eff) .le. 1.0 .and. zb .lt. 0.6*bld%Ht(ibuild))then
                              xf=-bld%lf(ibuild)*sqrt((1-(y_v/Wp_eff)**2)*(1-(zb/(0.6*bld%Ht(ibuild)))**2))
                              if(x_v .gt. xpf+xf)then
                                 quwinds%vo(i,j,k)=0.4*quwinds%vo(i,j,k)
                                 quwinds%wo(i,j,k)=0.
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=2
                                    endif
                                 endif
                                 if(zb .le. 0.5*bld%Ht(ibuild))then
                                    xf=0.4*xf
                                    if(x_v .gt. xpf+xf)then
                                       vel_mag=sqrt(uo_h**2.+vo_h**2.)
                                       quwinds%vo(i,j,k)=-vo_h*(-0.6*cos(((pi*zb)/(0.5*bld%Ht(ibuild))))+0.05)   &
                                                              *(-0.6*sin(((pi*abs(x_v-xpf))/(0.4*bld%lf(ibuild)))+0))
                                       quwinds%wo(i,j,k)=-vel_mag*(0.1*cos(((pi*abs(x_v-xpf))/(0.4*bld%lf(ibuild))))-0.05)
                                       if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                          if(icellflag(i,j,k) .ne. 0)then
                                             icellflag(i,j,k)=2
                                          endif
                                       endif
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
! end MAN 7/29/2005
! MAN 4/29/05 add far wake parameterization with building rotation
! u values
                     in_wake=0
                     if(y_u .ge. 0.)then
                        ynorm=yw4
                     else
                        ynorm=yw1
                     endif
                     if(y_u .le. yw2)then
                        xnorm=((xw2-xw1)/(yw2-yw1))*(y_u-yw1)+xw1
                     elseif(y_u .ge. yw3)then
                        xnorm=((xw3-xw4)/(yw3-yw4))*(y_u-yw4)+xw4
                     else
                        if(yw2 .eq. yw3)then
                           xnorm=xw2
                        else
                           xnorm=((xw2-xw3)/(yw2-yw3))*(y_u-yw2)+xw2
                        endif
                     endif
                     if(xw2 .eq. xw3)then
                        xnorm2=xw2
                     else
                        xnorm2=((xw2-xw3)/(yw2-yw3))*(-yw2)+xw2
                     endif
                     if(y_u .le. yw4 .and. y_u .ge. yw1 .and. k .le. bld%kend(ibuild))then
                        dNu=sqrt((1.-(y_u/ynorm)**2.)*   &
                                    (1.-((zb)/bld%Ht(ibuild))**2.)*L_pentwake**2)
                        if(x_u .gt. xnorm+dNu .and.   &
                                   theta_u .gt. thetaw1 .and. theta_u .lt. thetaw4)then
                           if(x_u-xnorm .le. farwake_factor*dNu)then
                              farwake_velocity=quwinds%uo(i,j,k)*(1.-(dNu/(x_u-xnorm+1.5*bld%Ht(ibuild)))**(farwake_exponent))
                              if(abs(farwake_velocity) .lt. abs(quwinds%uo(i,j,k)))then
                                 quwinds%uo(i,j,k)=farwake_velocity
                                 in_wake=1
                                 quwinds%wo(i,j,k)=0.
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=5
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
! Expanding Wake MAN 5/5/2005
                     if(in_wake .eq. 0)then
                        if(k .le. bld%kend(ibuild) .and. R_u .gt. radius_in)then
                           if(y_u .ge. 0.)then
                              xnorm3=xw4
                           else
                              xnorm3=xw1
                           endif
                           if(x_u .gt. xnorm3 .and. x_u .lt. farwake_factor*L_pentwake+xnorm2)then
                              if(x_u .gt. xnorm2)then
                                 Wp_eff=ynorm*sqrt(abs(1-((x_u-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                                 Ht_eff=bld%Ht(ibuild)*sqrt(abs(1-((x_u-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                              else
                                 Wp_eff=ynorm
                                 Ht_eff=bld%Ht(ibuild)
                              endif
                              u_ref=uz_ref(k)*(1.-(1/farwake_factor)**farwake_exponent)
                              zeta_ew=((y_u/(ew_denom_horz*2*Wp_eff)))*   &
                                                (((x_u+ew_start)/bld%Ht(ibuild))**(-0.5))!zb/(ew_denom_vert*Ht_eff)+
                              expanding_velocity=(u_ref-quwinds%uo(i,j,k))*exp(-(zeta_ew**2))+quwinds%uo(i,j,k)
                              if(abs(expanding_velocity) .le. abs(quwinds%uo(i,j,k)*0.98))then
                                 quwinds%uo(i,j,k)=expanding_velocity
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=7
                                    endif
                                 endif
                              endif
                           endif
                        elseif(k .gt. bld%kend(ibuild))then
                           if(x_u+ew_start .gt. 0. .and. x_u .lt. farwake_factor*L_pentwake+xnorm2)then
                              if(x_u .gt. xnorm2)then
                                 Wp_eff=ynorm*sqrt(abs(1-((x_u-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                                 Ht_eff=bld%Ht(ibuild)*sqrt(abs(1-((x_u-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                              else
                                 Wp_eff=ynorm
                                 Ht_eff=bld%Ht(ibuild)
                              endif
                              if(icellflag(i-1,j,bld%kend(ibuild)) .eq. 0 .or.   &
                                           icellflag(i,j,bld%kend(ibuild)) .eq. 0)then
                                 u_ref=uo_h*ret_coeff_roof
                              elseif(R_u .le. radius_out/2.5)then
                                 u_ref=quwinds%uo(i,j,bld%kend(ibuild))
                              else
                                 u_ref=uz_ref(bld%kend(ibuild))*(1.-(1/farwake_factor)**farwake_exponent)
                                 if(abs(u_ref) .gt. abs(quwinds%uo(i,j,bld%kend(ibuild))))then
                                    u_ref=quwinds%uo(i,j,bld%kend(ibuild))
                                 endif
                              endif
                              zeta_ew=(zb/(ew_denom_vert*Ht_eff)+   &
                                                (y_u/(ew_denom_horz*2*Wp_eff)))*   &
                                                (((x_u+ew_start)/bld%Ht(ibuild))**(-0.5))
                              expanding_velocity=(u_ref-quwinds%uo(i,j,k))*exp(-(zeta_ew**2))+quwinds%uo(i,j,k)
                              if(abs(expanding_velocity) .le. abs(quwinds%uo(i,j,k)*0.98))then
                                 quwinds%uo(i,j,k)=expanding_velocity
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=7
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
! end Expanding wake MAN 5/5/2005
! v values
                     in_wake=0
                     if(y_v .ge. 0.)then
                        ynorm=yw4
                     else
                        ynorm=yw1
                     endif
                     if(y_v .le. yw2)then
                        xnorm=((xw2-xw1)/(yw2-yw1))*(y_v-yw1)+xw1
                     elseif(y_u .ge. yw3)then
                        xnorm=((xw3-xw4)/(yw3-yw4))*(y_v-yw4)+xw4
                     else
                        if(yw2 .eq. yw3)then
                           xnorm=xw2
                        else
                           xnorm=((xw2-xw3)/(yw2-yw3))*(y_v-yw2)+xw2
                        endif
                     endif
                     if(xw2 .eq. xw3)then
                        xnorm2=xw2
                     else
                        xnorm2=((xw2-xw3)/(yw2-yw3))*(-yw2)+xw2
                     endif
                     if(y_v .le. yw4 .and. y_v .ge. yw1 .and. k .le. bld%kend(ibuild))then
                        dNv=sqrt((1.-(y_v/ynorm)**2.)*(1.-((zb)/bld%Ht(ibuild))**2.)*L_pentwake**2)
                        if(x_v .gt. xnorm+dNv .and. theta_v .gt. thetaw1 &
                                 .and. theta_v .lt. thetaw4)then
                           if(x_v-xnorm .le. farwake_factor*dNv)then
                              farwake_velocity=quwinds%vo(i,j,k)*(1.-(dNv/(x_v-xnorm+1.5*bld%Ht(ibuild)))**(farwake_exponent))
                              if(abs(farwake_velocity) .lt. abs(quwinds%vo(i,j,k)))then
                                 quwinds%vo(i,j,k)=farwake_velocity
                                 quwinds%wo(i,j,k)=0.
                                 in_wake=1
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=5
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
! Expanding Wake MAN 5/5/2005
                     if(in_wake .eq. 0)then
                        if(k .le. bld%kend(ibuild) .and. R_v .gt. radius_in)then
                           if(y_v .ge. 0.)then
                              xnorm3=xw4
                           else
                              xnorm3=xw1
                           endif
                           if(x_v .gt. xnorm3 .and. x_v .lt. farwake_factor*L_pentwake+xnorm2)then
                              if(x_v .gt. xnorm2)then
                                 Wp_eff=ynorm*sqrt(abs(1-((x_v-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                                 Ht_eff=bld%Ht(ibuild)*sqrt(abs(1-((x_v-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                              else
                                 Wp_eff=ynorm
                                 Ht_eff=bld%Ht(ibuild)
                              endif
                              v_ref=vz_ref(k)*(1.-(1/farwake_factor)**farwake_exponent)
                              zeta_ew=((y_v/(ew_denom_horz*2*Wp_eff)))*   &
                                                (((x_v+ew_start)/bld%Ht(ibuild))**(-0.5))
                              expanding_velocity=(v_ref-quwinds%vo(i,j,k))*exp(-(zeta_ew**2))+quwinds%vo(i,j,k)
                              if(abs(expanding_velocity) .le. abs(quwinds%vo(i,j,k)*0.98))then
                                 quwinds%vo(i,j,k)=expanding_velocity
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=7
                                    endif
                                 endif
                              endif
                           endif
                        elseif(k .gt. bld%kend(ibuild))then
                           if(x_v+ew_start .gt. 0. .and. x_v .lt. farwake_factor*L_pentwake+xnorm2)then
                              if(x_v .gt. xnorm2)then
                                 Wp_eff=ynorm*sqrt(abs(1-((x_v-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                                 Ht_eff=bld%Ht(ibuild)*sqrt(abs(1-((x_v-xnorm2)/   &
                                                   (farwake_factor*L_pentwake))**2.))
                              else
                                 Wp_eff=ynorm
                                 Ht_eff=bld%Ht(ibuild)
                              endif
                              if(icellflag(i-1,j,bld%kend(ibuild)) .eq. 0 .or.   &
                                           icellflag(i,j,bld%kend(ibuild)) .eq. 0)then
                                 v_ref=vo_h*ret_coeff_roof
                              elseif(R_v .le. radius_out/2.5)then
                                 v_ref=quwinds%vo(i,j,bld%kend(ibuild))
                              else
                                 v_ref=vz_ref(bld%kend(ibuild))*   &
                                                  (1.-(1/farwake_factor)**farwake_exponent)
                                 if(abs(v_ref) .gt. abs(quwinds%vo(i,j,bld%kend(ibuild))))then
                                    v_ref=quwinds%vo(i,j,bld%kend(ibuild))
                                 endif
                              endif
                              zeta_ew=(zb/(ew_denom_vert*Ht_eff)+   &
                                                (y_v/(ew_denom_horz*2*Wp_eff)))*   &
                                                (((x_v+ew_start)/bld%Ht(ibuild))**(-0.5))
                              expanding_velocity=(v_ref-quwinds%vo(i,j,k))*exp(-(zeta_ew**2))+quwinds%vo(i,j,k)
                              if(abs(expanding_velocity) .le. abs(quwinds%vo(i,j,k)*0.98))then
                                 quwinds%vo(i,j,k)=expanding_velocity
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=7
                                    endif
                                 endif
                              endif
                           endif
                        endif
                     endif
! end Expanding wake MAN 5/5/2005
! end MAN 5/2/05 far wake parameterization with building rotation
! Near wake parameterization with building rotation MAN 4/29/05
                     if(y_u .le. yw4 .and. y_u .ge. yw1 .and. k .le. bld%kend(ibuild))then
                        if(R_u .gt. radius_in .and.   &
                                   theta_u .gt. thetaw1 .and. theta_u .lt. thetaw4)then
                           if(x_u .le. dNu+xnorm .and. x_u .ge. xnorm)then
                              quwinds%wo(i,j,k)=0.1*sqrt(uz_ref(k)**2.+uz_ref(k)**2.)*(1-abs((x_u-xnorm)/dNu))
                              if(y_u .le. yw2)then
                                 vel_mag=sqrt(uo_h**2.+vo_h**2.)
                                 quwinds%uo(i,j,k)=(0.3+0.7*y_u/ynorm)*(ret_coeff_wake*vel_mag*   &
                                                      cos(normw1-pi/6.+0.67*theta_u)*(1-(zb/bld%Ht(ibuild))**2.)+   &
                                                      quwinds%uo(i,j,k)*((zb/bld%Ht(ibuild))**2.))
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=4
                                    endif
                                 endif
                              elseif(y_u .ge. yw3)then
                                 vel_mag=sqrt(uo_h**2.+vo_h**2.)
                                 quwinds%uo(i,j,k)=(0.3+0.7*y_u/ynorm)*(ret_coeff_wake*vel_mag*   &
                                                      cos(normw3+pi/6.+0.67*theta_u)*(1-(zb/bld%Ht(ibuild))**2.)+   &
                                                      quwinds%uo(i,j,k)*((zb/bld%Ht(ibuild))**2.))
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=4
                                    endif
                                 endif
                              else
                                 vel_mag=sqrt(uo_h**2.+vo_h**2.)*(1-((x_u-xnorm)/dNu)**2.)
                                 Ur_u=ret_coeff_wake*vel_mag*(-0.8-cos(theta_u*180./36.))   &
                                                  *(1.3-cos(theta_u*90./36.))
                                 Vr_u=-1.3*ret_coeff_wake*vel_mag*(sin(3.*theta_u)**3.0)
                                 quwinds%uo(i,j,k)=(Ur_u*cos(-upwind_dir)+Vr_u*sin(-upwind_dir))*   &
                                                       (1-(zb/bld%Ht(ibuild))**2.)+(0.25+0.75*y_u/ynorm)*((zb/bld%Ht(ibuild))**2.)   &
                                                       *quwinds%uo(i,j,k)*(1-(dNu/(dNu+1.5*bld%Ht(ibuild)))**farwake_exponent)
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                              !quwinds%wo(i,j,k)=0.
                              if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                 write(325,*)'Parameterized U exceeds max in pentagon wake',&
                                    quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                              endif
                           endif
                        endif
                     endif
                     if(y_v .le. yw4 .and. y_v .ge. yw1 .and. k .le. bld%kend(ibuild))then
                        if(R_v .gt. radius_in .and.   &
                                   theta_v .gt. thetaw1 .and. theta_v .lt. thetaw4)then
                           if(x_v .le. dNv+xnorm .and. x_v .ge. xnorm)then
                              quwinds%wo(i,j,k)=0.1*sqrt(uz_ref(k)**2.+uz_ref(k)**2.)*(1-abs((x_v-xnorm)/dNv))
                              if(y_v .le. yw2)then
                                 vel_mag=sqrt(uo_h**2.+vo_h**2.)
                                 quwinds%vo(i,j,k)=(0.3+0.7*y_v/ynorm)*(ret_coeff_wake*vel_mag*   &
                                                      sin(normw1-pi/6.+0.67*theta_v)*(1-(zb/bld%Ht(ibuild))**2.)+   &
                                                      quwinds%vo(i,j,k)*(zb/bld%Ht(ibuild))**2.)
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=4
                                    endif
                                 endif
                              elseif(y_v .ge. yw3)then
                                 vel_mag=sqrt(uo_h**2.+vo_h**2.)
                                 quwinds%vo(i,j,k)=(0.3+0.7*y_v/ynorm)*(ret_coeff_wake*vel_mag*   &
                                                      sin(normw3+pi/6.+0.67*theta_v)*(1-(zb/bld%Ht(ibuild))**2.)+   &
                                                      quwinds%vo(i,j,k)*(zb/bld%Ht(ibuild))**2.)
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=4
                                    endif
                                 endif
                              else
                                 vel_mag=sqrt(uo_h**2.+vo_h**2.)*(1-((x_v-xnorm)/dNv)**2.)
                                 Ur_v=ret_coeff_wake*vel_mag*(-0.8-cos(theta_v*180./36.))   &
                                                  *(1.3-cos(theta_v*90./36.))
                                 Vr_v=-1.3*ret_coeff_wake*vel_mag*(sin(3.*theta_v)**3.0)
                                 quwinds%vo(i,j,k)=(-Ur_v*sin(-upwind_dir)+Vr_v*cos(-upwind_dir))*   &
                                                       (1-(zb/bld%Ht(ibuild))**2.)+(0.25+0.75*y_v/ynorm)*((zb/bld%Ht(ibuild))**2.)   &
                                                       *quwinds%vo(i,j,k)*(1-(dNv/(dNv+1.5*bld%Ht(ibuild)))**farwake_exponent)
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                    if(icellflag(i,j,k) .ne. 0)then
                                       icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                              !quwinds%wo(i,j,k)=0.
                              if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                 write(325,*)'Parameterized V exceeds max in pentagon wake',&
                                    quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                              endif
                           endif
                        endif
                     endif
! end MAN 5/2/05 near wake parameterization with building rotation
! make entire solid pentagon with rotation angle bld%gamma MAN 4/25/05
                     if(j .ge. bld%jstart(ibuild) .and. j .le. bld%jend(ibuild)   &
                               .and. i .ge. bld%istart(ibuild) .and. i .le. bld%iend(ibuild)   &
                               .and. k .le. bld%kend(ibuild))then
                        xpent=((real(i)-0.5)*qugrid%dx - xco)*cos(bld%gamma(ibuild))+	&	!var dz   
                                      ((real(j)-0.5)*qugrid%dy - yco)*sin(bld%gamma(ibuild))		!var dz
                        ypent=-((real(i)-0.5)*qugrid%dx -xco)*sin(bld%gamma(ibuild))+	&	!var dz   
                                      ((real(j)-0.5)*qugrid%dy-yco)*cos(bld%gamma(ibuild))		!var dz
                        xpentu=((real(i)-1)*qugrid%dx - xco)*cos(bld%gamma(ibuild))+	&		!var dz   
                                       ((real(j)-0.5)*qugrid%dy - yco)*sin(bld%gamma(ibuild))		!var dz
                        ypentu=-((real(i)-1)*qugrid%dx - xco)*sin(bld%gamma(ibuild))+	&		!var dz   
                                       ((real(j)-0.5)*qugrid%dy - yco)*cos(bld%gamma(ibuild))		!var dz
                        xpentv=((real(i)-0.5)*qugrid%dx - xco)*cos(bld%gamma(ibuild))+	&	!var dz   
                                       ((real(j)-1)*qugrid%dy - yco)*sin(bld%gamma(ibuild))			!var dz
                        ypentv=-((real(i)-0.5)*qugrid%dx - xco)*sin(bld%gamma(ibuild))+ &	!var dz   
                                       ((real(j)-1)*qugrid%dy - yco)*cos(bld%gamma(ibuild))			!var dz
                        if(ypent .gt. y2 .and. ypent .lt. y1)then
                           if((xpent .lt. (x1-x2)/(y1-y2)*(ypent - y2)+x2) .and.   &
                                      (xpent .gt. (x4-x3)/(y4-y3)*(ypent - y3)+x3))then
                              icellflag(i,j,k)=0
                              ! ibldflag(i,j,k)=ibuild ! MAN 8/29/2007 building flags
                           endif
                        endif
                        if(ypent .gt. y1 .and. ypent .lt. y5)then
                           if((xpent .lt. (x5-x1)/(y5-y1)*(ypent - y1)+x1) .and.   &
                                      (xpent .gt. (x5-x4)/(y5-y4)*(ypent - y4)+x4))then
                              icellflag(i,j,k)=0
                              ! ibldflag(i,j,k)=ibuild ! MAN 8/29/2007 building flags
                           endif
                        endif
! begin new pentagon with central hole
! set up the wind parameterization in the courtyard MAN 4/25/05
                        in_court=0
                        if(ypent .gt. y2in .and. ypent .lt. y1in)then
                           if((xpent .lt. (x1in-x2in)/(y1in-y2in)*(ypent - y2in)+x2in) .and.   &
                                      (xpent .gt. (x4in-x3in)/(y4in-y3in)*(ypent - y3in)+x3in))then
                              icellflag(i,j,k)=1
                              ! ibldflag(i,j,k)=0 ! MAN 8/29/2007 building flags
                              in_court=1
                           endif
                        endif
                        if(ypent .gt. y1in .and. ypent .lt. y5in)then
                           if((xpent .lt. (x5in-x1in)/(y5in-y1in)*(ypent - y1in)+x1in) .and.   &
                                      (xpent .gt. (x5in-x4in)/(y5in-y4in)*(ypent - y4in)+x4in))then
                              icellflag(i,j,k)=1
                              ! ibldflag(i,j,k)=0 ! MAN 8/29/2007 building flags
                              in_court=1
                           endif
                        endif
                        if(in_court .eq. 1.)then
                           in_courtwake_u=0
                           in_courtwake_v=0
                           x_nearu=max(real(qugrid%nx)*qugrid%dx,real(qugrid%ny)*qugrid%dy)
                           x_nearv=x_nearu
                           if(cos(upwind_dir-norm1) .gt. 0.01)then
                              delta_y=ypentu-((y1in-y5in)/(x1in-x5in)*(xpentu - x5in)+y5in)
                              x_eff=abs(delta_y*sin(norm1-bld%gamma(ibuild))/cos(upwind_dir-norm1))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_u=1
                              endif
                              if(x_eff .lt. x_nearu)x_nearu=x_eff
                              delta_y=ypentv-((y1in-y5in)/(x1in-x5in)*(xpentv - x5in)+y5in)
                              x_eff=abs(delta_y*sin(norm1-bld%gamma(ibuild))/cos(upwind_dir-norm1))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_v=1
                              endif
                              if(x_eff .lt. x_nearv)x_nearv=x_eff
                           endif
                           if(cos(upwind_dir-norm2) .gt. 0.01)then
                              delta_y=ypentu-((y1in-y2in)/(x1in-x2in)*(xpentu - x2in)+y2in)
                              x_eff=abs(delta_y*sin(norm2-bld%gamma(ibuild))/cos(upwind_dir-norm2))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_u=1
                              endif
                              if(x_eff .lt. x_nearu)x_nearu=x_eff
                              delta_y=ypentv-((y1in-y2in)/(x1in-x2in)*(xpentv - x2in)+y2in)
                              x_eff=abs(delta_y*sin(norm2-bld%gamma(ibuild))/cos(upwind_dir-norm2))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_v=1
                              endif
                              if(x_eff .lt. x_nearv)x_nearv=x_eff
                           endif
                           if(cos(upwind_dir-norm3) .gt. 0.01)then
                              x_eff=(ypentu-y3in)/cos(upwind_dir-norm3)
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_u=1
                              endif
                              if(x_eff .lt. x_nearu)x_nearu=x_eff
                              x_eff=(ypentv-y3in)/cos(upwind_dir-norm3)
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_v=1
                              endif
                              if(x_eff .lt. x_nearv)x_nearv=x_eff
                           endif
                           if(cos(upwind_dir-norm4) .gt. 0.01)then
                              delta_y=ypentu-((y4in-y3in)/(x4in-x3in)*(xpentu - x4in)+y4in)
                              x_eff=abs(delta_y*sin(norm4-bld%gamma(ibuild))/cos(upwind_dir-norm4))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_u=1
                              endif
                              if(x_eff .lt. x_nearu)x_nearu=x_eff
                              delta_y=ypentv-((y4in-y3in)/(x4in-x3in)*(xpentv - x4in)+y4in)
                              x_eff=abs(delta_y*sin(norm4-bld%gamma(ibuild))/cos(upwind_dir-norm4))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_v=1
                              endif
                              if(x_eff .lt. x_nearv)x_nearv=x_eff
                           endif
                           if(cos(upwind_dir-norm5) .gt. 0.01)then
                              delta_y=ypentu-((y5in-y4in)/(x5in-x4in)*(xpentu - x4in)+y4in)
                              x_eff=abs(delta_y*sin(norm5-bld%gamma(ibuild))/cos(upwind_dir-norm5))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_u=1
                              endif
                              if(x_eff .lt. x_nearu)x_nearu=x_eff
                              delta_y=ypentv-((y5in-y4in)/(x5in-x4in)*(xpentv - x4in)+y4in)
                              x_eff=abs(delta_y*sin(norm5-bld%gamma(ibuild))/cos(upwind_dir-norm5))
                              if(zb .lt. bld%Ht(ibuild)*(1.-((x_eff/L_courtwake)**2.)))then
                                 in_courtwake_v=1
                              endif
                              if(x_eff .lt. x_nearv)x_nearv=x_eff
                           endif
                           if(in_courtwake_u .eq. 1)then
                              quwinds%uo(i,j,k)=-(uz_ref(k)*(1-abs(x_nearu/L_courtwake)**0.5)+uo_h*(abs(x_nearu/L_courtwake)**0.5))*   &
                                                  ret_coeff_courtwake*(1-(zb/bld%Ht(ibuild))**0.5)
                              quwinds%wo(i,j,k)=0.1*sqrt(uz_ref(k)**2+vz_ref(k)**2)*(1-(x_nearu/L_courtwake)**2)
                              if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                 if(icellflag(i,j,k) .ne. 0)then
                                    icellflag(i,j,k)=4
                                 endif
                              endif
                           else
                              if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                 if(icellflag(i,j,k) .ne. 0)then
                                    if(zb .ge. 0.8*(bld%Ht(ibuild)-bld%zfo(ibuild)) .or. &
                                          (icellflag(i,j,2) .eq. 4 .and. k .gt. 2))then
                                       icellflag(i,j,k)=9
                                    else
                                       icellflag(i,j,k)=5
                                    endif
                                 endif
                              endif
                              quwinds%uo(i,j,k)=ret_coeff_court*uo_h*(((zb)/bld%Ht(ibuild))**courtwake_exponent)	!var dz
                           endif
                           if(in_courtwake_v .eq. 1)then
                              quwinds%vo(i,j,k)=-(vz_ref(k)*(1-abs(x_nearv/L_courtwake)**0.5)+vo_h*(abs(x_nearv/L_courtwake)**0.5))*   &
                                                  ret_coeff_courtwake*(1-(zb/bld%Ht(ibuild))**0.5)
                              quwinds%wo(i,j,k)=0.1*sqrt(uz_ref(k)**2+vz_ref(k)**2)*(1-(x_nearv/L_courtwake)**2)
                              if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                 if(icellflag(i,j,k) .ne. 0)then
                                    icellflag(i,j,k)=4
                                 endif
                              endif
                           else
                              if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. nz)then
                                 if(icellflag(i,j,k) .ne. 0)then
                                    if(zb .ge. 0.8*(bld%Ht(ibuild)-bld%zfo(ibuild)) .or. &
                                          (icellflag(i,j,2) .eq. 4 .and. k .gt. 2))then
                                       icellflag(i,j,k)=9
                                    else
                                       icellflag(i,j,k)=5
                                    endif
                                 endif
                              endif
                              quwinds%vo(i,j,k)=ret_coeff_court*vo_h*(((zb)/bld%Ht(ibuild))**courtwake_exponent)		!var dz
                           endif
                        endif
                     endif
!end courtyard parameterization man 1/12/05
!end pentagon with central hole erp 11/22/04
                  enddo   lp001      
               enddo   lp002      
            enddo   lp003      
            deallocate(uz_ref,vz_ref)
         endif
!end routine to generate Pentagon shaped building
         return
      end

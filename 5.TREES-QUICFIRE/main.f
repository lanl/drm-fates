!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! program to populate a grid with an established fuel map
!
! Author: Alexander Josephson (11/19)
! Last Modified: 11/19
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program fuel_maps

      print *,'===================================='
      print *,' Running TREES to generate fuel     '
      print *,' files for FIRETEC or QUIC-Fire     '
      print *,'===================================='

      !-----Initialize
      call namelist_input
      call define_constant_variables
      call define_grid_variables
      print*, 'HERE ADAM 0' 
      !-----Establish baseline
      call baseline
      print*, 'HERE ADAM 1'
      !-----Perform fuel treatments
      call treatment
      print*, 'HERE ADAM 2'
      !-----Export data to binary files
      call output

      end

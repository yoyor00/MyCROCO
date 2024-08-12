!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
! A long character string to hold activated cpp-switches.
! Basically it is used to keep track of cpp-switches by placing
! them together and writing into history file.
!                                    !
      integer max_opt_size           ! NOTE: Parameter max_opt_size
      parameter (max_opt_size=4500)  ! must be equal to the length
      character*4500 Coptions,srcs   ! of character string.
      common /strings/ Coptions,srcs !

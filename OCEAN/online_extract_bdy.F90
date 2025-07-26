!======================================================================
      ! CROCO is a branch of ROMS developped at IRD, INRIA, 
      ! Ifremer, CNRS and Univ. Toulouse III  in France
      ! The two other branches from UCLA (Shchepetkin et al)
      ! and Rutgers University (Arango et al) are under MIT/X style
      ! license.
      ! CROCO specific routines (nesting) are under CeCILL-C license.
      !
      ! CROCO website : http://www.croco-ocean.org
      !======================================================================
      !
      module online_extract_bdy
!        
#include "cppdefs.h"

#ifdef ONLINE_EXTRACT
      ! Extract data at various positions and frequencies
      ! Only for things that need interpolating


      ! Auxiliary tools required:
      ! - Tools-Roms/scripts/pre/add_object.m -> add a data extraction object to a netcdf file
      ! - Tools-Roms/scripts/pre/bry_extract_obj.m -> Example showing
      !             how to add objects for a child grid boundary forcing file

      !
      ! STEPS:
      !
      ! 0) CREATE INPUT FILE OF EXTRACTION OBJECTS
      !
      !    - Positions (in fractional i,j locations of the grid)
      !    - Name of the object ('child1_south','mooring2', etc, etc,...
      !    - Output frequency of each object
      !    - Variables to be output
      !  
      !
      ! 1) Assign positions to subdomains
      !    - if ispresent(iobj) then see above
      !    - keep track of placement in global object arrays
      !    The ROMS parent simulation needs to know where the child boundary sits in its
      !    parent domain. To do this the child bry points are given i and j coords in
      !    terms of the parent grid. We created a matlab tool to do this found here:
      !    Tools-Roms/scripts/pre/

      ! 2) Create a single file per subdomain containing all the objects
      !    for which ispresent(iobj) is true
      !
      ! 3) Loop through all objects and all vars, and write when needed
      !    Add averaging capability at some point
      !
      !    Vel points always need to interpolate both u, and v in order
      !    to rotate the vector to a desired angle
      !
      !]

      use netcdf
      use module_extract
!
      implicit none
!
!!!# if defined MPI & !defined PARALLEL_FILES
# if defined MPI
      include 'mpif.h'
      integer :: status(MPI_STATUS_SIZE), blank
# endif     
# ifdef MPI
#  define LOCALLM Lmmpi
#  define LOCALMM Mmmpi
# else
#  define LOCALLM Lm
#  define LOCALMM Mm
# endif
# ifdef OUT_DOUBLE
#  define NF90_OUT NF90_DOUBLE
#  define NF90_FILL NF90_FILL_DOUBLE
# else
#  define NF90_OUT NF90_REAL
#  define NF90_FILL NF90_FILL_REAL
# endif      

      ! TODO: add averaging


      private

      integer :: total_rec=0 ! total records
      integer :: nobj  ! number of extraction objects

      logical                         :: extend_up !! flag to extend up,vp
      real,dimension(:,:),allocatable :: upe,vpe  !! buffer filled versions of up,vp

      integer, dimension(:), allocatable :: ncidextr !netcdf ids for all objects
      character(len=20), dimension(:), allocatable :: unique_set !name of each set
      integer, dimension(:), allocatable :: extrTime, &!variable ids for all objects
                                            extrZ_N , extrZ_S , extrZ_W , extrZ_E , & 
                                            extrUb_N, extrUb_S, extrUb_W, extrUb_E, &
                                            extrVb_N, extrVb_S, extrVb_W, extrVb_E, &
                                            extrU_N , extrU_S , extrU_W , extrU_E , &
                                            extrV_N , extrV_S , extrV_W , extrV_E , &
                                            extrT_N , extrT_S , extrT_W , extrT_E , &
                                            extrS_N , extrS_S , extrS_W , extrS_E 
      real, dimension(:,:), allocatable :: qext
      integer, save :: allc_ext_size=0, t_count=0, mskd_pts=0
      integer, save :: alloc_msk_size=0
      real, dimension(:), allocatable :: ptch
      integer, dimension(:,:), allocatable :: ijmsk
      integer, dimension(:,:), allocatable :: ijptch
      integer, dimension(:), allocatable :: nijmsk

      type extract_object  ! contains all information for a data_extraction object

        ! needed as input
        character(len=60) :: obj_name                   ! name of object
        character(len=20) :: set                        ! name of set that the object belongs to
        character(len=20) :: bnd                        ! name of boundary (for bry type data)
        character(len=20) :: vargrid                      ! grid of variable to output (u,v or rho)
        character(:),allocatable :: pre                 ! preamble for vars and dims
        character(:),allocatable :: dname               ! dimension name for object
        integer                  :: dsize               ! dimension of object
        logical                            :: scalar    ! scalar or vector
        !real                               :: period    ! output period (seconds)
        integer                            :: np        ! local number of locations
        real,   dimension(:)  ,allocatable :: ipos,jpos ! fractional index locations
        real,   dimension(:)  ,allocatable :: ang       ! desired angle for velocities
        integer                            :: start_idx ! starting position in the obj array


        !! Initializing record at nrpf will trigger the making of a file
        !! at the first time step
        !integer                            :: record    ! record number in file
        real                               :: otime=0   ! time since last output

        real,   dimension(:,:),pointer     :: vari      ! data for output

        ! these are only for scalars
        integer,dimension(:)  ,pointer     :: ip,jp     ! lower left hand point index
        real   ,dimension(:,:),pointer     :: coef      ! interpolation coefficients

        ! these are only for velocities
        real,   dimension(:)  ,allocatable :: cosa,sina ! for rotation of vectors
        integer,dimension(:)  ,pointer     :: ipu,jpu   ! only for vectors
        integer,dimension(:)  ,pointer     :: ipv,jpv   ! only for vectors
        real   ,dimension(:,:),pointer     :: cfu,cfv   ! only for vectors
        real   ,dimension(:,:),pointer     :: ui,vi     ! only for vectors 

        ! These logicals determine which variables are desired for an
        ! object. Only the ones listed below are currently functional
        ! False by default
        logical :: zeta = .false.
        logical :: ubar = .false.
        logical :: vbar = .false.
        logical :: u    = .false.
        logical :: v    = .false.
        logical :: temp = .false.
        logical :: salt = .false.
        logical :: up   = .false.
        logical :: vp   = .false.

      end type extract_object

      type(extract_object),dimension(:),allocatable :: obj

      interface interpolate
        module procedure  interpolate_2D, interpolate_3D
      end interface

      public :: do_extract_bdy, dealloc_closecdf 
      contains
! ----------------------------------------------------------------------
      subroutine init_extract_data  ![
      ! Allocate space and compute interpolation coefficients
      ! for rho,u,and v variables.
      implicit none

      !local
      character(len=30) :: preamb
      integer :: i,np,ierr,lpre
      real,dimension(:),allocatable :: angp  ! grid angle

      call read_extraction_objects
!!!      MPI_master_only write(stdout,*) 'unique_set',unique_set 

      !!!if (diag_pflx) then
      !!!  allocate( upe(GLOBAL_2D_ARRAY) ); upe = 0
      !!!  allocate( vpe(GLOBAL_2D_ARRAY) ); vpe = 0
      !!!endif

      do i = 1,nobj
        np = obj(i)%np
        if (np>0) then
!!!          write(stdout,'(A,2x,I2,2x,I2,2x,I2,A,I2)') 
!!!     &      'read_extr',i,mynode,np,'/',obj(i)%dsize

          preamb = trim(obj(i)%obj_name)
          lpre = len(trim(preamb))-1
          allocate(character(len=lpre) :: obj(i)%pre)
          obj(i)%pre = preamb(1:lpre)

          allocate(obj(i)%vari(np,N))
          allocate(obj(i)%coef(np,4))
          allocate(obj(i)%ip(np))                          ! DON'T NEED THESE FOR U or V objects
          allocate(obj(i)%jp(np))

          ! from absolute index to rho-index
          ! ipos,jpos are in 'absolute' index space: [0,nx]x[0,ny]
          obj(i)%ipos = obj(i)%ipos+0.5
          obj(i)%jpos = obj(i)%jpos+0.5
          call compute_coef(obj(i)%ipos, obj(i)%jpos, obj(i)%coef, obj(i)%ip, obj(i)%jp &
# ifdef MASKING            
            , rmask &
# endif
                           )

          if (.not.obj(i)%scalar) then
            allocate(obj(i)%cosa(np))
            allocate(obj(i)%sina(np))

            allocate(angp(np))
            call interpolate(angler, angp, obj(i)%coef, obj(i)%ip, obj(i)%jp)
            obj(i)%cosa = cos(angp-obj(i)%ang)
            obj(i)%sina = sin(angp-obj(i)%ang)

            deallocate(angp)

            allocate(obj(i)%ui(np,N))
            allocate(obj(i)%cfu(np,4))
            allocate(obj(i)%ipu(np))
            allocate(obj(i)%jpu(np))

            ! from rho-index to u-index 
            obj(i)%ipos = obj(i)%ipos+0.5
            call compute_coef(obj(i)%ipos, obj(i)%jpos, obj(i)%cfu, obj(i)%ipu, obj(i)%jpu &
# ifdef MASKING            
            , umask &
# endif
                             )
            allocate(obj(i)%vi(obj(i)%np,N))
            allocate(obj(i)%cfv(obj(i)%np,4))
            allocate(obj(i)%ipv(obj(i)%np))
            allocate(obj(i)%jpv(obj(i)%np))

            ! from u-index to v-index
            obj(i)%ipos = obj(i)%ipos-0.5
            obj(i)%jpos = obj(i)%jpos+0.5
            call compute_coef(obj(i)%ipos, obj(i)%jpos, obj(i)%cfv, obj(i)%ipv, obj(i)%jpv &
# ifdef MASKING            
            , vmask &
# endif
                             )
          endif                   
# ifdef NC4PAR
        else
            allocate(obj(i)%vari(np,N))
# endif
        endif
      enddo

      end subroutine init_extract_data !]
! ----------------------------------------------------------------------
      subroutine read_extraction_objects  ![
      ! Read all objects from file and determine local ranges 
      implicit none

      ! local
      integer               :: iobj,ncid,varid
      integer,dimension(2)  :: dimids
      character(len=20) :: dname
      character(len=30) :: objname
      character(len=80) :: variables
      real,dimension(:,:),allocatable :: object
      logical,dimension(:), allocatable :: setmask
      integer               :: n1,n2,i0,i1,lstr,lenstr
      integer ierr,sidx,sidx2,num
      !real :: period

      MPI_master_only write(stdout,'(7x,2A)') 'extract_data :: read objects: ',inp_extrname

      ierr = nf90_open(inp_extrname,nf90_nowrite,ncid)
      ierr = nf90_inquire(ncid,nVariables=nobj)

      allocate(obj(nobj))

      ! Read all objects from file. 
      do iobj = 1,nobj
        ierr = nf90_inquire_variable(ncid,iobj,name=objname,dimids=dimids)
        ierr = nf90_inquire_dimension(ncid,dimids(1),name=dname,len=n1)
        ierr = nf90_inquire_dimension(ncid,dimids(2),len=n2)

        lstr = len(trim(dname))
        allocate(character(len=lstr) :: obj(iobj)%dname)
        obj(iobj)%dname = trim(dname)
        obj(iobj)%dsize = n1
        
        ! scalar objects have i,j. Vector obj also have ang
        if (n2==2) obj(iobj)%scalar = .true.

        allocate(object(n1,n2))
        ierr = nf90_inq_varid(ncid,objname,varid)
        ierr = nf90_get_var(ncid,varid,object)

        call find_local_points(obj(iobj),object(:,1),object(:,2))

        obj(iobj)%obj_name = objname
        !lstr = len(trim(objname))
        lstr = lenstr(objname)
        if (findstr(objname,'_',sidx) ) then
          obj(iobj)%set = objname(1:sidx-1)
          if (findstr(objname(sidx+1:lstr),'_',sidx2) ) then
            obj(iobj)%bnd = objname(sidx+1:sidx+sidx2-1)
          else
            obj(iobj)%bnd = ' '
          endif
        else
          obj(iobj)%set = objname
        endif
!!!        write(stdout,*) 'set:',obj(iobj)%set,lenstr(obj(iobj)%set)

!       if (mynode==84) then
!         print *,i,obj(iobj)%obj_name,obj(iobj)%set,obj(iobj)%bnd
!         print *,trim(obj(iobj)%set)//'_'//vname//'_'//trim(obj(iobj)%bnd)
!         stop
!       endif

!        ierr = nf_get_att(ncid,iobj,'output_period',period)
!        obj(iobj)%period   = period

        ! Figure out which variables to output for this object
        ierr = nf90_get_att(ncid,iobj,'output_vars',variables)

        if (findstr(variables,'zeta') ) obj(iobj)%zeta = .True.
        if (findstr(variables,'temp') ) obj(iobj)%temp = .True.
        if (findstr(variables,'salt') ) obj(iobj)%salt = .True.
        if (findstr(variables,'ubar') ) obj(iobj)%ubar = .True.
        if (findstr(variables,'vbar') ) obj(iobj)%vbar = .True.
        if (findstr(variables,'u'   ) ) obj(iobj)%u    = .True.
        if (findstr(variables,'v'   ) ) obj(iobj)%v    = .True.

!!!      MPI_master_only write(stdout,*) 'findstr:',objname,obj(iobj)%zeta,
!!!     &     obj(iobj)%temp,obj(iobj)%salt,obj(iobj)%ubar,obj(iobj)%vbar,
!!!     &     obj(iobj)%u,obj(iobj)%v   

        if (findstr(variables,'zeta') .or. &
            findstr(variables,'temp') .or. &
            findstr(variables,'salt')) obj(iobj)%vargrid = 'rho'
        if (findstr(variables,'u') .or. &
            findstr(variables,'ubar')) obj(iobj)%vargrid = 'u'
        if (findstr(variables,'v') .or. &
            findstr(variables,'vbar')) obj(iobj)%vargrid = 'v'

!!!        MPI_master_only write(stdout,'(A,2x,A,2x,A,2x,A,2x,A)') 
!!!     &          'read_extr',
!!!     &          obj(iobj)%obj_name,obj(iobj)%set,
!!!     &          obj(iobj)%bnd,obj(iobj)%vargrid

        if (obj(iobj)%np>0) then
          !! only for objects with a presences in this subdomain

          if (n2==3) then
            allocate(obj(iobj)%ang(obj(iobj)%np))
            i0 = obj(iobj)%start_idx
            i1 = i0+obj(iobj)%np-1
            obj(iobj)%ang = object(i0:i1,3)
          endif

          !!!! DIAGNOSTICS needs to be activated and diag_pflx=.true.
          !!!if (obj(iobj)%up.and. .not.diag_pflx) then
          !!!  print *,'Fatal ERROR: up extraction, diag_pflx is not set'
          !!!  stop
          !!!endif
          !!!if (obj(iobj)%vp.and. .not.diag_pflx) then
          !!!  print *,'Fatal ERROR: vp extraction, diag_pflx is not set'
          !!!  stop
          !!!endif

          !!!if (obj(iobj)%up.or.obj(iobj)%vp) extend_up = .true.

        endif

        deallocate(object)
      enddo
      ierr = nf90_close(ncid)

      !identify the number of sets
      !different sets are written in separate files (ncids allocation)
      allocate(setmask(nobj))
      setmask(:) = .false.
      do iobj = 1,nobj
        num = count( obj(iobj)%set==obj(:)%set )
        if (num==1) then
          setmask(iobj) = .true.
        else
            if(.not. any(obj(iobj)%set==obj(:)%set .and. setmask) ) setmask(iobj) = .true.
        endif
      enddo
!      
      allocate(unique_set(count(setmask)))
      allocate(ncidextr(count(setmask)))
      ncidextr(:) = -1
      unique_set = pack( obj(:)%set, setmask )
      deallocate(setmask)

      end subroutine read_extraction_objects !]
! ----------------------------------------------------------------------
      subroutine find_local_points(obj,gobj_i,gobj_j) ![

      ! Find object index locations that are within the subdomain
      ! Assign start and lenght of the local points in the global array
      ! of the object
      ! Translate global index locations to local ones
      implicit none
      ! import/export
      type(extract_object),intent(inout) ::obj
      real,dimension(:)   ,intent(inout) ::gobj_i,gobj_j ! global indices

      ! local
      integer :: i,start_idx,end_idx
      integer :: np
      
      np = size(gobj_i)
# ifdef MPI      
      gobj_i = gobj_i-iminmpi+1
      gobj_j = gobj_j-jminmpi+1
# endif

      ! Assume that local ranges of objects are contiguous
      start_idx = 0
      do i = 1,np
        if ( gobj_i(i)>=0.and.gobj_i(i)<LOCALLM .and. &
             gobj_j(i)>=0.and.gobj_j(i)<LOCALMM ) then
          start_idx = i
          exit
        endif
      enddo

      end_idx = np
      if (start_idx>0) then
        do i = start_idx,np
          if (gobj_i(i)<0.or.gobj_i(i)>=LOCALLM.or. &
              gobj_j(i)<0.or.gobj_j(i)>=LOCALMM ) then
            end_idx = i-1
            exit
          endif
        enddo
        obj%np = end_idx - start_idx + 1
      else ! object not in local range
        obj%np = 0
      endif

      if (obj%np>0) then
        obj%start_idx = start_idx
        obj%np = end_idx-start_idx+1
        allocate(obj%ipos(obj%np))
        allocate(obj%jpos(obj%np))
        obj%ipos = gobj_i(start_idx:end_idx)
        obj%jpos = gobj_j(start_idx:end_idx)
# ifdef NC4PAR
      else
        obj%start_idx = 1
# endif
      endif

      end subroutine find_local_points  !]
! ----------------------------------------------------------------------
      subroutine compute_coef(ipos,jpos,coef,ip,jp,mask)  ![
      ! compute interpolation coefficients
      implicit none

      ! inport/export
      real   ,dimension(:)  ,intent(in) :: ipos,jpos
      real   ,dimension(:,:),intent(out):: coef
      integer,dimension(:)  ,intent(out):: ip,jp
      real   ,dimension(GLOBAL_2D_ARRAY),intent(in), optional :: mask
      real :: coef_sum

      ! local
      integer :: i,np
      real :: cfx,cfy

      np = size(ip,1)
      do i = 1,np
        ip(i)  = floor(ipos(i))
        cfx    = ipos(i)-ip(i)
        jp(i)  = floor(jpos(i))
        cfy    = jpos(i)-jp(i)
# ifndef ETCH_INTO_LAND
        if(present(mask)) then
          coef(i,1) = (1-cfx)*(1-cfy)*mask(ip(i)  ,jp(i)  ) 
          coef(i,2) = cfx    *(1-cfy)*mask(ip(i)+1,jp(i)  ) 
          coef(i,3) = (1-cfx)*   cfy *mask(ip(i)  ,jp(i)+1)
          coef(i,4) =    cfx *   cfy *mask(ip(i)+1,jp(i)+1)
        else
# endif
          coef(i,1) = (1-cfx)*(1-cfy)
          coef(i,2) = cfx    *(1-cfy)
          coef(i,3) = (1-cfx)*   cfy 
          coef(i,4) =    cfx *   cfy 
# ifndef ETCH_INTO_LAND
        endif
# endif
        !! possibly check for all masked ....
        coef_sum = sum(coef(i,:))
        if(coef_sum/=0.) then
          coef(i,:) = coef(i,:)/sum(coef(i,:))
        else
          coef(i,:) = 0.
        endif
      enddo

      end subroutine compute_coef !]
! ----------------------------------------------------------------------
      subroutine interpolate_2D(var,vari,coef,ip,jp)  ![
      ! Interpolate a scalar variable
      !use dimensions
      implicit none

      ! inputs
      real   ,dimension(GLOBAL_2D_ARRAY),intent(in) :: var  ! assumed size arrays always start at 1.
      real   ,dimension(:)              ,intent(out):: vari
      real   ,dimension(:,:)            ,intent(in) :: coef
      integer,dimension(:)              ,intent(in) :: ip
      integer,dimension(:)              ,intent(in) :: jp

      ! local
      integer :: i,k,np

      np = size(ip,1)
      do i = 1,np
          vari(i) = var(ip(i)  ,jp(i)  )*coef(i,1) + &
                    var(ip(i)+1,jp(i)  )*coef(i,2) + &
                    var(ip(i)  ,jp(i)+1)*coef(i,3) + &
                    var(ip(i)+1,jp(i)+1)*coef(i,4)
      enddo

      end subroutine interpolate_2D  !]
! ----------------------------------------------------------------------
      subroutine interpolate_3D(var,vari,coef,ip,jp)  ![
      ! Interpolate a variable
      !use dimensions
      implicit none

      ! inputs
      real   ,dimension(GLOBAL_2D_ARRAY,N),intent(in) :: var  ! assumed size arrays would always start at 1.
      real   ,dimension(:,:)               ,intent(out):: vari
      real   ,dimension(:,:)               ,intent(in) :: coef
      integer,dimension(:)                 ,intent(in) :: ip
      integer,dimension(:)                 ,intent(in) :: jp

      ! local
      integer :: i,k,np


!     if (mynode==1) then
!       print *,'interp'
!       print *,shape(var)
!       print *, ip(10),jp(10)
!       print *, nx,ny
!       print *, var(ip(10),jp(10),nz)
!     endif
      np = size(ip,1)
      do i = 1,np
        do k = 1,N
          vari(i,k) = var(ip(i)  ,jp(i)  ,k)*coef(i,1) + &
                      var(ip(i)+1,jp(i)  ,k)*coef(i,2) + &
                      var(ip(i)  ,jp(i)+1,k)*coef(i,3) + &
                      var(ip(i)+1,jp(i)+1,k)*coef(i,4)
        enddo
      enddo

      end subroutine interpolate_3D  !]
! ----------------------------------------------------------------------
      subroutine do_extract_bdy ![
      ! extract data for all objects, for all vars
      ! and write to file
      implicit none

      ! local
      integer :: i,itrc,ierr,ncid,k,record,ifile,nset
      character(len=30) :: obj_name
      character(len=99),save :: fname
      character(len=20)              :: tname
      character(len=40) :: oname
      integer :: lpre
      real, dimension(:,:),pointer :: vi
      real, dimension(:,:),pointer :: ui
      real, dimension(:,:),pointer :: coef, cfu, cfv
      integer,dimension(:),pointer :: ip,jp,ipu,ipv,jpu,jpv
      real,dimension(:),allocatable :: dummy
      integer,dimension(3) :: start2D, count2D
      integer,dimension(2) :: start1D, count1D
      real   ,dimension(GLOBAL_2D_ARRAY,N) :: var
      integer :: imin, imax, jmin, jmax

      if (.not.allocated(obj)) then
        call init_extract_data
        !call def_extract
      endif

# if defined MPI & !defined NC4PAR
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1, 1, MPI_COMM_WORLD, status, ierr)
      endif
# endif

      call def_extract
!
      nrecextr = max(1,nrecextr)
      record = nrecextr
!      
      nset = size(unique_set)

      imin = 1
      imax = LOCALLM
      jmin = 1
      jmax = LOCALMM

      if (WEST_INTER) then
         imin = imin - 1
      endif
      if (EAST_INTER) then
         imax = imax + 1
      endif
      if (SOUTH_INTER) then
         jmin = jmin - 1
      endif
      if (NORTH_INTER) then
         jmax = jmax + 1
      endif

      do ifile=1,nset
!      
        do i = 1,nobj
!
          if(trim(obj(i)%set) == trim(unique_set(ifile))) then

# ifndef NC4PAR
            if (obj(i)%np>0) then
# endif
              start1D = (/obj(i)%start_idx,record/)
              start2D = (/obj(i)%start_idx,1,record/)
              count1D = (/obj(i)%np,1/)
              count2D = (/obj(i)%np,N,1/)
!
              coef => obj(i)%coef
              ip   => obj(i)%ip
              jp   => obj(i)%jp
  
              if (.not.obj(i)%scalar) then
                ui => obj(i)%ui
                vi => obj(i)%vi
                cfu=> obj(i)%cfu
                cfv=> obj(i)%cfv
                ipu=> obj(i)%ipu
                ipv=> obj(i)%ipv
                jpu=> obj(i)%jpu
                jpv=> obj(i)%jpv
              endif

              ierr = nf90_put_var(ncidextr(ifile),extrTime(ifile), &
# ifdef USE_CALENDAR
                                  (/time - origin_date_in_sec/), &
# else
                                  (/time*sec2day/), &
# endif
                                  (/record/),(/1/) )
              call check(ierr,'problem in writing the variable: time')

              if (obj(i)%zeta) then
# ifdef NC4PAR
                if (obj(i)%np>0) then
# endif
# ifdef ETCH_INTO_LAND
                var(imin:imax,jmin:jmax,1) = zeta(imin:imax,jmin:jmax,fast_indx_out)
                call etch_into_land(int(rmask(imin:imax,jmin:jmax)), var(imin:imax,jmin:jmax,1))
                call interpolate(var(:,:,1), obj(i)%vari(:,1),coef, ip, jp)
# else
                call interpolate(zeta(:,:,fast_indx_out), obj(i)%vari(:,1), coef, ip, jp)
# endif
                where(sum(coef,2)==0.) obj(i)%vari(:,1)=NF90_FILL_DOUBLE
# ifdef NC4PAR
                endif                
# endif
                select case(trim(obj(i)%bnd))
                case('north')
                  ierr = nf90_put_var(ncidextr(ifile), extrZ_N(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'zeta_'//trim(obj(i)%bnd))
                case('south')
                  ierr = nf90_put_var(ncidextr(ifile), extrZ_S(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'zeta_'//trim(obj(i)%bnd))
                case('west')
                  !print*,'zeta_west',mynode,start1D,count1D
                  ierr = nf90_put_var(ncidextr(ifile), extrZ_W(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'zeta_'//trim(obj(i)%bnd))
                case('east')
                  ierr = nf90_put_var(ncidextr(ifile), extrZ_E(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'zeta_'//trim(obj(i)%bnd))
                end select
              endif
!
              if (obj(i)%ubar) then
# ifdef NC4PAR
                if (obj(i)%np>0) then
# endif
# ifdef ETCH_INTO_LAND
                var(imin:imax,jmin:jmax,1) = ubar(imin:imax,jmin:jmax,fast_indx_out)
                call etch_into_land(int(umask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,1))
                call interpolate(var(:,:,1), ui(:,1), cfu, ipu, jpu)
                var(imin:imax,jmin:jmax,1) = vbar(imin:imax,jmin:jmax,fast_indx_out)
                call etch_into_land(int(vmask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,1))
                call interpolate(var(:,:,1), vi(:,1), cfu, ipu, jpu)
# else
                call interpolate(ubar(:,:,fast_indx_out), ui(:,1), cfu, ipu, jpu)
                call interpolate(vbar(:,:,fast_indx_out), vi(:,1), cfv, ipv, jpv)
# endif
                obj(i)%vari(:,1) = obj(i)%cosa*ui(:,1) &
                                 - obj(i)%sina*vi(:,1)
                where(sum(coef,2)==0.) obj(i)%vari(:,1)=NF90_FILL_DOUBLE
# ifdef NC4PAR
                endif                
# endif
                select case(trim(obj(i)%bnd))
                case('north')
                  ierr = nf90_put_var(ncidextr(ifile), extrUb_N(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'ubar_'//trim(obj(i)%bnd))
                case('south')
                  ierr = nf90_put_var(ncidextr(ifile), extrUb_S(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'ubar_'//trim(obj(i)%bnd))
                case('west')
                  ierr = nf90_put_var(ncidextr(ifile), extrUb_W(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'ubar_'//trim(obj(i)%bnd))
                case('east')
                  ierr = nf90_put_var(ncidextr(ifile), extrUb_E(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'ubar_'//trim(obj(i)%bnd))
                end select
              endif
!
              if (obj(i)%vbar) then
# ifdef NC4PAR
                if (obj(i)%np>0) then
# endif
# ifdef ETCH_INTO_LAND
                var(imin:imax,jmin:jmax,1) = ubar(imin:imax,jmin:jmax,fast_indx_out)
                call etch_into_land(int(umask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,1))
                call interpolate(var(:,:,1), ui(:,1), cfu, ipu, jpu)
                var(imin:imax,jmin:jmax,1) = vbar(imin:imax,jmin:jmax,fast_indx_out)
                call etch_into_land(int(vmask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,1))
                call interpolate(var(:,:,1), vi(:,1), cfu, ipu, jpu)
# else
                call interpolate(ubar(:,:,fast_indx_out), ui(:,1), cfu, ipu, jpu)
                call interpolate(vbar(:,:,fast_indx_out), vi(:,1), cfv, ipv, jpv)
# endif
                obj(i)%vari(:,1) = obj(i)%sina*ui(:,1) &
                                 + obj(i)%cosa*vi(:,1)
                where(sum(coef,2)==0.) obj(i)%vari(:,1)=NF90_FILL_DOUBLE
# ifdef NC4PAR
                endif                
# endif
                select case(trim(obj(i)%bnd))
                case('north')
                  ierr = nf90_put_var(ncidextr(ifile), extrVb_N(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'vbar_'//trim(obj(i)%bnd))
                case('south')
                  ierr = nf90_put_var(ncidextr(ifile), extrVb_S(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'vbar_'//trim(obj(i)%bnd))
                case('west')
                  ierr = nf90_put_var(ncidextr(ifile), extrVb_W(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'vbar_'//trim(obj(i)%bnd))
                case('east')
                  ierr = nf90_put_var(ncidextr(ifile), extrVb_E(ifile), obj(i)%vari(:,1), start1D, count1D)
                  call check(ierr,'problem in writing variable: '//'vbar_'//trim(obj(i)%bnd))
                end select
              endif
!
# ifdef SOLVE3D              
              if (obj(i)%temp) then
#  ifdef NC4PAR
                if (obj(i)%np>0) then
#  endif
#  ifdef ETCH_INTO_LAND
                var(imin:imax,jmin:jmax,:) = t(imin:imax,jmin:jmax,:,nstp,itemp)
                do k=1,N
                  call etch_into_land(int(rmask(imin:imax,jmin:jmax)), var(imin:imax,jmin:jmax,k))
                enddo  
                call interpolate(var(:,:,:), obj(i)%vari,coef, ip, jp)
#  else
                call interpolate(t(:,:,:,nstp,itemp), obj(i)%vari, coef, ip, jp)
#  endif
                do k=1,N
                  where(sum(coef,2)==0.) obj(i)%vari(:,k)=NF90_FILL_DOUBLE
                enddo
#  ifdef NC4PAR
                endif                
#  endif
                select case(trim(obj(i)%bnd))
                case('north')
                  ierr = nf90_put_var(ncidextr(ifile), extrT_N(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'temp_'//trim(obj(i)%bnd))
                case('south')
                  ierr = nf90_put_var(ncidextr(ifile), extrT_S(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'temp_'//trim(obj(i)%bnd))
                case('west')
                  ierr = nf90_put_var(ncidextr(ifile), extrT_W(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'temp_'//trim(obj(i)%bnd))
                case('east')
                  ierr = nf90_put_var(ncidextr(ifile), extrT_E(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'temp_'//trim(obj(i)%bnd))
                end select
              endif
!              
#  ifdef SALINITY
              if (obj(i)%salt) then
#   ifdef NC4PAR
                if (obj(i)%np>0) then
#   endif
#   ifdef ETCH_INTO_LAND
                var(imin:imax,jmin:jmax,:) = t(imin:imax,jmin:jmax,:,nstp,isalt)
                do k=1,N
                  call etch_into_land(int(rmask(imin:imax,jmin:jmax)), var(imin:imax,jmin:jmax,k))
                enddo  
                call interpolate(var(:,:,:), obj(i)%vari,coef, ip, jp)
#   else
                call interpolate(t(:,:,:,nstp,isalt), obj(i)%vari, coef, ip, jp)
#   endif
                do k=1,N
                  where(sum(coef,2)==0.) obj(i)%vari(:,k)=NF90_FILL_DOUBLE
                enddo
#   ifdef NC4PAR
                endif                
#   endif
                select case(trim(obj(i)%bnd))
                case('north')
                  ierr = nf90_put_var(ncidextr(ifile), extrS_N(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'salt_'//trim(obj(i)%bnd))
                case('south')
                  ierr = nf90_put_var(ncidextr(ifile), extrS_S(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'salt_'//trim(obj(i)%bnd))
                case('west')
                  ierr = nf90_put_var(ncidextr(ifile), extrS_W(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'salt_'//trim(obj(i)%bnd))
                case('east')
                  ierr = nf90_put_var(ncidextr(ifile), extrS_E(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'salt_'//trim(obj(i)%bnd))
                end select
              endif
#  endif
              if (obj(i)%u) then
#  ifdef NC4PAR
                if (obj(i)%np>0) then
#  endif
#  ifdef ETCH_INTO_LAND
                var(imin:imax,jmin:jmax,:) = u(imin:imax,jmin:jmax,:,nstp)
                do k=1,N
                  call etch_into_land(int(umask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,k))
                enddo
                call interpolate(var(:,:,:),ui,cfu,ipu,jpu)
                var(imin:imax,jmin:jmax,:) = v(imin:imax,jmin:jmax,:,nstp)
                do k=1,N
                  call etch_into_land(int(vmask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,k))
                enddo
                call interpolate(var(:,:,:),vi,cfv,ipv,jpv)
#  else
                call interpolate(u(:,:,:,nstp),ui,cfu,ipu,jpu)
                call interpolate(v(:,:,:,nstp),vi,cfv,ipv,jpv)
                !obj(i)%vari = ui
#  endif                
                do k=1,N
                  obj(i)%vari(:,k) = obj(i)%cosa*ui(:,k) &
                                   - obj(i)%sina*vi(:,k)
                  where(sum(coef,2)==0.) obj(i)%vari(:,k)=NF90_FILL_DOUBLE               
                enddo
#  ifdef NC4PAR
                endif                
#  endif
                select case(trim(obj(i)%bnd))
                case('north')
                  ierr = nf90_put_var(ncidextr(ifile), extrU_N(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'u_'//trim(obj(i)%bnd))
                case('south')
                  ierr = nf90_put_var(ncidextr(ifile), extrU_S(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'u_'//trim(obj(i)%bnd))
                case('west')
                  ierr = nf90_put_var(ncidextr(ifile), extrU_W(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'u_'//trim(obj(i)%bnd))
                case('east')
                  ierr = nf90_put_var(ncidextr(ifile), extrU_E(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'u_'//trim(obj(i)%bnd))
                end select
              endif
!              
              if (obj(i)%v) then
#  ifdef NC4PAR
                if (obj(i)%np>0) then
#  endif
#  ifdef ETCH_INTO_LAND
                var(imin:imax,jmin:jmax,:) = u(imin:imax,jmin:jmax,:,nstp)
                do k=1,N
                  call etch_into_land(int(umask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,k))
                enddo
                call interpolate(var(:,:,:),ui,cfu,ipu,jpu)
                var(imin:imax,jmin:jmax,:) = v(imin:imax,jmin:jmax,:,nstp)
                do k=1,N
                  call etch_into_land(int(vmask(imin:imax,jmin:jmax)),var(imin:imax,jmin:jmax,k))
                enddo
                call interpolate(var(:,:,:),vi,cfv,ipv,jpv)
#  else
                call interpolate(u(:,:,:,nstp),ui,cfu,ipu,jpu)
                call interpolate(v(:,:,:,nstp),vi,cfv,ipv,jpv)
                !obj(i)%vari = ui
#  endif 
                do k=1,N
                  obj(i)%vari(:,k) = obj(i)%sina*ui(:,k) &
                                   + obj(i)%cosa*vi(:,k)
                  where(sum(coef,2)==0.) obj(i)%vari(:,k)=NF90_FILL_DOUBLE
                enddo
#  ifdef NC4PAR
                endif                
#  endif
                select case(trim(obj(i)%bnd))
                case('north')
                  ierr = nf90_put_var(ncidextr(ifile), extrV_N(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'v_'//trim(obj(i)%bnd))
                case('south')
                  ierr = nf90_put_var(ncidextr(ifile), extrV_S(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'v_'//trim(obj(i)%bnd))
                case('west')
                  ierr = nf90_put_var(ncidextr(ifile), extrV_W(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'v_'//trim(obj(i)%bnd))
                case('east')
                  ierr = nf90_put_var(ncidextr(ifile), extrV_E(ifile), obj(i)%vari, start2D, count2D)
                  call check(ierr,'problem in writing variable: '//'v_'//trim(obj(i)%bnd))
                end select
              endif
# endif
!
# ifndef NC4PAR
            endif !obj(i)%np>0
# endif
!            
          endif
        enddo

# if defined MPI & !defined NC4PAR
        ierr = nf90_close(ncidextr(ifile))
# else
        ierr = nf90_sync(ncidextr(ifile))
# endif

      enddo
!
      MPI_master_only write(stdout,'(6x,A,2(A,I4,1x),A,I3)') &
                           'WRT_EXTR -- wrote ', &
                           'boundary fields into time record =', &
                            record, '/', nrecextr  MYID
!      
# if defined MPI & !defined NC4PAR
      if (mynode .lt. NNODES-1) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1, 1, MPI_COMM_WORLD, ierr)
      endif
# endif
      
      return

      end subroutine do_extract_bdy  !]
! ----------------------------------------------------------------------
      subroutine def_extract ![
      implicit none

      !local
      integer :: ierr
      logical :: create_new_file
      integer :: lenstr, lstr, lvar, ifile, iobj, nset
      character*180 :: fname
      integer,dimension(1) :: timedim
      integer,dimension(2) :: r2dx, r2dy, &
                              u2dx, u2dy, &
                              v2dx, v2dy 
# ifdef SOLVE3D
      integer,dimension(3) :: r3dx, r3dy, &
                              u3dx, u3dy, &
                              v3dx, v3dy 
# endif
# ifdef NC4PAR
      integer :: cmode
# endif
!
      nset = size(unique_set)
      if (.not.allocated(extrTime)) then
        allocate(extrTime(nset), &
                 extrZ_N(nset) , extrZ_S(nset) , &
                 extrZ_W(nset) , extrZ_E(nset) , &
                 extrUb_N(nset), extrUb_S(nset), &
                 extrUb_W(nset), extrUb_E(nset), &
                 extrVb_N(nset), extrVb_S(nset), &
                 extrVb_W(nset), extrVb_E(nset), &
                 extrU_N(nset) , extrU_S(nset) , &
                 extrU_W(nset) , extrU_E(nset) , &
                 extrV_N(nset) , extrV_S(nset) , &
                 extrV_W(nset) , extrV_E(nset) , &
                 extrT_N(nset) , extrT_S(nset) , &
                 extrT_W(nset) , extrT_E(nset) , &
                 extrS_N(nset) , extrS_S(nset) , &
                 extrS_W(nset) , extrS_E(nset)  )
      endif
!      
      do ifile=1,nset
!
!!!        write(stdout,*) 'nset= ',nset, mynode, 
!!!     &                   lenstr(unique_set(ifile))
        create_new_file = .true.
        lstr=lenstr(out_extrname)
        if(nset.eq.1) then
          fname=out_extrname(1:lstr)
        else
          fname=out_extrname(1:lstr)//'.'//trim(unique_set(ifile))
        endif
        lstr=lenstr(fname)
        if (ncidextr(ifile).ne.-1) create_new_file=.false.
# if defined MPI & !defined NC4PAR
        if (mynode.gt.0) create_new_file = .false.
# endif
!
! Create output file 
! -----------------------------------------------------------------------
        if (create_new_file) then
# ifndef NC4PAR         
          ierr=nf90_create(fname(1:lstr),NF90_CLOBBER, ncidextr(ifile))
# else
          cmode = ior(nf90_netcdf4, nf90_mpiio)
          ierr=nf90_create(fname(1:lstr), cmode, ncidextr(ifile), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
# endif          
          if (ierr .ne. nf90_noerr) then
            MPI_master_only write(stdout,'(/2(1x,A)/)') &
                                  'ERROR in def_extract: Cannot create netCDF file:', &
                                   fname(1:lstr)
            goto 99                                         !--> ERROR
          endif
! 
! Define dimensions
! -----------------------------------------------------------------------          
          ierr=nf90_def_dim (ncidextr(ifile), 'bry_time', nf90_unlimited, timedim(1))
          do iobj = 1,nobj
!
            if(trim(obj(iobj)%set) == trim(unique_set(ifile))) then
!          
              if(trim(obj(iobj)%vargrid) == 'rho') then
                if(trim(obj(iobj)%bnd) == 'south' &
              .or. trim(obj(iobj)%bnd) == 'north') then
                  ierr = nf90_inq_dimid(ncidextr(ifile),'xi_rho',r2dx(1))
                  if(ierr/=0) then
                    ierr=nf90_def_dim(ncidextr(ifile),'xi_rho',obj(iobj)%dsize,r2dx(1))
                  endif
                else
                  ierr = nf90_inq_dimid(ncidextr(ifile),'eta_rho',r2dy(1))
                  if(ierr/=0) then
                    ierr=nf90_def_dim(ncidextr(ifile),'eta_rho',obj(iobj)%dsize,r2dy(1))
                  endif
                endif  
!
              elseif(trim(obj(iobj)%vargrid) == 'u') then
                if(trim(obj(iobj)%bnd) == 'south' &
              .or. trim(obj(iobj)%bnd) == 'north') then
                  ierr = nf90_inq_dimid(ncidextr(ifile),'xi_u',u2dx(1))
                  if(ierr/=0) then
                    ierr=nf90_def_dim(ncidextr(ifile),'xi_u',obj(iobj)%dsize,u2dx(1))
                  endif
                else
                  ierr = nf90_inq_dimid(ncidextr(ifile),'eta_rho',u2dy(1))
                  if(ierr/=0) then
                    ierr=nf90_def_dim(ncidextr(ifile),'eta_rho',obj(iobj)%dsize,u2dy(1))
                  endif
                endif  
!              
              elseif(trim(obj(iobj)%vargrid) == 'v') then
                if(trim(obj(iobj)%bnd) == 'south' &
              .or. trim(obj(iobj)%bnd) == 'north') then
                  ierr = nf90_inq_dimid(ncidextr(ifile),'xi_rho',v2dx(1))
                  if(ierr/=0) then
                    ierr=nf90_def_dim(ncidextr(ifile),'xi_rho',obj(iobj)%dsize,v2dx(1))
                  endif
                else
                  ierr = nf90_inq_dimid(ncidextr(ifile),'eta_v',v2dy(1))
                  if(ierr/=0) then
                    ierr=nf90_def_dim(ncidextr(ifile),'eta_v',obj(iobj)%dsize,v2dy(1))
                  endif
                endif 
              endif
            endif  
          enddo  
!
          r2dx(2) = timedim(1)
          u2dx(2) = timedim(1)
          v2dx(2) = timedim(1)
          r2dy(2) = timedim(1)
          u2dy(2) = timedim(1)
          v2dy(2) = timedim(1)
# ifdef SOLVE3D
          r3dx(1) = r2dx(1)
          r3dy(1) = r2dy(1)
          u3dx(1) = u2dx(1)
          u3dy(1) = u2dy(1)
          v3dx(1) = v2dx(1)
          v3dy(1) = v2dy(1)

          ierr=nf90_def_dim (ncidextr(ifile), 's_rho', N, r3dx(2))
          r3dy(2) = r3dx(2)
          u3dx(2) = r3dx(2)
          u3dy(2) = r3dx(2)
          v3dx(2) = r3dx(2)
          v3dy(2) = r3dx(2)

          r3dx(3) = timedim(1)
          u3dx(3) = timedim(1)
          v3dx(3) = timedim(1)
          r3dy(3) = timedim(1)
          u3dy(3) = timedim(1)
          v3dy(3) = timedim(1)
# endif          
!
! Define variables
!------------------------------------------------------------------------
          ierr=nf90_def_var(ncidextr(ifile), 'bry_time', NF90_DOUBLE, timedim, extrTime(ifile))
# ifdef NC4PAR
          ierr=nf90_var_par_access(ncidextr(ifile),extrTime(ifile),nf90_collective)
          !ierr=nf90_var_par_access(ncidextr(ifile),extrTime(ifile),nf90_independent)
# endif         
          lvar=lenstr(vname(2,indxTime))
          ierr=nf90_put_att(ncidextr(ifile),extrTime(ifile),'long_name',vname(2,indxTime)(1:lvar) )
# ifdef USE_CALENDAR
          lvar=lenstr(vname(3,indxTime))
          ierr=nf90_put_att (ncidextr(ifile),extrTime(ifile),'units',vname(3,indxTime)(1:lvar) )
# else
          ierr=nf90_put_att (ncidextr(ifile),extrTime(ifile),'units','days' )
# endif          
          lvar=lenstr(vname(4,indxTime))
          ierr=nf90_put_att (ncidextr(ifile),extrTime(ifile),'field',vname(4,indxTime)(1:lvar) )
!          
          do iobj = 1,nobj
!          
            if(trim(obj(iobj)%set) == trim(unique_set(ifile))) then
!
              !if(trim(obj(iobj)%bnd) == 'north') then
              select case(trim(obj(iobj)%bnd))
              case('north')  
!                  
                if(obj(iobj)%zeta) then
                  call create_var(ncidextr(ifile),'zeta_north',r2dx,extrZ_N(ifile),indxZ)
                endif  
!                
                if(obj(iobj)%ubar) then
                  call create_var(ncidextr(ifile),'ubar_north',u2dx,extrUb_N(ifile),indxUb)
                endif  
!                
                if(obj(iobj)%vbar) then
                  call create_var(ncidextr(ifile),'vbar_north',v2dx,extrVb_N(ifile),indxVb)
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  call create_var(ncidextr(ifile),'u_north',u3dx,extrU_N(ifile),indxU)
                endif
!                
                if(obj(iobj)%v) then
                  call create_var(ncidextr(ifile),'v_north',v3dx,extrV_N(ifile),indxV)
                endif
!                
                if(obj(iobj)%temp) then
                  call create_var(ncidextr(ifile),'temp_north',r3dx,extrT_N(ifile),indxT)
                endif
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  call create_var(ncidextr(ifile),'salt_north',r3dx,extrS_N(ifile),indxS)
                endif 
#  endif
# endif                

              !elseif(trim(obj(iobj)%bnd) == 'south') then
              case('south')
!                  
                if(obj(iobj)%zeta) then
                  call create_var(ncidextr(ifile),'zeta_south',r2dx,extrZ_S(ifile),indxZ)
                endif  
!                
                if(obj(iobj)%ubar) then
                  call create_var(ncidextr(ifile),'ubar_south',u2dx,extrUb_S(ifile),indxUb)
                endif  
!                
                if(obj(iobj)%vbar) then
                  call create_var(ncidextr(ifile),'vbar_south',v2dx,extrVb_S(ifile),indxVb)
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  call create_var(ncidextr(ifile),'u_south',u3dx,extrU_S(ifile),indxU)
                endif
!                
                if(obj(iobj)%v) then
                  call create_var(ncidextr(ifile),'v_south',v3dx,extrV_S(ifile),indxV)
                endif
!                
                if(obj(iobj)%temp) then
                  call create_var(ncidextr(ifile),'temp_south',r3dx,extrT_S(ifile),indxT)
                endif
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  call create_var(ncidextr(ifile),'salt_south',r3dx,extrS_S(ifile),indxS)
                endif 
#  endif
# endif
              !elseif(trim(obj(iobj)%bnd) == 'west') then
              case('west')
!                  
                if(obj(iobj)%zeta) then
                  call create_var(ncidextr(ifile),'zeta_west',r2dy,extrZ_W(ifile),indxZ)
                endif  
!                
                if(obj(iobj)%ubar) then
                  call create_var(ncidextr(ifile),'ubar_west',u2dy,extrUb_W(ifile),indxUb)
                endif  
!                
                if(obj(iobj)%vbar) then
                  call create_var(ncidextr(ifile),'vbar_west',v2dy,extrVb_W(ifile),indxVb)
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  call create_var(ncidextr(ifile),'u_west',u3dy,extrU_W(ifile),indxU)
                endif
!                
                if(obj(iobj)%v) then
                  call create_var(ncidextr(ifile),'v_west',v3dy,extrV_W(ifile),indxV)
                endif
!                
                if(obj(iobj)%temp) then
                  call create_var(ncidextr(ifile),'temp_west',r3dy,extrT_W(ifile),indxT)
                endif
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  call create_var(ncidextr(ifile),'salt_west',r3dy,extrS_W(ifile),indxS)
                endif 
#  endif
# endif
              !elseif(trim(obj(iobj)%bnd) == 'east') then
              case('east')
!                  
                if(obj(iobj)%zeta) then
                  call create_var(ncidextr(ifile),'zeta_east',r2dy,extrZ_E(ifile),indxZ)
                endif  
!                
                if(obj(iobj)%ubar) then
                  call create_var(ncidextr(ifile),'ubar_east',u2dy,extrUb_E(ifile),indxUb)
                endif  
!                
                if(obj(iobj)%vbar) then
                  call create_var(ncidextr(ifile),'vbar_east',v2dy,extrVb_E(ifile),indxVb)
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  call create_var(ncidextr(ifile),'u_east',u3dy,extrU_E(ifile),indxU)
                endif
!                
                if(obj(iobj)%v) then
                  call create_var(ncidextr(ifile),'v_east',v3dy,extrV_E(ifile),indxV)
                endif
!                
                if(obj(iobj)%temp) then
                  call create_var(ncidextr(ifile),'temp_east',r3dy,extrT_E(ifile),indxT)
                endif
!                
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  call create_var(ncidextr(ifile),'salt_east',r3dy,extrS_E(ifile),indxS)
                endif 
#  endif
# endif
              !endif
              end select


            endif
          enddo
          ierr = nf90_enddef(ncidextr(ifile))
          MPI_master_only write(stdout,'(6x,4A,1x,A,i4)') &
                        'DEF_EXTRACT - Created ', &
                        'new netCDF file ''', fname(1:lstr), '''.' &
                         MYID
!
        elseif(ncidextr(ifile)==-1) then
!
# ifndef NC4PAR            
          ierr = nf90_open(fname(1:lstr),NF90_WRITE, ncidextr(ifile))
# else
          ierr = nf90_open(fname(1:lstr),IOR(nf90_write, nf90_mpiio), ncidextr(ifile), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL)
# endif
!          
          do iobj = 1,nobj
!          
            if(trim(obj(iobj)%set) == trim(unique_set(ifile))) then
!
              if(trim(obj(iobj)%bnd) == 'north') then
!                  
                if(obj(iobj)%zeta) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'zeta_north',extrZ_N(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrZ_N(ifile),nf90_collective)
# endif
                endif  
!                
                if(obj(iobj)%ubar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'ubar_north',extrUb_N(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrUb_N(ifile),nf90_collective)
# endif
                 
                endif  
!                
                if(obj(iobj)%vbar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'vbar_north',extrVb_N(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrVb_N(ifile),nf90_collective)
# endif
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'u_north',extrU_N(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrU_N(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%v) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'v_north',extrV_N(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrV_N(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%temp) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'temp_north',extrT_N(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrT_N(ifile),nf90_collective)
#  endif
                endif
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'salt_north',extrS_N(ifile))
#   ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrS_N(ifile),nf90_collective)
#   endif
                endif 
#  endif
# endif                

              elseif(trim(obj(iobj)%bnd) == 'south') then
!                  
                if(obj(iobj)%zeta) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'zeta_south',extrZ_S(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrZ_S(ifile),nf90_collective)
# endif
                endif  
!                
                if(obj(iobj)%ubar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'ubar_south',extrUb_S(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrUb_S(ifile),nf90_collective)
# endif
                endif  
!                
                if(obj(iobj)%vbar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'vbar_south',extrVb_S(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrVb_S(ifile),nf90_collective)
# endif
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'u_south',extrU_S(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrU_S(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%v) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'v_south',extrV_S(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrV_S(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%temp) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'temp_south',extrT_S(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrT_S(ifile),nf90_collective)
#  endif
                endif
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'salt_south',extrS_S(ifile))
#   ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrS_S(ifile),nf90_collective)
#   endif
                endif 
#  endif
# endif
              elseif(trim(obj(iobj)%bnd) == 'west') then
!                  
                if(obj(iobj)%zeta) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'zeta_west',extrZ_W(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrZ_W(ifile),nf90_collective)
# endif
                endif  
!                
                if(obj(iobj)%ubar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'ubar_west',extrUb_W(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrUb_W(ifile),nf90_collective)
# endif
                endif  
!                
                if(obj(iobj)%vbar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'vbar_west',extrVb_W(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrVb_W(ifile),nf90_collective)
# endif
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'u_west',extrU_W(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrU_W(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%v) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'v_west',extrV_W(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrV_W(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%temp) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'temp_west',extrT_W(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrT_W(ifile),nf90_collective)
#  endif
                endif
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'salt_west',extrS_W(ifile))
#   ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrS_W(ifile),nf90_collective)
#   endif
                endif 
#  endif
# endif
              elseif(trim(obj(iobj)%bnd) == 'east') then
!                  
                if(obj(iobj)%zeta) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'zeta_east',extrZ_E(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrZ_E(ifile),nf90_collective)
# endif
                endif  
!                
                if(obj(iobj)%ubar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'ubar_east',extrUb_E(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrUb_E(ifile),nf90_collective)
# endif
                endif  
!                
                if(obj(iobj)%vbar) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'vbar_east',extrVb_E(ifile))
# ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrVb_E(ifile),nf90_collective)
# endif
                endif  
# ifdef SOLVE3D
                if(obj(iobj)%u) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'u_east',extrU_E(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrU_E(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%v) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'v_east',extrV_E(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrV_E(ifile),nf90_collective)
#  endif
                endif
!                
                if(obj(iobj)%temp) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'temp_east',extrT_E(ifile))
#  ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrT_E(ifile),nf90_collective)
#  endif
                endif
!                
#  ifdef SALINITY                
                if(obj(iobj)%salt) then
                  ierr = nf90_inq_varid(ncidextr(ifile),'salt_east',extrS_E(ifile))
#   ifdef NC4PAR
                  ierr = nf90_var_par_access(ncidextr(ifile),extrS_E(ifile),nf90_collective)
#   endif
                endif 
#  endif
# endif
              endif


            endif
          enddo

# if defined MPI & !defined NC4PAR
        else
          ierr = nf90_open(fname(1:lstr),NF90_WRITE, ncidextr(ifile))
          if (ierr .ne. nf90_noerr) then
            MPI_master_only write(stdout,'(/1x,4A,2x,A,I4/)') &
                                  'DEF_EXTRACT ERROR: ', &
                                  'Cannot open file ''', fname(1:lstr), '''.' &
                                   MYID
            goto 99
          endif
# endif
        endif !create_new_file
!
        ierr=nf90_set_fill (ncidextr(ifile), nf90_nofill, lvar)
!
      enddo  

  99  return

      end subroutine def_extract !]
!
!----------------------------------------------------------------------------------
      subroutine create_var(ncid,varname,dims,varid,indx)
        
        implicit none

        integer, intent(in) :: ncid
        character(len=*), intent(in) :: varname
        integer, dimension(:), intent(in) :: dims
        integer, intent(out) :: varid
        integer, intent(in) :: indx
        integer :: lvar, lenstr, ierr

        ierr=nf90_def_var(ncid,varname,NF90_OUT,dims,varid)
        call check(ierr,'unable to create variable: '//varname)
# ifdef NC4PAR
        ierr=nf90_var_par_access(ncid,varid,nf90_collective)
        !ierr=nf90_var_par_access(ncid,varid,nf90_independent)
        call check(ierr,'unable to create variable: '//varname)
# endif 
        lvar=lenstr(vname(2,indx))
        ierr=nf90_put_att(ncid,varid,'long_name',vname(2,indx)(1:lvar) )
        call check(ierr,'unable to create variable: '//varname)
        lvar=lenstr(vname(3,indx))
        ierr=nf90_put_att (ncid,varid,'units',vname(3,indx)(1:lvar) )
        call check(ierr,'unable to create variable: '//varname)
        lvar=lenstr(vname(4,indx))
        ierr=nf90_put_att (ncid,varid,'field',vname(4,indx)(1:lvar) )
        call check(ierr,'unable to create variable: '//varname)
        ierr=nf90_put_att (ncid,varid,'_FillValue',NF90_FILL)

        return
      end subroutine 
!
! ------------------------------------------------------------------------------------------
# ifdef MASKING
#  define spv -9.9E+9
#  define spv_set -1.D+10
!
      subroutine etch_into_land(mask, qfld)
      implicit none
      integer nx,ny
      integer, dimension(:,:) :: mask
      real, dimension(:,:)    :: qfld
!>
      integer ntrds,trd, nmsk,nmsk_new, ncst, npths, &
              istr,iend,jstr,jend, j,j0, jskip, ierr
      integer, parameter :: jsize=3
      nx = size(mask,1)
      ny = size(mask,2)
      ntrds=1 ; trd=0
      j0=trd*jsize ; jskip=ntrds*jsize ; istr=1 ; iend=nx

      mskd_pts=1!<-- set to start while loop
      if ((nx+2)*(ny+2) > allc_ext_size) then
        allc_ext_size=(nx+2)*(ny+2)
        if (allocated(qext)) deallocate(qext)
        allocate(qext(nx+2,ny+2))
!!!        write(*,'(1x,2A,F16.8,1x,A)') 'etch_into_land_thread :: ',
!!!     &             'allocated',  dble((nx+2)*(ny+2))/dble(262144),
!!!     &                                'MB shared workspace array'
      endif

      nmsk=0
      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call copy_extend_tile(istr,iend,jstr,jend, nx,ny, mask, qfld, qext, nmsk)
        call set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
      enddo

      !print*, 'nmsk:', mynode, nmsk
      if(nmsk>0) then
        if (nmsk > alloc_msk_size) then
          alloc_msk_size=nmsk
          if (allocated(ptch)) deallocate(ptch,ijptch, nijmsk,ijmsk)
          allocate(  ijmsk(2,alloc_msk_size),  nijmsk(alloc_msk_size), &
                    ijptch(2,alloc_msk_size),    ptch(alloc_msk_size), &
                                                            stat=ierr )
  !!!        if (ierr == 0) then
  !!!C$OMP CRITICAL(etch_cr_rgn)
  !!!          write(*,'(1x,2A,F16.8,1x,A,I3)') 'etch_into_land_thread :: ',
  !!!     &               'allocated',  dble(4*alloc_msk_size)/dble(262144),
  !!!     &                             'MB private workspace, trd =', trd
  !!!C$OMP END CRITICAL(etch_cr_rgn)
  !!!        else
  !!!          write(*,*) '### ERROR: etch_into_land_thread :: ',
  !!!     &                           'memory allocation error.'
  !!!        endif
        endif
  
  
        nmsk=0 !<-- number of masked points
        do j=j0,ny,jskip
          jstr=max(1,j) ; jend=min(j+jsize-1,ny)
          call init_ijmsk_tile(istr,iend,jstr,jend, nx,ny, qext, alloc_msk_size, ijmsk, nmsk)
        enddo
  
        do while (mskd_pts > 0)
          ncst=0 !<-- number of coastal points
          call select_coast_pts(nx,ny, qext, nmsk, ijmsk, ncst,nijmsk)
  
          npths=0 !<-- number of points to be patched
          call comp_patch_pts(nx,ny, qext, nmsk,ijmsk, ncst,nijmsk, npths, ijptch,ptch)
  !
          call apply_patch_pts(npths, ijptch,ptch, nx,ny, qext)
          do j=j0,ny,jskip
            jstr=max(1,j) ; jend=min(j+jsize-1,ny)
            call set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
          enddo
          nmsk_new=0
          call shortlist(nmsk, ijmsk, nmsk_new)
          nmsk=nmsk_new
  !        
          if (t_count == 0) mskd_pts=0
          t_count=t_count+1
          mskd_pts=mskd_pts+nmsk
          if (t_count == ntrds) then
            t_count=0
  !!!#define VERBOSE
  !!!#ifdef VERBOSE
  !!!          write(*,*) 'rmaining masked points =', mskd_pts
  !!!#endif
          endif
        enddo  !<-- while
      endif

      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call strip_halo_tile(istr,iend,jstr,jend, nx,ny, qext,qfld)
      enddo
      end

      subroutine copy_extend_tile(istr,iend,jstr,jend, nx,ny, mask, qsrc,qext, nmsk)

! Copy array "qsrc" into "qext" with has one row of ghost points all
! around, while setting land points to special values. The logic here
! is designed to work both ways: either there a non-trivial land mask
! array, or masked data is already set to some special value, while
! mask(i,j) is trivially set to all-water mask(:,:)=1 status and makes
! no effect. The secondary goal is to determine the total number of
! special-value points encountered by this thread so an appropiate
! sized arrays can be allocated to hold list of indices.

      implicit none
      integer istr,iend,jstr,jend, nmsk, i,j
      integer nx,ny
      integer, dimension(nx,ny) :: mask
      real, dimension(nx,ny) :: qsrc
      real, dimension(0:nx+1,0:ny+1) :: qext
!>    write(*,*) 'enter copy_extend_tile'
      do j=jstr,jend
        do i=istr,iend
          if (mask(i,j) > 0 .and. abs(qsrc(i,j)) < abs(spv)) then
            qext(i,j)=qsrc(i,j)
          else
            qext(i,j)=spv_set ; nmsk=nmsk+1
          endif
        enddo
      enddo
      end

      subroutine strip_halo_tile(istr,iend,jstr,jend, nx, ny, qext,qsrc)
      implicit none
      integer istr,iend,jstr,jend, nx,ny, i,j
      real, dimension(0:nx+1,0:ny+1) :: qext
      real, dimension(nx,ny)         :: qsrc
!>    write(*,*) 'enter strip_halo_tile'
      do j=jstr,jend
        do i=istr,iend
          qsrc(i,j)=qext(i,j)
        enddo
      enddo
      end

      subroutine set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
      implicit none
      integer istr,iend,jstr,jend, i,j
      integer nx,ny
      real, dimension(0:nx+1,0:ny+1) :: qext
!>    write(*,*) 'enter set_qext_bc_tile'
      if (WESTERN_EDGE) then
        do j=jstr,jend
          qext(istr-1,j)=qext(istr,j)
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          qext(iend+1,j)=qext(iend,j)
        enddo
      endif
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          qext(i,jstr-1)=qext(i,jstr)
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          qext(i,jend+1)=qext(i,jend)
        enddo
      endif
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        qext(istr-1,jstr-1)=qext(istr,jstr)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        qext(istr-1,jend+1)=qext(istr,jend)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        qext(iend+1,jstr-1)=qext(iend,jstr)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        qext(iend+1,jend+1)=qext(iend,jend)
      endif
      end

      subroutine init_ijmsk_tile(istr,iend,jstr,jend, nx,ny, src, max_pts, ijmsk, nmsk)
      implicit none
      integer istr,iend,jstr,jend, nx,ny, max_pts, nmsk, i,j
      real, dimension(0:nx+1,0:ny+1)    :: src
      integer, dimension(2,max_pts) :: ijmsk
!>    write(*,*) 'enter init_ijmsk_tile'
      do j=jstr,jend                       ! Form list of indices of
        do i=istr,iend                     ! points with special value.
          if (src(i,j) < spv) then
            nmsk=nmsk+1
            ijmsk(1,nmsk)=i
            ijmsk(2,nmsk)=j
          endif
        enddo
      enddo
      end

      subroutine select_coast_pts(nx,ny, src, nmsk,ijmsk, ncst,nijmsk)
      implicit none
      integer nx,ny, nmsk, ncst, i,j,n     ! Take list of previously
      real, dimension(0:nx+1,0:ny+1)    :: src      ! identified masked points
      integer, dimension(2,nmsk) :: ijmsk        ! and select the ones among
      integer, dimension(nmsk)   :: nijmsk         ! them which have at least
!>    write(*,*) 'enter select_coast_pts'  ! one water neighbor.
      do n=1,nmsk
        i=ijmsk(1,n) ; j=ijmsk(2,n)
        if (src(i,j) < spv) then
          if (src(i+1,j) > spv .or. src(i,j+1) > spv .or. &
              src(i-1,j) > spv .or. src(i,j-1) > spv ) then
            ncst=ncst+1 ; nijmsk(ncst)=n
          endif
        endif
      enddo
      end

      subroutine comp_patch_pts(nx,ny, src, nmsk,ijmsk, ncst, nijmsk, npths, ijptch, ptch)
      implicit none
      integer nx, ny, nmsk, ncst, npths, i,j,n
      real, dimension(0:nx+1,0:ny+1)    :: src
      real, dimension(ncst)      :: ptch
      real                    :: wgt,vlu
      integer, dimension(2,nmsk) :: ijmsk
      integer, dimension(2,ncst) :: ijptch
      integer, dimension(nmsk)   :: nijmsk

      real, parameter :: grad=1./3.,  corn=0.707106781186547, corngrad=0.5*corn*grad, threshold=2.4

!>    write(*,*) 'enter comp_patch_pts, ncst =',ncst

      do n=1,ncst                         ! check surrounding points:
        i=ijmsk(1,nijmsk(n)) ; wgt=0.     ! counterclockwise direction
        j=ijmsk(2,nijmsk(n)) ; vlu=0.     ! starting from the east.

        if (src(i+1,j) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i+1,j)
          if (i < nx) then
            if (src(i+2,j) > spv) vlu=vlu+grad*(src(i+1,j)-src(i+2,j))
          endif
        endif
#  define DIAGONAL
#  ifdef DIAGONAL
        if (src(i+1,j+1) > spv) then
          wgt=wgt+corn ;  vlu=vlu + corn*src(i+1,j+1)
          if (i < nx) then
            if (src(i+2,j+1) > spv) vlu=vlu + corngrad*( src(i+1,j+1)-src(i+2,j+1) )
          endif
          if (j < ny) then
            if (src(i+1,j+2) > spv) vlu=vlu  + corngrad*(src(i+1,j+1)-src(i+1,j+2) )
          endif
        endif
#  endif

        if (src(i,j+1) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i,j+1)
          if (j < ny) then
            if (src(i,j+2) > spv) vlu=vlu+grad*(src(i,j+1)-src(i,j+2))
          endif
        endif

#  ifdef DIAGONAL
        if (src(i-1,j+1) > spv) then
          wgt=wgt+corn ;   vlu=vlu + corn*src(i-1,j+1)
          if (j < ny) then
            if (src(i-1,j+2) > spv) vlu=vlu  + corngrad*(src(i-1,j+1)-src(i-1,j+2))
          endif
          if (i > 1) then
            if (src(i-2,j+1) > spv) vlu=vlu  + corngrad*(src(i-1,j+1)-src(i-2,j+1))
          endif
        endif
#  endif

        if (src(i-1,j) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i-1,j)
          if (i > 1) then
            if (src(i-2,j) > spv) vlu=vlu+grad*(src(i-1,j)-src(i-2,j))
          endif
        endif

#  ifdef DIAGONAL
        if (src(i-1,j-1) > spv) then
          wgt=wgt+corn ;   vlu=vlu + corn*src(i-1,j-1)
          if (i > 1) then
            if (src(i-2,j-1) > spv) vlu=vlu  + corngrad*(src(i-1,j-1)-src(i-2,j-1))
          endif
          if (j > 1) then
            if (src(i-1,j-2) > spv) vlu=vlu  + corngrad*(src(i-1,j-1)-src(i-1,j-2))
          endif
        endif
#  endif

        if (src(i,j-1) > spv) then
          wgt=wgt+1.  ;  vlu=vlu + src(i,j-1)
          if (j > 1) then
            if (src(i,j-2) > spv) vlu=vlu+grad*(src(i,j-1)-src(i,j-2))
          endif
        endif

#  ifdef DIAGONAL
        if (src(i+1,j-1) > spv) then
          wgt=wgt+corn ;  vlu=vlu + corn*src(i+1,j-1)
          if (j > 1) then
            if (src(i+1,j-2) > spv) vlu=vlu  + corngrad*(src(i+1,j-1)-src(i+1,j-2))
          endif
          if (i < nx) then
            if (src(i+2,j-1) > spv) vlu=vlu  + corngrad*(src(i+1,j-1)-src(i+2,j-1))
          endif
        endif
#  endif
        if (wgt > threshold) then
          npths=npths+1              ! At the end set "ijmsk" i-index
          ptch(npths)=vlu/wgt        ! to zero to signal that the point
          ijptch(1,npths)=i          ! is no longer a special value.
          ijptch(2,npths)=j
          ijmsk(1,nijmsk(n))=0
        endif
      enddo
      end

      subroutine apply_patch_pts(npths, ijptch,ptch, nx,ny, src)
      implicit none
      integer npths, nx,ny, i,j,n
      integer, dimension(2,npths) :: ijptch
      real, dimension(npths)      :: ptch
      real, dimension(0:nx+1,0:ny+1)    :: src

!>    write(*,*) 'enter apply_patch_pts'
      do n=1,npths
        i=ijptch(1,n) ; j=ijptch(2,n)
        src(i,j)=ptch(n)
      enddo
      end


      subroutine shortlist(nmsk, ijmsk, nmsk_new)
      implicit none
      integer nmsk,nmsk_new,n
      integer, dimension(2,nmsk) :: ijmsk

      do n=1,nmsk
        if (ijmsk(1,n) > 0) then
          nmsk_new=nmsk_new+1
          ijmsk(1,nmsk_new)=ijmsk(1,n)
          ijmsk(2,nmsk_new)=ijmsk(2,n)
        endif
      enddo
!>    write(*,*) 'shortlist, nmsk=', nmsk, ' -->', nmsk_new
      end
# endif 
!
! ------------------------------------------------------------------------------------------
      subroutine check(status,msg)
        integer, intent ( in) :: status
        character(len=*),optional :: msg
      
        if(status /= nf90_noerr) then
          MPI_master_only write(stdout,*) trim(nf90_strerror(status))
          if(present(msg)) then
            MPI_master_only write(stdout,*) msg
          endif  
        end if
      end subroutine check
!
!-------------------------------------------------------------------------------------------
      subroutine dealloc_closecdf
        integer :: nset, i, ierr

        deallocate(obj)
        deallocate(unique_set)

        nset = size(unique_set)
        do i=1,nset
          if(ncidextr(i).ne.-1) then
            ierr = nf90_close(ncidextr(i))
          endif
        enddo
        deallocate(ncidextr)
        deallocate(extrTime, extrZ_N , extrZ_S ,extrZ_W , extrZ_E , &
                             extrUb_N, extrUb_S,extrUb_W, extrUb_E, &
                             extrVb_N, extrVb_S,extrVb_W, extrVb_E, &
                             extrU_N , extrU_S ,extrU_W , extrU_E , &
                             extrV_N , extrV_S ,extrV_W , extrV_E , &
                             extrT_N , extrT_S ,extrT_W , extrT_E , &
                             extrS_N , extrS_S ,extrS_W , extrS_E  )
 
        end subroutine
!
! ------------------------------------------------------------------------------------------
      logical function findstr(string,pattern,istart) ![
      ! Should be in roms_read_write and then that module should be available to
      ! Tools-Roms/ as well as src/
      ! didn't do it now because roms_read_write depends on other modules so would
      ! require a bit of work (do once ncvars is removed)
      implicit none

      !input/output
      character(len=*),intent(in)  :: string        ! string
      character(len=*),intent(in)  :: pattern       ! desired pattern to find within string
      integer,optional,intent(out) :: istart

      !local
      integer :: nl,nlv,i

      nl  = len(trim(string))
      nlv = len(pattern)

      findstr = .false.
      do i = 1,nl-nlv+1
         if (string(i:i+nlv-1) == pattern) then
          findstr = .true.
          exit
         endif
      enddo

      if (present(istart)) then
        if (findstr) then
          istart=i                                  ! return string starting index
        else
          istart=0
        endif
      endif

      end function findstr !]

#endif

      end module online_extract_bdy

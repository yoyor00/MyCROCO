#ifdef CROCO_TIMINGS
!======================================================================
!
!  Simplified wall-clock timing facility for CROCO.
!
!  The module is deliberately small and self contained: it only uses
!  Fortran 2003 features, has no C helper code and does not depend on
!  any external library.  Regions are created on demand, so the number
!  of regions is not limited by a compile-time constant.
!
!  Annotations in the source files are comments:
!
!     !$croco timings use
!     !$croco timings init
!     !$croco timings start my_region
!     !$croco timings stop  my_region
!     !$croco timings finalize
!
!  They are converted to calls to this module by timings_convert.py.
!  Without -DCROCO_TIMINGS the complete module is compiled away.
!
!======================================================================
module timings
   implicit none
   private

   public :: timings_init
   public :: timings_start
   public :: timings_stop
   public :: timings_finalize

   ! Named stops use (name, id).  The id-first form is kept as a
   ! compatibility path for existing callers.
   interface timings_stop
      module procedure timings_stop_by_name
      module procedure timings_stop_by_id
   end interface

   ! One row of the growable region table.
   type :: timing_region
      character(len=:), allocatable :: name
      logical :: defined = .false.
      logical :: root = .false.
      integer(8) :: calls = 0
      real(8) :: total_seconds = 0.0d0
      real(8) :: min_seconds = huge(1.0d0)
      real(8) :: max_seconds = 0.0d0
   end type timing_region

   ! Error records are collected rather than silently ignored.
   type :: timing_error
      character(len=:), allocatable :: kind
      character(len=:), allocatable :: region
      character(len=:), allocatable :: message
   end type timing_error

   type(timing_region), allocatable :: regions(:)
   type(timing_error), allocatable :: errors(:)
   integer :: nregions = 0
   integer :: region_capacity = 0
   integer :: nerrors = 0
   integer :: error_capacity = 0

   ! Stack of currently open regions (innermost last).
   integer, allocatable :: stack_id(:)
   real(8), allocatable :: stack_t0(:)
   integer :: stack_depth = 0
   integer :: stack_capacity = 0

   logical :: initialized = .false.
   logical :: finished = .false.
   character(len=:), allocatable :: output_base
   character(len=64) :: first_start = ""

contains

!======================================================================
!  Initialisation
!======================================================================
   subroutine timings_init()
      implicit none

      character(len=:), allocatable :: value
      logical :: found

      if (initialized) return

      call get_env_value('CROCO_TIMINGS_FILE', value, found)
      if (found .and. len(value) > 0) then
         output_base = value
      else
         output_base = 'timings'
      end if
      if (allocated(value)) deallocate (value)

      call timestamp_now(first_start)

      initialized = .true.
      finished = .false.
   end subroutine timings_init

!======================================================================
!  Start a region.  The name is registered the first time it is seen;
!  later starts of the same name reuse the existing entry.  The optional
!  id selects an explicit region id; -1 (or no argument) uses the name.
!======================================================================
   subroutine timings_start(name, id)
      implicit none
      character(len=*), intent(in) :: name
      integer, intent(in), optional :: id

      integer :: region_id, requested_id
      real(8) :: t

      if (.not. initialized) call timings_init()
      if (finished) return

      requested_id = -1
      if (present(id)) requested_id = id

      call clock_seconds(t)
      region_id = find_or_register(name, requested_id)

      call push_region(region_id, t)
   end subroutine timings_start

!======================================================================
!  Stop a region.  The name is the first argument and the optional id the
!  second.  With id == -1 (or no id) the matching open region is searched
!  by name, starting at the innermost level.  A non-negative id is
!  matched directly against the open stack.
!======================================================================
   subroutine timings_stop_by_name(name, id)
      implicit none
      character(len=*), intent(in) :: name
      integer, intent(in), optional :: id

      integer :: requested_id

      requested_id = -1
      if (present(id)) requested_id = id

      call stop_region(requested_id, name)
   end subroutine timings_stop_by_name

!======================================================================
!  Compatibility form: timings_stop(id [, name]).
!======================================================================
   subroutine timings_stop_by_id(region_id, region_name)
      implicit none
      integer, intent(in) :: region_id
      character(len=*), intent(in), optional :: region_name

      if (present(region_name)) then
         call stop_region(region_id, region_name)
      else
         call stop_region(region_id, '')
      end if
   end subroutine timings_stop_by_id

   subroutine stop_region(id, name)
      implicit none
      integer, intent(in) :: id
      character(len=*), intent(in) :: name

      integer :: k, p, target
      real(8) :: t

      if (.not. initialized) call timings_init()
      if (finished) return

      call clock_seconds(t)

      if (id >= 0) then
         target = id
      else if (len_trim(name) > 0) then
         target = find_open_by_name(name)
         if (target < 0) then
            call add_error('unbalanced_stop', name, &
     &           'stop without matching start')
            return
         end if
      else
         call add_error('invalid_stop', '', &
     &        'stop called without an id or a region name')
         return
      end if

      p = stack_depth
      do while (p >= 1)
         if (stack_id(p) == target) exit
         p = p - 1
      end do

      if (p < 1) then
         if (target >= 0 .and. target < nregions) then
            if (regions(target + 1)%defined) then
               call add_error('unbalanced_stop', regions(target + 1)%name, &
     &              'stop without matching start')
            else
               call add_error('unbalanced_stop', name, &
     &              'stop without matching start')
            end if
         else
            call add_error('unbalanced_stop', name, &
     &           'stop without matching start')
         end if
         return
      end if

      ! Regions opened after the target are closed implicitly, but the
      ! imbalance is reported in the JSON output.
      do k = stack_depth, p + 1, -1
         call add_error('unbalanced_stop', regions(stack_id(k) + 1)%name, &
     &        'region closed out of order')
         call close_region(stack_id(k), stack_t0(k), t)
      end do
      call close_region(target, stack_t0(p), t)
      stack_depth = p - 1
   end subroutine stop_region

!======================================================================
!  Write the JSON file and a short human readable summary.
!======================================================================
   subroutine timings_finalize()
      implicit none

      integer :: i, unit, rank, size
      logical :: is_mpi
      character(len=:), allocatable :: filename
      real(8) :: root_total, percent, tnow

      if (.not. initialized .or. finished) return

      ! Any region still open is reported and closed here.
      call clock_seconds(tnow)
      do i = stack_depth, 1, -1
         call add_error('unclosed_region', regions(stack_id(i) + 1)%name, &
     &        'region still open at finalize')
         call close_region(stack_id(i), stack_t0(i), tnow)
      end do
      stack_depth = 0

      call detect_mpi(rank, size, is_mpi)
      call make_filename(output_base, rank, is_mpi, filename)

      open (newunit=unit, file=filename, status='replace', action='write', &
     &     iostat=i)
      if (i /= 0) then
         write (*, '(a,a)') '[timings] cannot write ', filename
      else
         call write_json(unit, rank, size, is_mpi)
         close (unit)
         write (*, '(a,a)') '[timings] wrote ', filename
      end if

      root_total = 0.0d0
      do i = 1, nregions
         if (.not. regions(i)%defined) cycle
         if (regions(i)%root) root_total = root_total + regions(i)%total_seconds
      end do

      do i = 1, nregions
         if (.not. regions(i)%defined) cycle
         percent = 0.0d0
         if (root_total > 0.0d0) then
            percent = 100.0d0*regions(i)%total_seconds/root_total
         end if
         write (*, '(a,1x,a,1x,i0,1x,es12.5,1x,a,1x,f6.2,a)') &
     &        '[timings]', trim(regions(i)%name), regions(i)%calls, &
     &        regions(i)%total_seconds, 's', percent, '%'
      end do

      if (nerrors > 0) then
         write (*, '(a,i0,a)') '[timings] ', nerrors, ' error(s) recorded'
      end if

      finished = .true.
   end subroutine timings_finalize

!======================================================================
!  Region table handling
!======================================================================
   function find_or_register(name, requested_id) result(id)
      implicit none
      character(len=*), intent(in) :: name
      integer, intent(in) :: requested_id
      integer :: id

      if (requested_id >= 0) then
         call ensure_region_capacity(requested_id + 1)
         if (.not. regions(requested_id + 1)%defined) then
            id = requested_id
            call define_region(id, name)
            return
         end if
         if (same_name(regions(requested_id + 1)%name, name)) then
            id = requested_id
            return
         end if

         call add_error('id_conflict', name, &
     &        'id is already used by another region; using name lookup')
      end if

      id = find_registered_by_name(name)
      if (id >= 0) return

      call ensure_region_capacity(nregions + 1)
      id = nregions
      call define_region(id, name)
   end function find_or_register

   subroutine define_region(id, name)
      implicit none
      integer, intent(in) :: id
      character(len=*), intent(in) :: name

      integer :: n

      n = len_trim(name)
      allocate (character(len=n) :: regions(id + 1)%name)
      regions(id + 1)%name = name(1:n)
      regions(id + 1)%defined = .true.
      regions(id + 1)%root = .false.
      regions(id + 1)%calls = 0
      regions(id + 1)%total_seconds = 0.0d0
      regions(id + 1)%min_seconds = huge(1.0d0)
      regions(id + 1)%max_seconds = 0.0d0
      if (id + 1 > nregions) nregions = id + 1
   end subroutine define_region

   function find_registered_by_name(name) result(id)
      implicit none
      character(len=*), intent(in) :: name
      integer :: id

      integer :: i

      id = -1
      do i = 1, nregions
         if (.not. regions(i)%defined) cycle
         if (same_name(regions(i)%name, name)) then
            id = i - 1
            return
         end if
      end do
   end function find_registered_by_name

   function same_name(stored, requested) result(equal)
      implicit none
      character(len=*), intent(in) :: stored, requested
      logical :: equal

      equal = (trim(stored) == trim(requested))
   end function same_name

   subroutine ensure_region_capacity(needed)
      implicit none
      integer, intent(in) :: needed

      type(timing_region), allocatable :: tmp(:)
      integer :: newcap

      if (needed <= region_capacity) return

      newcap = max(8, region_capacity)
      do while (newcap < needed)
         newcap = 2*newcap
      end do

      allocate (tmp(newcap))
      if (nregions > 0) tmp(1:nregions) = regions(1:nregions)
      call move_alloc(tmp, regions)
      region_capacity = newcap
   end subroutine ensure_region_capacity

!======================================================================
!  Open-region stack
!======================================================================
   subroutine push_region(id, t)
      implicit none
      integer, intent(in) :: id
      real(8), intent(in) :: t

      integer, allocatable :: itmp(:)
      real(8), allocatable :: rtmp(:)
      integer :: newcap

      if (stack_depth + 1 > stack_capacity) then
         newcap = max(8, 2*stack_capacity)
         allocate (itmp(newcap))
         allocate (rtmp(newcap))
         if (stack_depth > 0) then
            itmp(1:stack_depth) = stack_id(1:stack_depth)
            rtmp(1:stack_depth) = stack_t0(1:stack_depth)
         end if
         call move_alloc(itmp, stack_id)
         call move_alloc(rtmp, stack_t0)
         stack_capacity = newcap
      end if

      stack_depth = stack_depth + 1
      stack_id(stack_depth) = id
      stack_t0(stack_depth) = t
      if (stack_depth == 1) regions(id + 1)%root = .true.
   end subroutine push_region

   function find_open_by_name(name) result(id)
      implicit none
      character(len=*), intent(in) :: name
      integer :: id

      integer :: i, n

      n = len_trim(name)
      do i = stack_depth, 1, -1
         if (len(regions(stack_id(i) + 1)%name) == n) then
            if (regions(stack_id(i) + 1)%name == name(1:n)) then
               id = stack_id(i)
               return
            end if
         end if
      end do
      id = -1
   end function find_open_by_name

   subroutine close_region(id, t0, t1)
      implicit none
      integer, intent(in) :: id
      real(8), intent(in) :: t0, t1

      real(8) :: elapsed

      elapsed = t1 - t0
      if (elapsed < 0.0d0) elapsed = 0.0d0

      regions(id + 1)%calls = regions(id + 1)%calls + 1
      regions(id + 1)%total_seconds = regions(id + 1)%total_seconds + elapsed
      if (elapsed < regions(id + 1)%min_seconds) then
         regions(id + 1)%min_seconds = elapsed
      end if
      if (elapsed > regions(id + 1)%max_seconds) then
         regions(id + 1)%max_seconds = elapsed
      end if
   end subroutine close_region

!======================================================================
!  Error table
!======================================================================
   subroutine add_error(kind, region, message)
      implicit none
      character(len=*), intent(in) :: kind, region, message

      type(timing_error), allocatable :: tmp(:)
      integer :: newcap

      if (nerrors + 1 > error_capacity) then
         newcap = max(8, 2*error_capacity)
         allocate (tmp(newcap))
         if (nerrors > 0) tmp(1:nerrors) = errors(1:nerrors)
         call move_alloc(tmp, errors)
         error_capacity = newcap
      end if

      nerrors = nerrors + 1
      allocate (character(len=len(kind)) :: errors(nerrors)%kind)
      allocate (character(len=len(region)) :: errors(nerrors)%region)
      allocate (character(len=len(message)) :: errors(nerrors)%message)
      errors(nerrors)%kind = kind
      errors(nerrors)%region = region
      errors(nerrors)%message = message
   end subroutine add_error

!======================================================================
!  Clock and environment helpers
!======================================================================
   subroutine clock_seconds(seconds)
      implicit none
      real(8), intent(out) :: seconds

      integer :: count, rate

      call system_clock(count, rate)
      seconds = real(count, 8)/real(rate, 8)
   end subroutine clock_seconds

   subroutine timestamp_now(stamp)
      implicit none
      character(len=*), intent(out) :: stamp

      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone

      call date_and_time(date, time, zone)
      stamp = date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'T'// &
     &        time(1:2)//':'//time(3:4)//':'//time(5:6)// &
     &        trim(zone)
   end subroutine timestamp_now

   subroutine get_env_value(name, value, found)
      implicit none
      character(len=*), intent(in) :: name
      character(len=:), allocatable, intent(out) :: value
      logical, intent(out) :: found

      character(len=2048) :: buffer
      integer :: length, status

      length = 0
      call get_environment_variable(name, buffer, length, status)
      if (status == 0 .and. length > 0) then
         allocate (character(len=length) :: value)
         value = buffer(1:length)
         found = .true.
      else
         allocate (character(len=0) :: value)
         found = .false.
      end if
   end subroutine get_env_value

   subroutine parse_integer(text, value, found)
      implicit none
      character(len=*), intent(in) :: text
      integer, intent(out) :: value
      logical, intent(out) :: found

      integer :: ios

      value = 0
      found = .false.
      if (len_trim(text) == 0) return
      read (text, *, iostat=ios) value
      if (ios == 0) found = .true.
   end subroutine parse_integer

   subroutine detect_mpi(rank, size, is_mpi)
      implicit none
      integer, intent(out) :: rank, size
      logical, intent(out) :: is_mpi

      character(len=:), allocatable :: value
      logical :: found

      rank = 0
      size = 1
      is_mpi = .false.

      call get_env_value('OMPI_COMM_WORLD_RANK', value, found)
      if (found) then
         is_mpi = .true.
         call parse_integer(value, rank, found)
      end if
      if (allocated(value)) deallocate (value)

      call get_env_value('OMPI_COMM_WORLD_SIZE', value, found)
      if (found) then
         is_mpi = .true.
         call parse_integer(value, size, found)
      end if
      if (allocated(value)) deallocate (value)

      call read_mpi_env('PMI_RANK', rank, found)
      if (found) is_mpi = .true.
      call read_mpi_env('SLURM_PROCID', rank, found)
      if (found) is_mpi = .true.
      call read_mpi_env('MPI_LOCALRANKID', rank, found)
      if (found) is_mpi = .true.

      call read_mpi_env('PMI_SIZE', size, found)
      if (found) is_mpi = .true.
      call read_mpi_env('SLURM_NTASKS', size, found)
      if (found) is_mpi = .true.
      call read_mpi_env('MPI_LOCALNRANKS', size, found)
      if (found) is_mpi = .true.
   end subroutine detect_mpi

   subroutine read_mpi_env(name, value, found)
      implicit none
      character(len=*), intent(in) :: name
      integer, intent(inout) :: value
      logical, intent(out) :: found

      character(len=:), allocatable :: text
      logical :: present
      integer :: parsed

      call get_env_value(name, text, present)
      found = .false.
      if (present) then
         call parse_integer(text, parsed, found)
         if (found) value = parsed
      end if
      if (allocated(text)) deallocate (text)
   end subroutine read_mpi_env

   subroutine make_filename(base, rank, is_mpi, filename)
      implicit none
      character(len=*), intent(in) :: base
      integer, intent(in) :: rank
      logical, intent(in) :: is_mpi
      character(len=:), allocatable, intent(out) :: filename

      character(len=32) :: suffix
      integer :: n

      n = len_trim(base)
      if (is_mpi) then
         write (suffix, '(a,i4.4,a)') '_mpi_rank_', rank, '.json'
         if (n >= 5) then
            if (base(n - 4:n) == '.json') then
               filename = base(1:n - 5)//trim(suffix)
               return
            end if
         end if
         filename = base(1:n)//trim(suffix)
      else
         if (n >= 5) then
            if (base(n - 4:n) == '.json') then
               filename = base(1:n)
               return
            end if
         end if
         filename = base(1:n)//'.json'
      end if
   end subroutine make_filename

!======================================================================
!  JSON writer
!======================================================================
   subroutine write_json(unit, rank, size, is_mpi)
      implicit none
      integer, intent(in) :: unit, rank, size
      logical, intent(in) :: is_mpi

      integer :: i, ndefined, nwritten
      logical :: complete
      real(8) :: root_total, percent

      ndefined = 0
      root_total = 0.0d0
      do i = 1, nregions
         if (.not. regions(i)%defined) cycle
         ndefined = ndefined + 1
         if (regions(i)%root) root_total = root_total + regions(i)%total_seconds
      end do
      complete = (nerrors == 0)

      write (unit, '(a)') '{'
      write (unit, '(a)') '  "program": "croco",'
      write (unit, '(a,i0,a)') '  "mpi_rank": ', rank, ','
      write (unit, '(a,i0,a)') '  "mpi_size": ', size, ','
      write (unit, '(a,a,a)') '  "mpi": ', json_bool(is_mpi), ','
      write (unit, '(a)') '  "clock": "system_clock",'
      write (unit, '(a,a,a)') '  "first_start": ', json_quote(trim(first_start)), ','
      write (unit, '(a,a,a)') '  "complete": ', json_bool(complete), ','

      write (unit, '(a)') '  "regions": ['
      nwritten = 0
      do i = 1, nregions
         if (.not. regions(i)%defined) cycle
         nwritten = nwritten + 1
         percent = 0.0d0
         if (root_total > 0.0d0) then
            percent = 100.0d0*regions(i)%total_seconds/root_total
         end if
         write (unit, '(a)', advance='no') '    {"id": '
         write (unit, '(i0)', advance='no') i - 1
         write (unit, '(a)', advance='no') ', "name": '
         write (unit, '(a)', advance='no') json_quote(regions(i)%name)
         write (unit, '(a)', advance='no') ', "region": '
         write (unit, '(a)', advance='no') json_quote(regions(i)%name)
         write (unit, '(a,i0)', advance='no') ', "calls": ', regions(i)%calls
         write (unit, '(a,i0)', advance='no') ', "counter": ', regions(i)%calls
         write (unit, '(a)', advance='no') ', "total_seconds": '
         write (unit, '(es23.15e3)', advance='no') regions(i)%total_seconds
         write (unit, '(a)', advance='no') ', "min_seconds": '
         if (regions(i)%calls > 0) then
            write (unit, '(es23.15e3)', advance='no') regions(i)%min_seconds
         else
            write (unit, '(es23.15e3)', advance='no') 0.0d0
         end if
         write (unit, '(a)', advance='no') ', "max_seconds": '
         write (unit, '(es23.15e3)', advance='no') regions(i)%max_seconds
         write (unit, '(a)', advance='no') ', "total_percent": '
         write (unit, '(f8.3)', advance='no') percent
         if (nwritten < ndefined) then
            write (unit, '(a)') '},'
         else
            write (unit, '(a)') '}'
         end if
      end do
      write (unit, '(a)') '  ],'

      write (unit, '(a)') '  "errors": ['
      do i = 1, nerrors
         write (unit, '(a)', advance='no') '    {"kind": '
         write (unit, '(a)', advance='no') json_quote(errors(i)%kind)
         write (unit, '(a)', advance='no') ', "region": '
         write (unit, '(a)', advance='no') json_quote(errors(i)%region)
         write (unit, '(a)', advance='no') ', "message": '
         write (unit, '(a)', advance='no') json_quote(errors(i)%message)
         if (i < nerrors) then
            write (unit, '(a)') '},'
         else
            write (unit, '(a)') '}'
         end if
      end do
      write (unit, '(a)') '  ]'
      write (unit, '(a)') '}'
   end subroutine write_json

   function json_bool(value) result(text)
      implicit none
      logical, intent(in) :: value
      character(len=:), allocatable :: text

      if (value) then
         text = 'true'
      else
         text = 'false'
      end if
   end function json_bool

   function json_quote(text) result(quoted)
      implicit none
      character(len=*), intent(in) :: text
      character(len=:), allocatable :: quoted

      character(len=:), allocatable :: tmp
      integer :: i, k

      allocate (character(len=2*len(text) + 2) :: quoted)
      quoted(1:1) = '"'
      k = 1
      do i = 1, len(text)
         select case (text(i:i))
         case ('"')
            k = k + 1; quoted(k:k) = '\'
            k = k + 1; quoted(k:k) = '"'
         case ('\')
            k = k + 1; quoted(k:k) = '\'
            k = k + 1; quoted(k:k) = '\'
         case (char(10))
            k = k + 1; quoted(k:k) = '\'
            k = k + 1; quoted(k:k) = 'n'
         case (char(13))
            k = k + 1; quoted(k:k) = '\'
            k = k + 1; quoted(k:k) = 'r'
         case (char(9))
            k = k + 1; quoted(k:k) = '\'
            k = k + 1; quoted(k:k) = 't'
         case default
            if (iachar(text(i:i)) < 32) then
               k = k + 1; quoted(k:k) = ' '
            else
               k = k + 1; quoted(k:k) = text(i:i)
            end if
         end select
      end do
      k = k + 1
      quoted(k:k) = '"'
      tmp = quoted(1:k)
      call move_alloc(tmp, quoted)
   end function json_quote

end module timings
#endif

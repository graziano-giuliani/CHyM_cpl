module mod_chym_iface

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

  use, intrinsic :: iso_fortran_env
  use mpi
  use mod_chym_param
  use mod_chym_io
  use mod_chym_model

  implicit none

  logical :: model_standalone

  private

!-----------------------------------------------------------------------
!     Public subroutines
!-----------------------------------------------------------------------

  public :: chym_init
  public :: chym_start
  public :: chym_run
  public :: chym_finalize

  contains

    subroutine chym_init(mpiCommunicator)
      implicit none
      integer, intent(in), optional :: mpiCommunicator
      integer :: iprov, ierr, irank, nproc
      character(len=256) :: namelistfile

!-----------------------------------------------------------------------
!     Read configuration parameters
!-----------------------------------------------------------------------

      write(output_unit, *) 'This is CHYM - coupled runoff version'
      write(output_unit, *) 'Version coupled in RegCM-ES-1.0'
      write(output_unit, *) 'CHYM - Reading namelist files'

      if (present(mpiCommunicator)) then
        call mpi_comm_dup(mpiCommunicator,mycomm,ierr)
        if ( ierr /= 0 ) then
          write(error_unit,*) 'Cannot get communicator!'
          call mpi_abort(mpiCommunicator,100,ierr)
        end if
        model_standalone = .false.
      else
        call mpi_init_thread(mpi_thread_single,iprov,ierr)
        if ( ierr /= mpi_success ) then
          write(error_unit,*) 'Cannot initilize MPI'
          stop
        end if
        call mpi_comm_dup(MPI_COMM_WORLD,mycomm,ierr)
        if ( ierr /= 0 ) then
          write(error_unit,*) 'Cannot get communicator!'
          call mpi_abort(mpiCommunicator,100,ierr)
        end if
        model_standalone = .true.
      end if

      call mpi_comm_rank(mycomm, irank, ierr)
      if ( ierr /= 0 ) then
        write(error_unit,*) 'mpi_comm_rank Failure!'
        call mpi_abort(mycomm,3,ierr)
      end if
      call mpi_comm_size(mycomm, nproc, ierr)
      if ( ierr /= 0 ) then
        write(error_unit,*) 'mpi_comm_size Failure!'
        call mpi_abort(mycomm,4,ierr)
      end if

      if ( nproc > 1 ) then
        write(error_unit,*) 'CHyM for coupling runs on a single MPI process'
        call mpi_abort(mycomm,5,ierr)
      end if

#ifdef CPL
      namelistfile = 'chym.namelist'
#else
      call getarg(1, namelistfile)
#endif
      call read_config(trim(namelistfile))

!-----------------------------------------------------------------------
!     Initialize variables
!     Read parameter, restart and masking files
!-----------------------------------------------------------------------

      write(output_unit, *) 'CHYM - Initialize variables'
      call read_init()

!-----------------------------------------------------------------------
!     Open output files
!-----------------------------------------------------------------------

      write(output_unit, *) 'CHYM - Initialize output files'
      if (iswrit /= 0) then
        call chym_out_init()
        call chym_rst_init()
      end if

    end subroutine chym_init

    subroutine chym_start(imon,iday)
      implicit none
      integer, intent(in) :: imon, iday

      write(output_unit, *) 'CHYM - Initializing data space'

      if ( isread == 1 ) then
        write(output_unit, *) 'CHYM - WARM RESTART.'
        call coupler(imon,iday)
      else
        write(output_unit, *) 'CHYM - COLD START.'
        chym_dis(:,:) = 0.0
      end if
      write(output_unit, *) 'CHYM - Ready for time integration !!!!'

    end subroutine chym_start

!-----------------------------------------------------------------------
!       Run the model
!-----------------------------------------------------------------------

    subroutine chym_run(istart, iend, restarted, imon, iday)
      implicit none
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      logical, intent(in) :: restarted
      integer, intent(in) :: imon , iday

      integer :: istep, idir, i, j, ii, jj

!-----------------------------------------------------------------------
!     Run the model
!-----------------------------------------------------------------------

      write(output_unit, *) 'CHYM - Time integration step'
      do istep = istart, iend

!-----------------------------------------------------------------------
!       Get input
!-----------------------------------------------------------------------

#ifndef CPL
        ! information comes from input file
        if (istep == istart .and. istart /= 1) then
          write(output_unit,*) "CHYM - Restarting the model ..."
        end if
#endif

        call chymmodel(chym_runoff, chym_dis, imon, iday)

        if (iswrit /= 0) then
          if (mod(istep,iswrit) == 0) then
            call chym_rst(istep)
            call chym_out(istep)
            write(output_unit,fmt='(1x,A,I8)') &
                    'CHYM - Out and restart data are written', istep
          end if
        end if
      end do
      write(output_unit, *) 'CHYM - Time integration done!!!'
    end subroutine chym_run

!-----------------------------------------------------------------------
!   Model finalize
!-----------------------------------------------------------------------

    subroutine chym_finalize()
      implicit none
      integer :: ierr
      call chym_dispose( )
      if ( model_standalone ) then
        call mpi_finalize(ierr)
      end if
      write(output_unit, *) 'CHyM model finalized. Goodbye!'
    end subroutine chym_finalize

end module mod_chym_iface

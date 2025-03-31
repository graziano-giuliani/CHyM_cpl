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
  private

!-----------------------------------------------------------------------
!     Public subroutines
!-----------------------------------------------------------------------

  public :: chym_init
  public :: chym_start
  public :: chym_run
  public :: chym_finalize

  contains

    subroutine chym_init
      implicit none

      character(len=256) :: namelistfile

!-----------------------------------------------------------------------
!     Read configuration parameters
!-----------------------------------------------------------------------

      write(output_unit, *) 'This is CHYM - coupled runoff version'
      write(output_unit, *) 'Version coupled in RegCM-ES-1.0'
      write(output_unit, *) 'CHYM - Reading namelist files'

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

!-----------------------------------------------------------------------
!     Allocate variables
!-----------------------------------------------------------------------

      write(output_unit, *) 'CHYM - Allocate data space'

      if (.not. allocated(chym_runoff)) &
           allocate(chym_runoff(chym_nlon,chym_nlat))
      if (.not. allocated(chym_dis))    &
           allocate(chym_dis(chym_nlon,chym_nlat))

      chym_runoff(:,:) = 0.0
      chym_dis(:,:) = 0.0

    end subroutine chym_init

    subroutine chym_start(imon,iday)
      implicit none
      integer, intent(in) :: imon, iday
      integer :: i,j,idir,ilnd

      write(output_unit, *) 'CHYM - Initializing data space'

      if ( isread == 1 ) then

        do j = 2 , chym_nlat-1
          do i = 2 , chym_nlon-1
            idir = fmap(i,j)
            if ( idir >= 1 .and. idir <= 8 ) then
              ilnd = luse(i+ir(idir),j+jr(idir))
              if ( chym_drai(i,j) > thrriv .and. ilnd == ocean ) then
                chym_dis(i,j) = port(i,j)
              end if
            end if
          end do
        end do

#ifdef AZOV
        chym_dis(ikerch,jkerch) = 0.0
        do j = 2 , chym_nlat-1
          do i = 2 , chym_nlon-1
            if ( chym_lat(i,j) > 45.00 .and. chym_lon(i,j) > 34.0 ) then
              if ( chym_dis(i,j) > 0.0 ) then
                chym_dis(ikerch,jkerch) = chym_dis(ikerch,jkerch) + &
                               chym_dis(i,j)
                chym_dis(i,j) = 0.0
              end if
            end if
          end do
        end do
#endif

#ifdef NILE
        chym_dis(idamietta,jdamietta) = 0.5 * mval(nile_fresh_flux,imon,iday)
        chym_dis(irosetta,jrosetta) = 0.5 * mval(nile_fresh_flux,imon,iday)
#endif

#ifdef BLACKSEA
        chym_dis(idardanelli,jdardanelli) = mval(bs_fresh_flux,imon,iday)
#endif
      end if
      write(output_unit, *) 'CHYM - Ready to start !!!!'
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
!       Fake runoff to test

        chym_runoff = 1.0e-7
#endif
        call chymmodel(chym_runoff, chym_dis, imon, iday)

        do j = 2 , nbc-1
          do i = 2 , nlc-1
            if ( chym_lsm(i,j) > 0.0 ) then
              idir = fmap(i,j)
              if ( idir == 0 ) then
                call find_nearest_ocean(i,j,ii,jj)
                chym_dis(ii,jj) = chym_dis(i,j)
              else
                chym_dis(i+ir(idir),j+jr(idir)) = chym_dis(i,j)
              end if
            end if
          end do
        end do

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
      call chym_dispose( )
      write(output_unit, *) 'CHyM model finalized. Goodbye!'
    end subroutine chym_finalize

end module mod_chym_iface

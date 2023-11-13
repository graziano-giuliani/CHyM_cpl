      module mod_chym_iface
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use mpi
      use mod_chym_param
      use mod_chym_io
      use mod_chym_model
!
      implicit none
      private
!
!-----------------------------------------------------------------------
!     Public subroutines
!-----------------------------------------------------------------------
!
      public :: chym_init
      public :: chym_run
      public :: chym_finalize
!
      contains
!
      subroutine chym_init(imon, iday)
      implicit none
      integer, intent(in) :: imon, iday
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: i,j,idir,ilnd
      integer :: conta
!
!-----------------------------------------------------------------------
!     Read configuration parameters
!-----------------------------------------------------------------------
!
      call read_config('chymini.inp')
!
!-----------------------------------------------------------------------
!     Initialize variables
!     Read parameter, restart and masking files
!-----------------------------------------------------------------------
!
      call read_init()
!
!-----------------------------------------------------------------------
!     Open output files
!-----------------------------------------------------------------------
!
      call chym_out_init()
!
      if (iswrit /= 0) then
        call chym_rst_init()
      end if
!
!-----------------------------------------------------------------------
!     Allocate variables
!-----------------------------------------------------------------------
!
      if (.not. allocated(chym_runoff)) &
           allocate(chym_runoff(chym_nlon,chym_nlat))
      if (.not. allocated(chym_dis))    &
           allocate(chym_dis(chym_nlon,chym_nlat))
      if (.not. allocated(chym_surf))   &
           allocate(chym_surf(chym_nlon,chym_nlat))

      conta = 0

      if ( isread == 1 ) then
        chym_dis(:,:) = -1.0
        do j = 2 , chym_nlat-1
          do i = 2 , chym_nlon-1
            idir = int(fmap(i,j))
            if ( idir >= 1 .and. idir <= 8 ) then
              ilnd = int(luse(i+ir(idir),j+jr(idir)))
              if ( chym_drai(i,j) > thrriv .and. ilnd == mare ) then
                chym_dis(i,j) = port(i,j)
                conta = conta + 1
              end if
            end if
#ifdef NILE
            ! Booca 1 : Lat 31.52, Lon 31.84 (Damietta) 50%
            ! Bocca 2 : Lat 31.47, Lon 30.36 (Rosetta)  50%
            if ( sqrt(((lat1(j)-31.525)**2 + &
                       (lon1(i)-31.843)**2)) < 0.03 .or. &
                 sqrt(((lat1(j)-31.467)**2 + &
                       (lon1(i)-30.366)**2)) < 0.03) then
              chym_dis(i,j) = 0.5 * mval(nile_fresh_flux,imon,iday)
              conta = conta + 1
            end if
#endif
#ifdef BLACKSEA
            if ( sqrt(((lat1(j)-40.00)**2 + &
                       (lon1(i)-26.16)**2)) < 0.5*fscal ) then
              chym_dis(i,j) = mval(bs_fresh_flux,imon,iday)
              conta = conta + 1
            end if
#endif
          end do
        end do
        print*,"Number of discharge points accounted for : ",conta
      else
#ifdef CPL
        chym_dis(:,:) = -1.0
        do j = 2 , chym_nlat-1
          do i = 2 , chym_nlon-1
            idir = int(fmap(i,j))
            if ( idir >= 1 .and. idir <= 8 ) then
              ilnd = int(luse(i+ir(idir),j+jr(idir)))
              if ( chym_drai(i,j) > thrriv .and. ilnd == mare ) then
                chym_dis(i,j) = 1.0
                conta = conta + 1
              end if
            end if
#ifdef NILE
            ! Booca 1 : Lat 31.52, Lon 31.84 (Damietta) 50%
            ! Bocca 2 : Lat 31.47, Lon 30.36 (Rosetta)  50%
            if ( sqrt(((lat1(j)-31.525)**2 + &
                       (lon1(i)-31.843)**2)) < 0.03 .or. &
                 sqrt(((lat1(j)-31.467)**2 + &
                       (lon1(i)-30.366)**2)) < 0.03) then
              chym_dis(i,j) = 1.0
              conta = conta + 1
            end if
#endif
#ifdef BLACKSEA
            if ( sqrt(((lat1(j)-40.00)**2 + &
                       (lon1(i)-26.16)**2)) < 0.5*fscal ) then
              chym_dis(i,j) = 1.0
              conta = conta + 1
            end if
#endif
          end do
        end do
        print*,"Number of discharge points at beginning : ",conta
#endif
      end if
      end subroutine chym_init
!
      subroutine chym_run(istart, iend, restarted, imon, iday)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      logical, intent(in) :: restarted
      integer, intent(in) :: imon , iday
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: istep, icount
!
!-----------------------------------------------------------------------
!     Run the model
!-----------------------------------------------------------------------
!
      if (isread /= 0 .and. iswrit /= 0) then
        icount = mod(istart, iswrit)
      else
        icount = 0
      end if
!
      do istep = istart, iend
!
!-----------------------------------------------------------------------
!     Get input
!-----------------------------------------------------------------------
!
#ifdef CPL
      ! initial run
      if (istep == 1 .and. .not. restarted) then
        chym_runoff = 0.0
        chym_surf = 0.0
      end if
#else
      ! information comes from input file
      if (istep == istart .and. istart /= 1) then
        print*, "restarting the model ..."
      end if
#endif
!
!-----------------------------------------------------------------------
!     Run the model
!-----------------------------------------------------------------------
!
      call chymmodel(chym_runoff, chym_surf, chym_dis, imon, iday)
!
!-----------------------------------------------------------------------
!     Write to restart file
!-----------------------------------------------------------------------
!
      if (iswrit /= 0) then
        icount = icount+1
        if (iswrit == icount) then
          call chym_rst(istep)
          write(*,fmt='(A,I8)') 'restart data are written', istep
          icount = 0
        end if
      end if
!-----------------------------------------------------------------------
!     Write output to file
!-----------------------------------------------------------------------
!
      call chym_out(istep)
!
      end do
!
      end subroutine chym_run
!
!-----------------------------------------------------------------------
!     Model finalize
!-----------------------------------------------------------------------
!
      subroutine chym_finalize()
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     Close input files
!-----------------------------------------------------------------------
!
      write(*,fmt='(A)') 'CHyM model finalized '
!
      end subroutine chym_finalize
!
      end module mod_chym_iface

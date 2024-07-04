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
      real :: dlatlon
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
      if (.not. allocated(oro)) &
           allocate(oro(chym_nlon,chym_nlat))

      conta = 0
      chym_dis(:,:) = -1.0

      dlatlon = lon1(2)-lon1(1)

#ifdef AZOV
      ikerch = int((lon_kerch-lon1(1))/dlatlon)
      jkerch = int((lat_kerch-lat1(1))/dlatlon)
#endif

#ifdef NILE
      ! Booca 1 : Lat 31.52, Lon 31.84 (Damietta) 50%
      ! Bocca 2 : Lat 31.47, Lon 30.36 (Rosetta)  50%
      idamietta = int((lon_damietta-lon1(1))/dlatlon)
      jdamietta = int((lat_damietta-lat1(1))/dlatlon)
      irosetta = int((lon_rosetta-lon1(1))/dlatlon)
      jrosetta = int((lat_rosetta-lat1(1))/dlatlon)
#endif

#ifdef BLACKSEA
      idardanelli = int((lon_dardanelli-lon1(1))/dlatlon)
      jdardanelli = int((lat_dardanelli-lat1(1))/dlatlon)
#endif

      if ( isread == 1 ) then

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
          end do
        end do

#ifdef AZOV
        chym_dis(ikerch,jkerch) = 0.0
        conta = conta + 1
        do j = 2 , chym_nlat-1
          do i = 2 , chym_nlon-1
            if ( lat1(j) > 45.00 .and. lon1(i) > 34.0 ) then
              if ( chym_dis(i,j) > 0.0 ) then
                chym_dis(ikerch,jkerch) = chym_dis(ikerch,jkerch) + &
                               chym_dis(i,j)
                chym_dis(i,j) = -1.0
                conta = conta - 1
              end if
            end if
          end do
        end do
#endif

#ifdef NILE
        chym_dis(idamietta,jdamietta) = 0.5 * mval(nile_fresh_flux,imon,iday)
        chym_dis(irosetta,jrosetta) = 0.5 * mval(nile_fresh_flux,imon,iday)
        conta = conta + 2
#endif

#ifdef BLACKSEA
        chym_dis(idardanelli,jdardanelli) = mval(bs_fresh_flux,imon,iday)
        conta = conta + 1
#endif

        print*,"Number of discharge points accounted for : ",conta

      else

#ifdef CPL

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
          end do
        end do

#ifdef AZOV
        do j = 2 , chym_nlat-1
          do i = 2 , chym_nlon-1
            if ( lat1(j) > 45.00 .and. lon1(i) > 34.0 ) then
              if ( chym_dis(i,j) > 0.0 ) then
                chym_dis(i,j) = -1.0
                conta = conta - 1
              end if
            end if
          end do
        end do
        chym_dis(ikerch,jkerch) = 1.0
        conta = conta + 1
#endif

#ifdef NILE
        chym_dis(idamietta,jdamietta) = 1.0
        chym_dis(irosetta,jrosetta) = 1.0
        conta = conta + 2
#endif

#ifdef BLACKSEA
        chym_dis(idardanelli,jdardanelli) = 1.0
        conta = conta + 1
#endif

#endif
        print*,"Number of discharge points at beginning : ",conta
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

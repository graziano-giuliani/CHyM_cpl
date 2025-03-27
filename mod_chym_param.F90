module mod_chym_param
  implicit none

!-----------------------------------------------------------------------
!     CHYM model grid dimensions

!     nlc: number of longitudes
!     nbc: number of latitudes
!     fscal: resolution of the model
!-----------------------------------------------------------------------

  integer :: nlc
  integer :: nbc
  integer :: chym_nlon
  integer :: chym_nlat
  real, parameter :: umfang = 360.0
  real, allocatable :: areadom(:)
  integer, allocatable :: pointn(:)
  integer, allocatable :: npoint(:)
  real, allocatable :: portata(:)
  integer ir(9),jr(9)
  data ir /-1, 0, 1, 1, 1, 0,-1,-1,0/
  data jr / 1, 1, 1, 0,-1,-1,-1, 0,0/

!-----------------------------------------------------------------------
!     CHYM model I/O unit params
!-----------------------------------------------------------------------

  integer, parameter :: lu = 11
  integer, parameter :: lurun = 12
  integer, parameter :: lubas = 13
  integer, parameter :: lures = 14
  integer, parameter :: lupar = 15

!-----------------------------------------------------------------------
!     CHYM model netCDF data type
!-----------------------------------------------------------------------

  type CHYM_IO
    integer :: ncid
    integer, allocatable :: dimid(:)
    integer, allocatable :: varid(:)
    integer :: nrec
  end type CHYM_IO

  type(CHYM_IO) :: chymout, chymrst

!-----------------------------------------------------------------------
!     CHYM model config parameters
!-----------------------------------------------------------------------

  real :: thrriv
  logical :: restarted = .false.
  integer :: isread, iswrit
  integer :: jahr1, jahr2, jahr3, jahr4, nstep, pstep
  integer(8) :: sdate , edate
  character(len=32) :: calendar
  character(len=256) :: dnres, dnini, dnout, dnstt

!-----------------------------------------------------------------------
!     CHYM model parameters
!-----------------------------------------------------------------------

  ! area of the grid cells
  real, allocatable :: area(:)
  ! area of the chym grid cells
  real, allocatable :: chym_area(:,:)
  ! drainage (?)
  real, allocatable :: chym_drai(:,:)
  ! CHyM land-sea mask
  real, allocatable :: chym_lsm(:,:)
  ! Volume of river water
  real, allocatable :: h2o(:,:)
  ! Effective water flow velocity (m/s)
  real, allocatable :: alfa(:,:)
  ! Flow directions map
  integer, allocatable :: fmap(:,:)
  ! Acclivity map
  real, allocatable :: accl(:,:)
  ! Land use from CHyM
  integer, allocatable :: luse(:,:)
  ! BWET CHyM
  real, allocatable :: bwet(:,:)
  ! port CHyM
  real, allocatable :: port(:,:)
  ! wkm1 CHyM
  real, allocatable :: wkm1(:,:)
  ! Manning coeff
  real, allocatable :: manning(:)
  ! CHyM dx
  real, allocatable :: chym_dx(:,:)
  ! CHyM lat
  real, allocatable :: chym_lat(:,:)
  ! CHyM lon
  real, allocatable :: chym_lon(:,:)
  ! CHyM corner lat
  real, allocatable :: corner_lat(:,:,:)
  ! CHyM corner lon
  real, allocatable :: corner_lon(:,:,:)


  integer, parameter :: lntypes = 110
  integer, parameter :: ocean = 15
  integer, parameter :: nmemrf = 5
  integer, parameter :: nmemlf = 1
  integer, parameter :: mm = 4

  real, allocatable :: chym_runoff(:,:)
  real, allocatable :: chym_surf(:,:)
  real, allocatable :: oro(:,:)
  real, allocatable :: chym_dis(:,:)

#ifdef NILE
  real, parameter, dimension(12) :: nile_fresh_flux = &
      [ 336,  396,  407,  399,  472,  615,  &
        634,  557,  412,  373,  372,  352]
  integer :: idamietta , jdamietta
  integer :: irosetta , jrosetta
  real , parameter :: lat_damietta = 31.52
  real , parameter :: lon_damietta = 31.83
  real , parameter :: lat_rosetta = 31.46
  real , parameter :: lon_rosetta = 30.36
#endif
#ifdef BLACKSEA
  ! Kourafalou, V. H. and Barbopoulos, K.: High resolution
  ! simulations on the North Aegean Sea seasonal circulation, Ann.
  ! Geophys., 21, 251–265,
  ! https://doi.org/10.5194/angeo-21-251-2003, 2003.
  !
  ! Monthly Outflow from Table 2 (MOT2):
  !      5700, 7500,10000,12500,14300,15000,   => Mean 10000
  !     14300,12500,10000, 7500, 5700, 5000
  ! Salinity (SAT2):
  !     26.27,26.88,27.50,24.05,23.53,23.01,   => Mean 24.64
  !     22.50,23.13,23.77,24.41,25.05,25.66
  !
  ! García-García, D., Vigo, M.I., Trottini, M. et al. Hydrological
  ! cycle of the Mediterranean-Black Sea system. Clim Dyn 59,
  ! 1919–1938 (2022). https://doi.org/10.1007/s00382-022-06188-2
  !
  ! Total outflow from Black Sea to Mediterranean from EPR balance
  !        275 ± 59 km3/year => 8714 ± 1871 m^3/s
  !
  ! Mediterranean average salinity : 38.40 (Encyclopedia Britannica)
  !
  !  Used formula for FRESHWATER FLUX:
  !
  !        round(8714 * MOT2/10000 * (38.40-SAT2)/38.40)
  !
  ! Should we add some variability? How?
  !
  real, parameter, dimension(12) :: bs_fresh_flux = &
      [1600, 2000, 2500, 4100, 4800, 5200, &
       5200, 4300, 3300, 2400, 1700, 1400]
  real, parameter :: lat_dardanelli = 40.0
  real, parameter :: lon_dardanelli = 26.18
  integer :: idardanelli , jdardanelli
#endif

#ifdef AZOV
  integer :: ikerch , jkerch
  ! Control factor for discharge from azov sea
  real , parameter :: azovfac = 1.00
  real , parameter :: lon_kerch = 36.53
  real , parameter :: lat_kerch = 45.11
#endif

  contains

  logical function is_inbox(lat,lon,clat,clon)
    implicit none
    real, intent(in) :: lat, lon
    real, dimension(4), intent(in) :: clat, clon
    real :: m1, b1, m2, b2, m3, b3, m4, b4
    logical :: l1, l2, l3, l4
    is_inbox = .false.

    ! ASSUME NON PATOLOGICAL CASES, i.e. clon, clat proper quadrilater.
    if ( clat(1) /= clat(2) ) then
      m1 = (clat(1)-clat(2))/(clon(1)-clon(2))
      b1 = (clon(1)*clat(2)-clon(2)*clat(1))/(clon(1)-clon(2))
      l1 = (lat > m1*lon+b1)
    else
      l1 = lat > clat(1)
    end if
    if ( clat(3) /= clat(4) ) then
      m2 = (clat(3)-clat(4))/(clon(3)-clon(4))
      b2 = (clon(3)*clat(4)-clon(4)*clat(3))/(clon(3)-clon(4))
      l2 = lat < m2*lon+b2
    else
      l2 = lat < clat(3)
    end if
    if ( clon(1) /= clon(4) ) then
      m3 = (clon(1)-clon(4))/(clat(1)-clat(4))
      b3 = (clat(1)*clon(4)-clat(4)*clon(1))/(clat(1)-clat(4))
      l3 = lon > m3*lat+b3
    else
      l3 = lon > clon(1)
    end if
    if ( clon(2) /= clon(3) ) then
      m4 = (clon(2)-clon(3))/(clat(2)-clat(3))
      b4 = (clat(2)*clon(3)-clat(3)*clon(2))/(clat(2)-clat(3))
      l4 = lon < m4*lat+b4
    else
      l4 = lon < clon(3)
    end if
    is_inbox = l1 .and. l2 .and. l3 .and. l4
  end function is_inbox

  real function mval(series,imon,iday)
    implicit none
    integer , intent(in) :: imon, iday
    real, dimension(12) :: series
    integer :: imonp, imonn
    real :: w1, w2, w3
    real :: fm, hm
    real, dimension(12), parameter :: mpd = &
     [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
    imonp = imon - 1
    imonn = imon + 1
    if ( imonp < 1 ) imonp = 12
    if ( imonn > 12 ) imonn = 1
    fm = mpd(imon)
    hm = fm/2.0
    if ( iday < hm ) then
      w1 = 0.5 - (iday-1+0.5)/fm
      w2 = 1.0 - w1
      w3 = 0.0
    else
      w1 = 0.0
      w2 = 1.0 - (iday-hm+0.5)/fm
      w3 = 1.0 - w2
    end if
    mval = ( w1 * series(imonp) + &
             w2 * series(imon)  + &
             w3 * series(imonn) )
  end function mval

end module mod_chym_param

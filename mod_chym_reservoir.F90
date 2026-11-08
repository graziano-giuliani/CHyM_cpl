!
!    Simple reservoir module for CHyM (coupled version, serial).
!
!    Release rule of Hanasaki, Kanae and Oki (2006, J. Hydrol. 327, 22-41)
!    for non-irrigation reservoirs:
!
!      c  = C / (annual mean inflow volume)
!      k  = S(start of operational year) / (alpha C),   alpha = 0.85
!      r' = k * imean                                   (target release)
!      r  = r'                               if c >= 0.5
!      r  = (c/0.5)^2 r' + (1-(c/0.5)^2) i   if c <  0.5  (i: current inflow)
!
!    plus spill of any storage above capacity. The reservoir storage is the
!    storage (h2o) of the CHyM cell holding the dam. Losses: net open-water
!    evaporation from the surface (scaled with the filling S/C and with a
!    fixed monthly weight) and a constant consumptive withdrawal.
!
!    The table of reservoirs is a text file, one reservoir per line:
!      name lon lat upa_km2 cap_km3 surf_km2 evap_mm_yr use_m3s year0 opmon
!      qmean_m3s
!    upa_km2   : upstream area used to pick the dam cell (0: largest river
!                within 25 km)
!    year0     : first year of operation (the cell is a plain river before)
!    opmon     : first month of the operational year (k is reset then)
!    qmean_m3s : mean inflow (0: the mean discharge 'dis' of the static file;
!                in both cases it is then updated with a 1-year running mean)
!    Lines starting with '#' are comments.
!
module mod_chym_reservoir
  use iso_fortran_env, only : output_unit
  use mod_chym_param
  implicit none
  private

  type reservoir
    character(len=40) :: name
    real(8) :: lon, lat, upa, cap, surf, evap, use, qmean0
    integer :: year0, opmon
    integer :: ic = -1, jc = -1      ! dam cell (lon index, lat index)
    logical :: active = .false.
    real(8) :: imean = 0.0d0         ! running mean inflow (m3/s)
    real(8) :: iprev = 0.0d0         ! inflow over the last step (m3/s)
    real(8) :: kfac = 1.0d0          ! Hanasaki release coefficient
    real(8) :: release = 0.0d0       ! release for the current step (m3/s)
    ! accumulated over one input step, on the rank owning the dam cell
    real(8) :: vin = 0.0d0, vout = 0.0d0, vloss = 0.0d0, store = 0.0d0
    real :: dmlast = 0.0
  end type reservoir

  type(reservoir), allocatable, public :: res(:)
  integer, allocatable, public :: resid(:,:)
  integer, public :: nres = 0
  real(8), public :: qres = 0.0d0

  real(8), parameter :: alpha = 0.85d0
  real(8), parameter :: yrsec = 365.25d0*86400.0d0
  ! Monthly share of the annual open-water evaporation (mid-latitudes)
  real(8), parameter, dimension(12) :: evw = &
    [ 0.02d0, 0.02d0, 0.04d0, 0.07d0, 0.11d0, 0.14d0, &
      0.16d0, 0.15d0, 0.12d0, 0.08d0, 0.05d0, 0.04d0 ]

  public :: res_init, res_fill, res_active, res_outflow, res_update
  public :: res_daily, res_report, res_start

  contains

  ! Read the table and locate the dam cells. Called on every rank.
  subroutine res_init(fname, year)
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in) :: year
    character(len=512) :: line
    integer :: lun, ios, n, i, j, k
    real(8) :: d, best, score, up
    type(reservoir) :: r

    if (.not. allocated(resid)) allocate(resid(chym_nlon,chym_nlat))
    resid = 0
    nres = 0
    if ( len_trim(fname) == 0 ) return
    open(newunit=lun, file=trim(fname), status='old', action='read', &
         iostat=ios)
    if ( ios /= 0 ) then
      write(output_unit,*) 'Cannot open reservoir table ', trim(fname)
      error stop 'CHYM - reservoir table error'
    end if
    n = 0
    do
      read(lun,'(a)',iostat=ios) line
      if ( ios /= 0 ) exit
      if ( len_trim(line) == 0 .or. adjustl(line(1:1)) == '#' ) cycle
      n = n + 1
    end do
    allocate(res(n))
    rewind(lun)
    k = 0
    do
      read(lun,'(a)',iostat=ios) line
      if ( ios /= 0 ) exit
      if ( len_trim(line) == 0 .or. adjustl(line(1:1)) == '#' ) cycle
      k = k + 1
      read(line,*) r%name, r%lon, r%lat, r%upa, r%cap, r%surf, r%evap, &
                   r%use, r%year0, r%opmon, r%qmean0
      r%cap = r%cap * 1.0d9          ! km3 -> m3
      r%surf = r%surf * 1.0d6        ! km2 -> m2
      r%evap = r%evap * 1.0d-3       ! mm/yr -> m/yr
      res(k) = r
    end do
    close(lun)
    nres = n

    do k = 1, nres
      best = huge(1.0d0)
      do j = 2, chym_nlat-1
        do i = 2, chym_nlon-1
          if ( luse(i,j) == ocean .or. fmap(i,j) < 1 .or. &
               fmap(i,j) > 8 ) cycle
          d = gcdist(res(k)%lat,res(k)%lon, &
                     real(chym_lat(i,j),8),real(chym_lon(i,j),8))
          if ( d > 25.0d0 ) cycle
          up = real(chym_drai(i,j),8) + real(chym_area(i,j),8)
          if ( res(k)%upa > 0.0d0 ) then
            score = abs(log(up/res(k)%upa))
          else
            score = -up
          end if
          if ( score < best ) then
            best = score
            res(k)%ic = i
            res(k)%jc = j
          end if
        end do
      end do
      if ( res(k)%ic < 0 ) then
        write(output_unit,*) 'Reservoir ', trim(res(k)%name), &
          ': no river cell within 25 km, ignored'
        cycle
      end if
      if ( resid(res(k)%ic,res(k)%jc) /= 0 ) then
        write(output_unit,*) 'Reservoir ', trim(res(k)%name), &
          ': cell already used by another reservoir, ignored'
        res(k)%ic = -1
        cycle
      end if
      resid(res(k)%ic,res(k)%jc) = k
      res(k)%imean = res(k)%qmean0
      if ( res(k)%imean <= 0.0d0 ) then
        if ( .not. has_omdis ) then
          write(output_unit,*) 'Reservoir ', trim(res(k)%name), &
            ': no qmean in the table and no dis in the static file'
          error stop 'CHYM - reservoir table error'
        end if
        res(k)%imean = max(real(omdis(res(k)%ic,res(k)%jc),8),1.0d0)
      end if
      res(k)%iprev = res(k)%imean
      res(k)%active = ( year >= res(k)%year0 )
      if ( .true. ) then
        write(output_unit,'(1x,a,a20,a,2f8.3,a,f10.0,a,f9.0,a,f7.2,a,f8.1,a,l2)') &
          'CHYM - Reservoir ', trim(res(k)%name), ' cell ', &
          chym_lon(res(k)%ic,res(k)%jc), chym_lat(res(k)%ic,res(k)%jc), &
          ' upa ', real(chym_drai(res(k)%ic,res(k)%jc)) + &
                   real(chym_area(res(k)%ic,res(k)%jc)), &
          ' (table ', res(k)%upa, ') cap km3 ', res(k)%cap*1.0d-9, &
          ' qmean ', res(k)%imean, ' active', res(k)%active
      end if
    end do
  end subroutine res_init

  ! Make sure active reservoirs hold at least frac*C (cold start, or a
  ! restart written before the reservoirs existed). Works on any array
  ! indexed with global (lon, lat) indices; returns the volume added.
  real(8) function res_fill(h, frac) result(added)
    implicit none
    real, intent(inout) :: h(:,:)
    real, intent(in) :: frac
    integer :: k
    added = 0.0d0
    do k = 1, nres
      if ( res(k)%ic < 0 .or. .not. res(k)%active ) cycle
      if ( h(res(k)%ic,res(k)%jc) < 0.01*res(k)%cap ) then
        added = added + frac*res(k)%cap - h(res(k)%ic,res(k)%jc)
        h(res(k)%ic,res(k)%jc) = real(frac*res(k)%cap)
      end if
    end do
  end function res_fill

  ! After the storages are final (cold start fill or restart): set k and
  ! the first release from the actual storage, identically on every rank.
  subroutine res_start(h)
    implicit none
    real, intent(in) :: h(:,:)
    integer :: k
    real(8) :: c, w
    do k = 1, nres
      if ( res(k)%ic < 0 ) cycle
      res(k)%store = real(h(res(k)%ic,res(k)%jc),8)
      res(k)%kfac = res(k)%store / (alpha*res(k)%cap)
      c = res(k)%cap / (res(k)%imean*yrsec)
      w = min((c/0.5d0)**2, 1.0d0)
      res(k)%release = w*res(k)%kfac*res(k)%imean + (1.0d0-w)*res(k)%iprev
    end do
  end subroutine res_start

  logical function res_active(k)
    implicit none
    integer, intent(in) :: k
    res_active = .false.
    if ( k > 0 ) res_active = res(k)%active
  end function res_active

  ! Outflow volume over deltat from a dam cell with storage h (m3)
  real function res_outflow(k, h, deltat) result(dm)
    implicit none
    integer, intent(in) :: k
    real, intent(in) :: h, deltat
    dm = real(min(res(k)%release*deltat, real(h,8)))
    ! spill whatever exceeds the capacity after the release
    if ( h - dm > res(k)%cap ) dm = real(h - res(k)%cap)
    dm = max(dm, 0.0)
  end function res_outflow

  ! After the storage update of the dam cell: losses and accounting.
  ! inflow = wkm1 + own outflow + local runoff; owned: this rank owns it.
  subroutine res_update(k, h, inflow, dm, imon, deltat, owned)
    implicit none
    integer, intent(in) :: k, imon
    real, intent(inout) :: h
    real, intent(in) :: inflow, dm, deltat
    logical, intent(in) :: owned
    real(8) :: fill, loss
    fill = min(max(real(h,8)/res(k)%cap,0.0d0),1.0d0)
    loss = (res(k)%surf*fill*res(k)%evap*evw(imon)*12.0d0/yrsec + &
            res(k)%use) * deltat
    loss = min(loss, real(h,8))
    h = real(h - loss)
    if ( owned ) then
      res(k)%vin = res(k)%vin + inflow
      res(k)%vout = res(k)%vout + dm
      res(k)%vloss = res(k)%vloss + loss
      qres = qres + loss
    end if
  end subroutine res_update

  ! Once per input step, after the routing sub-steps: collect inflow and
  ! storage from the owning ranks, update the running mean, reset k at the
  ! start of the operational year and set the release for the next step.
  subroutine res_daily(h, dtin, year, month, day)
    implicit none
    integer, intent(in) :: year, month, day
    real, intent(in) :: h(:,:)
    real(8), intent(in) :: dtin
    integer :: k
    real(8) :: c, w, rp

    if ( nres == 0 ) return
    do k = 1, nres
      if ( res(k)%ic < 0 ) cycle
      res(k)%store = real(h(res(k)%ic,res(k)%jc),8)
      if ( .not. res(k)%active ) then
        if ( year >= res(k)%year0 ) then
          res(k)%active = .true.
          res(k)%kfac = res(k)%store / (alpha*res(k)%cap)
        end if
      end if
      if ( res(k)%active ) then
        res(k)%iprev = res(k)%vin / dtin
        res(k)%imean = res(k)%imean + &
                       (res(k)%iprev - res(k)%imean) * dtin / yrsec
        if ( month == res(k)%opmon .and. day == 1 ) then
          res(k)%kfac = res(k)%store / (alpha*res(k)%cap)
        end if
        c = res(k)%cap / (res(k)%imean*yrsec)
        rp = res(k)%kfac * res(k)%imean
        if ( c >= 0.5d0 ) then
          res(k)%release = rp
        else
          w = (c/0.5d0)**2
          res(k)%release = w*rp + (1.0d0-w)*res(k)%iprev
        end if
        res(k)%release = max(res(k)%release, 0.0d0)
      end if
    end do
  end subroutine res_daily

  ! Per-step report (root): one line per active reservoir to a text file,
  ! one summary line to the log. Resets the accumulators.
  subroutine res_report(lun, now, dtin)
    implicit none
    integer, intent(in) :: lun
    character(len=*), intent(in) :: now
    real(8), intent(in) :: dtin
    integer :: k
    real(8) :: s, qi, qo
    if ( nres == 0 ) return
    s = 0.0d0 ; qi = 0.0d0 ; qo = 0.0d0
    do k = 1, nres
      if ( res(k)%ic < 0 .or. .not. res(k)%active ) cycle
      if ( lun /= 0 ) then
        write(lun,'(a,1x,a20,f9.4,4f11.1,f8.3)') trim(now), &
          trim(res(k)%name), res(k)%store*1.0d-9, res(k)%vin/dtin, &
          res(k)%vout/dtin, res(k)%vloss/dtin, res(k)%imean, res(k)%kfac
      end if
      s = s + res(k)%store
      qi = qi + res(k)%vin/dtin
      qo = qo + res(k)%vout/dtin
      res(k)%vin = 0.0d0
      res(k)%vout = 0.0d0
      res(k)%vloss = 0.0d0
    end do
    write(output_unit,'(1x,a,f10.3,a,f10.1,a,f10.1)') &
      'CHYM - Reservoirs: storage km3', s*1.0d-9, '  inflow m3/s', qi, &
      '  outflow m3/s', qo
  end subroutine res_report

  real(8) function gcdist(la1, lo1, la2, lo2)
    implicit none
    real(8), intent(in) :: la1, lo1, la2, lo2
    real(8), parameter :: d2r = 0.017453292519943295d0
    real(8) :: a
    a = sin((la2-la1)*d2r/2)**2 + &
        cos(la1*d2r)*cos(la2*d2r)*sin((lo2-lo1)*d2r/2)**2
    gcdist = 2.0d0*6371.0d0*asin(sqrt(min(a,1.0d0)))
  end function gcdist

end module mod_chym_reservoir

module mod_chym_model

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

  use, intrinsic :: iso_fortran_env
  use mod_chym_param
  use mod_chym_reservoir

  implicit none
  private

!-----------------------------------------------------------------------
!     Public subroutines
!-----------------------------------------------------------------------

  public :: chymmodel , coupler , mouth_snapshot

  contains

    subroutine chymmodel(chym_runoff,chym_dis,imon,iday)
      implicit none

!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------

      real, intent(inout) :: chym_runoff(:,:) ! m/s
      real, intent(inout) :: chym_dis(:,:)
      integer, intent(in) :: imon, iday

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

      integer :: i, j, ii, jj, imin, ilnd, idir
      real :: dm, area, deltat, h2oload, irtmp, hmin
      ! Daily water budget (m^3), double precision accumulators
      real(8) :: qin, qsea, qout, qsink, qirr, qmin, store0, store1
      integer, save :: lastmon = -1, lastday = -1, lunres = 0
      character(len=10) :: resdate
      real(8) :: daysec

      ! chym_dis accumulates, over the day, the volume (m^3) leaving each
      ! land cell into a sea cell of (chym_mask == 1).
      ! It is converted to a daily mean discharge (m^3/s) at the end.
      chym_dis(:,:) = 0.0
      where (chym_mask /= 2)
        chym_runoff = 0.0
      end where

      deltat = 86400.0/real(model_nsteps)
      daysec = real(deltat,8)*real(model_nsteps,8)
      ! Storage floor (m^3): zero, as in CHyM-roff. Must be set: it was
      ! left uninitialised when minh2o was removed.
      hmin = 0.0
      qin = 0.0d0
      qsea = 0.0d0
      qout = 0.0d0
      qsink = 0.0d0
      qirr = 0.0d0
      qmin = 0.0d0
      qres = 0.0d0
      store0 = routed_storage( )
      ! Model year, for reservoirs coming into operation during the run
      if ( imon == 1 .and. iday == 1 .and. &
           .not. (lastmon == 1 .and. lastday == 1) .and. lastmon > 0 ) &
        chym_year = chym_year + 1
      lastmon = imon
      lastday = iday

      do imin = 1, model_nsteps
        wkm1(:,:) = 0.0
        do j = 2, chym_nlat-1
          do i = 2, chym_nlon-1
            idir = fmap(i,j)
            ilnd = luse(i,j)
            if ( ilnd .ne. ocean .and. idir .ge. 1 .and. idir .le. 8 ) then
              ii = i+ir(idir)
              jj = j+jr(idir)
              if ( res_active(resid(i,j)) ) then
                ! Dam: release rule instead of the linear reservoir
                dm = res_outflow(resid(i,j), h2o(i,j), deltat)
                res(resid(i,j))%dmlast = dm
              else
                dm = min(port(i,j)*deltat*efficiency, h2o(i,j)*efficiency)
              end if
              wkm1(i,j) = wkm1(i,j) - dm
              if ( luse(ii,jj) == ocean ) then
                if ( chym_mask(ii,jj) == 1 ) then
                  ! Into the coupled sea: delivered to the outlet this
                  ! mouth is aggregated into (see aggregate_mouths)
                  chym_dis(colx(ii,jj),coly(ii,jj)) = &
                    chym_dis(colx(ii,jj),coly(ii,jj)) + dm
                  qsea = qsea + dm
                else
                  ! Out of the catchments
                  qout = qout + dm
                end if
              else if ( fmap(ii,jj) < 1 .or. fmap(ii,jj) > 8 .or. &
                        ii == 1 .or. ii == chym_nlon .or.         &
                        jj == 1 .or. jj == chym_nlat ) then
                ! Endorheic sink (fdm = 9) or domain frame: water leaves
                qsink = qsink + dm
              else
                wkm1(ii,jj) = wkm1(ii,jj) + dm
              end if
            end if
          end do
        end do
        do j = 2, chym_nlat-1
          do i = 2, chym_nlon-1
            idir = fmap(i,j)
            ilnd = luse(i,j)
            if ( ilnd .ne. ocean .and. idir .ge. 1 .and. idir .le. 8 ) then
               ! m^3 of water recharge in the grid cell
               ! Area in the input file is in km^2, we put it in m^2
               ! m^2 * m/s * s = m^3
               area = chym_area(i,j)*1.0e6
               h2oload = area*chym_runoff(i,j)*deltat
               if (h2oload > 200000.0) then
                 write(error_unit,fmt='(A,F8.2)') &
                         "CHYM - *WARNING VERY BIG H2O RELOAD* : ", h2oload
               endif
               h2o(i,j) = h2o(i,j) + wkm1(i,j) + h2oload
               qin = qin + h2oload
               if ( h2o(i,j) < hmin ) then
                 ! Water created to keep the minimum storage
                 qmin = qmin + (hmin - h2o(i,j))
                 h2o(i,j) = hmin
               end if
               if ( res_active(resid(i,j)) ) then
                 ! Dam cell: evaporation and withdrawal from the
                 ! reservoir; port reports the release.
                 call res_update(resid(i,j), h2o(i,j),                  &
                        wkm1(i,j)+res(resid(i,j))%dmlast+h2oload,       &
                        res(resid(i,j))%dmlast, imon, deltat, .true.)
                 bwet(i,j) = h2o(i,j) / chym_dx(i,j)
                 port(i,j) = res(resid(i,j))%dmlast / deltat
                 cycle
               end if
               if ( farm(i,j) ) then
                 irtmp = irmonfac(imon)/deltat
                 ! Some water is not transmitted for irrigation
                 bwet(i,j) = (1.0-irtmp) * h2o(i,j) / chym_dx(i,j)
                 ! Evaporative loss of the water for irrigation
                 qirr = qirr + irloss*irtmp*h2o(i,j)
                 h2o(i,j) = (1.0-irloss*irtmp) * h2o(i,j)
               else
                 bwet(i,j) = h2o(i,j) / chym_dx(i,j)
               end if
               port(i,j) = alfa(i,j) * bwet(i,j)
            end if
          end do
        end do
      end do

      ! Daily mean discharge into each sea cell (m^3/s)
      chym_dis(:,:) = real(real(chym_dis(:,:),8)/daysec)
      store1 = routed_storage( )
      if ( nres > 0 ) then
        call res_daily(h2o, daysec, chym_year, imon, iday)
        if ( lunres == 0 ) then
          open(newunit=lunres, file='chym_reservoirs.txt', &
               status='replace', action='write')
          write(lunres,'(a)') '# year-mm-dd name storage_km3 inflow_m3s'// &
            ' outflow_m3s losses_m3s mean_inflow_m3s k'
        end if
        write(resdate,'(i4.4,a,i2.2,a,i2.2)') chym_year,'-',imon,'-',iday
        call res_report(lunres, resdate, daysec)
      end if

      call coupler(imon,iday)

      ! Daily water budget of the routing, as mean fluxes in m^3/s.
      ! residual = in + created - (sea + out + sinks + irrigation
      !            + reservoir losses) - dS/dt
      write(output_unit,fmt='(1x,A)') 'CHYM - Daily budget (m3/s): '// &
        '   runoff_in    to_sea   out_catch     sinks  irrigat.  '// &
        'created     dS/dt  residual reservoir'
      write(output_unit,fmt='(1x,A,9F10.1)') 'CHYM - Budget             : ', &
        qin/daysec, qsea/daysec, qout/daysec, qsink/daysec, qirr/daysec, &
        qmin/daysec, (store1-store0)/daysec, &
        (qin+qmin-qsea-qout-qsink-qirr-qres-(store1-store0))/daysec, &
        qres/daysec

      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Discharge max value: ", &
                             maxval(port)
      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Bwet max value     : ", &
                             maxval(bwet)
      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Alfa max value     : ", &
                             maxval(alfa)
      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - H2o max value      : ", &
                             maxval(h2o)
      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Wkm1 max value     : ", &
                             maxval(wkm1)
      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Dx max value       : ", &
                             maxval(chym_dx)
      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Area max value     : ", &
                             maxval(chym_area)
      write(output_unit,fmt='(1x,A,F16.2,A)')"CHYM - Runoff max value   : ", &
                             maxval(chym_runoff)*86400*1000.0, ' mm/day'
    end subroutine chymmodel

    real(8) function routed_storage( ) result(s)
      implicit none
      integer :: i, j, idir
      s = 0.0d0
      do j = 2, chym_nlat-1
        do i = 2, chym_nlon-1
          idir = fmap(i,j)
          if ( luse(i,j) /= ocean .and. idir >= 1 .and. idir <= 8 ) then
            s = s + real(h2o(i,j),8)
          end if
        end do
      end do
    end function routed_storage

    ! Instantaneous discharge from the current state, used only to provide
    ! a first field at a warm restart before any routing step is done.
    subroutine mouth_snapshot( )
      implicit none
      integer :: i, j, ii, jj, idir
      chym_dis(:,:) = 0.0
      do j = 2 , chym_nlat-1
        do i = 2 , chym_nlon-1
          idir = fmap(i,j)
          if ( luse(i,j) /= ocean .and. idir >= 1 .and. idir <= 8 ) then
            ii = i + ir(idir)
            jj = j + jr(idir)
            if ( luse(ii,jj) == ocean .and. chym_mask(ii,jj) == 1 ) then
              chym_dis(colx(ii,jj),coly(ii,jj)) = &
                chym_dis(colx(ii,jj),coly(ii,jj)) + port(i,j)
            end if
          end if
        end do
      end do
    end subroutine mouth_snapshot

    ! Apply the prescribed discharges on top of the routed ones.
    ! The routed discharge itself is computed in chymmodel.
    subroutine coupler(imon,iday)
      implicit none
      integer, intent(in) :: imon, iday
      integer :: i, j
      real :: tmp

#ifdef AZOV
      tmp = 0.0
      do j = 2 , chym_nlat-1
        do i = 2 , chym_nlon-1
          if ( chym_lat(i,j) > 45.00 .and. chym_lon(i,j) > 34.0 ) then
            if ( chym_dis(i,j) > 0.0 ) then
              tmp = tmp + chym_dis(i,j)
              chym_dis(i,j) = 0.0
            end if
          end if
        end do
      end do
      chym_dis(ikerch,jkerch) = tmp
      write(output_unit,fmt='(1x,A,F16.2)') &
              "CHYM - Azov sea added up discharge value:   ",tmp
#endif

#ifdef NILE
      tmp = 0.5 * mval(nile_fresh_flux,imon,iday)
      write(output_unit,fmt='(1x,A,F16.2)') &
        "CHYM - Routed water replaced at Nile outlets:", &
        chym_dis(idamietta,jdamietta) + chym_dis(irosetta,jrosetta)
      chym_dis(idamietta,jdamietta) = tmp
      chym_dis(irosetta,jrosetta) = tmp
      write(output_unit,fmt='(1x,A,F16.2)') &
              "CHYM - Prescribed Nile discharge :   ",tmp*2
#endif

#ifdef BLACKSEA
      tmp = mval(bs_fresh_flux,imon,iday)
      chym_dis(idardanelli,jdardanelli) = tmp
      write(output_unit,fmt='(1x,A,F16.2)') &
              "CHYM - Prescribed Dardanelli BS output :   ",tmp
#endif

    end subroutine coupler

end module mod_chym_model

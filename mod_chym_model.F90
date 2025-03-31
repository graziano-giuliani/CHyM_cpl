module mod_chym_model

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

  use, intrinsic :: iso_fortran_env
  use mod_chym_param

  implicit none
  private

!-----------------------------------------------------------------------
!     Public subroutines
!-----------------------------------------------------------------------

  public :: chymmodel , coupler

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
      real :: dm, area, deltat, rainload

      chym_dis(:,:) = 0.0
      where (chym_mask > 0)
        chym_runoff = 0.0
      end where

      deltat = 86400.0/real(model_nsteps)
      do imin = 1, model_nsteps
        wkm1(:,:) = 0.0
        do j = 2, chym_nlat-1
          do i = 2, chym_nlon-1
            idir = fmap(i,j)
            ilnd = luse(i,j)
            if ( ilnd .ne. ocean .and. idir .ge. 1 .and. idir .le. 8 ) then
              ii = i+ir(idir)
              jj = j+jr(idir)
              dm = min(port(i,j)*deltat,h2o(i,j))
              wkm1(i,j) = wkm1(i,j) - dm
              wkm1(ii,jj) = wkm1(ii,jj) + dm
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
               area = chym_area(i,j)*1.0e+06
               rainload = area*chym_runoff(i,j)*deltat
               if (rainload > 200000.0) then
                 write(error_unit,fmt='(A,F8.2)') &
                         "CHYM - *WARNING VERY BIG RAINLOAD*:", rainload
               endif
               h2o(i,j) = h2o(i,j) + wkm1(i,j) + rainload
               h2o(i,j) = max(h2o(i,j),0.0)
               bwet(i,j) = h2o(i,j) / chym_dx(i,j)
               port(i,j) = alfa(i,j) * bwet(i,j)
            end if
          end do
        end do
      end do

      call coupler(imon,iday)

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

    subroutine coupler(imon,iday)
      implicit none
      integer, intent(in) :: imon, iday
      integer :: i, j, ii, jj, ilnd, idir
      real :: tmp

      do j = 2 , chym_nlat-1
        do i = 2 , chym_nlon-1
          idir = fmap(i,j)
          if ( idir >= 1 .and. idir <= 8 ) then
            ii = i + ir(idir)
            jj = j + jr(idir)
            ilnd = luse(ii,jj)
            if ( chym_drai(i,j) > thrriv .and. ilnd == ocean ) then
              chym_dis(ii,jj) = port(i,j)
            end if
          end if
        end do
      end do

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

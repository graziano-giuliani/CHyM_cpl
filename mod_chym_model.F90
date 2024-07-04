module mod_chym_model
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
    use mod_chym_param
!
    implicit none
    private
!
!-----------------------------------------------------------------------
!     Public subroutines
!-----------------------------------------------------------------------
!
    public :: chymmodel
!
    contains
!
    subroutine chymmodel(chym_runoff,chym_surf,chym_dis,imon,iday)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      real, intent(in) :: chym_runoff(:,:) ! m/s
      real, intent(in) :: chym_surf(:,:)   ! m/s
      real, intent(inout) :: chym_dis(:,:)
      integer, intent(in) :: imon, iday
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: i, j, ii, jj, imin, ilnd, idir
      real :: dm, step, area, deltat, rainload, tmp

      chym_dis(:,:) = -1.0

      step = 600.0         ! Number of step per days
      deltat = 86400.0/step
      oro = chym_runoff + chym_surf
      do imin = 1, step
        wkm1(:,:) = 0.0
        do j = 2, chym_nlat-1
          do i = 2, chym_nlon-1
            idir = int(fmap(i,j))
            ilnd = int(luse(i,j))
            if ( ilnd .ne. mare .and. idir .ge. 1 .and. idir .le. 8 ) then
              ii = i+ir(idir)
              jj = j+jr(idir)
              dm = port(i,j)*deltat
!             write(6,'(12x,2i4,2f9.4)') i,j,dm,port(i,j)
              wkm1(i,j) = wkm1(i,j) - dm
!             write(6,'(12x,2i4,2f9.4)') i,j,wkm1(i,j),port(i,j)
              wkm1(ii,jj) = wkm1(ii,jj) + dm
            end if
          end do
        end do
        do j = 2, chym_nlat-1
          do i = 2, chym_nlon-1
            idir = int(fmap(i,j))
            ilnd = int(luse(i,j))
            if ( ilnd .ne. mare .and. idir .ge. 1 .and. idir .le. 8 ) then
               ! m^3 of water recharge in the grid cell
               ! Area in the input file is in km^2, we put it in m^2
               ! m^2 * m/s * s = m^3
               area = chym_area(i,j)*1.0e+06
               rainload = area*(chym_runoff(i,j)+chym_surf(i,j))*deltat
               if (rainload > 200000.0) then
                 write(*,fmt='(A,F8.2)')"*WARNING VERY BIG RAINLOAD*:",&
                                       rainload
               endif
               h2o(i,j) = h2o(i,j) + wkm1(i,j) + rainload
               h2o(i,j) = max(h2o(i,j),0.0)
               bwet(i,j) = h2o(i,j) / chym_dx(i,j)
               port(i,j) = alfa(i,j) * bwet(i,j)
            end if
          end do
        end do
      end do

      do j = 2 , chym_nlat-1
        do i = 2 , chym_nlon-1
          idir = int(fmap(i,j))
          if ( idir >= 1 .and. idir <= 8 ) then
            ilnd = int(luse(i+ir(idir),j+jr(idir)))
            !5400 km^2 ~ 6 grid points, with resolution of 33km, drained
            if ( chym_drai(i,j) > thrriv .and. ilnd == mare ) then
              chym_dis(i,j) = port(i,j)
            end if
          end if
        end do
      end do

#ifdef AZOV
      tmp = 0.0
      do j = 2 , chym_nlat-1
        do i = 2 , chym_nlon-1
          if ( lat1(j) > 45.00 .and. lon1(i) > 34.0 ) then
            if ( chym_dis(i,j) > 0.0 ) then
              tmp = tmp + chym_dis(i,j)
              chym_dis(i,j) = -1.0
            end if
          end if
        end do
      end do
      chym_dis(ikerch,jkerch) = tmp
      write(*,fmt='(A,F16.2)') "Azov sea added up discharge value:   ",tmp
#endif

#ifdef NILE
      tmp = 0.5 * mval(nile_fresh_flux,imon,iday)
      chym_dis(idamietta,jdamietta) = tmp
      chym_dis(irosetta,jrosetta) = tmp
      write(*,fmt='(A,F16.2)') "Prescribed Nile discharge :   ",tmp*2
#endif

#ifdef BLACKSEA
      tmp = mval(bs_fresh_flux,imon,iday)
      chym_dis(idardanelli,jdardanelli) = tmp
      write(*,fmt='(A,F16.2)') "Prescribed Dardanelli BS output :   ",tmp
#endif

      write(*,fmt='(A,F16.2)')"Discharge max value:   ",maxval(port)
      write(*,fmt='(A,F16.2)')"Bwet max value:   ",maxval(bwet)
      write(*,fmt='(A,F16.2)')"Alfa max value:   ",maxval(alfa)
      write(*,fmt='(A,F16.2)')"h2o max value:   ",maxval(h2o)
      write(*,fmt='(A,F16.2)')"wkm1 max value:   ",maxval(wkm1)
      write(*,fmt='(A,F16.2)')"dx max value:   ",maxval(chym_dx)
      write(*,fmt='(A,F16.2)')"chym_area max value:    ",              &
                             maxval(chym_area)
      write(*,fmt='(A,F16.2)')"chym_runoff max value:    ",            &
                             maxval(chym_runoff)
      write(*,fmt='(A,F16.2)')"chym_surf max value:    ",              &
                             maxval(chym_surf)
    end subroutine chymmodel

end module mod_chym_model

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
      real, intent(in) :: chym_runoff(:,:)
      real, intent(in) :: chym_surf(:,:)
      real, intent(inout) :: chym_dis(:,:)
      integer, intent(in) :: imon, iday
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: i, j, imin, ilnd, deltat, rainload,idir,step
      real :: dm

      chym_dis = -1.0
      step = 600           !Number of step per days
      deltat = 86400/step
      do imin = 1, step
        wkm1 = 0.0
        do j = 2, chym_nlat-1
          do i = 2, chym_nlon-1
            idir = int(fmap(i,j))
            ilnd = int(luse(i,j))
            if ( ilnd .ne. mare .and. idir .ge. 1 .and. idir .le. 8 ) then
              dm = port(i,j)*deltat
!              write(6,'(12x,2i4,2f9.4)') i,j,dm,port(i,j)
              if (dm.gt.h2o(i,j)) dm=h2o(i,j)
              wkm1(i,j)=wkm1(i,j)-dm
!              write(6,'(12x,2i4,2f9.4)') i,j,wkm1(i,j),port(i,j)
              wkm1(i+ir(idir),j+jr(idir))=wkm1(i+ir(idir),j+jr(idir))+dm
            endif
          enddo
        enddo
        do ji = 2, chym_nlat-1
          do i = 2, chym_nlon-1
            idir = int(fmap(i,j))
            ilnd = int(luse(i,j))
            if ( ilnd .ne. mare .and. idir .ge. 1 .and. idir .le. 8 ) then
               ! m3 of water recharge in the  grid cell
               rainload = chym_area(i,j)*1.0e+06*(chym_runoff(i,j)+    &
                  chym_surf(i,j))*deltat
               if (rainload.gt.200000.0) then
                write(*,fmt='(A,F8.2)')"*WARNING VERY BIG RAINLOAD*:",&
                                       rainload
               endif
               h2o(i,j) = h2o(i,j) + wkm1(i,j) + rainload
               bwet(i,j) = h2o(i,j) / chym_dx(i,j)
               port(i,j) = alfa(i,j) * bwet(i,j)
            endif
          enddo
        enddo
      enddo

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
#ifdef NILE
          ! Booca 1 : Lat 31.52, Lon 31.84 (Damietta) 50%
          ! Bocca 2 : Lat 31.47, Lon 30.36 (Rosetta)  50%
          if ( sqrt(((lat1(j)-31.525)**2 + &
                     (lon1(i)-31.843)**2)) < 0.5*fscal .or. &
               sqrt(((lat1(j)-31.467)**2 + &
                     (lon1(i)-30.366)**2)) < 0.5*fscal ) then
            chym_dis(i,j) = 0.5 * mval(nile_fresh_flux,imon,iday)
          end if
#endif
#ifdef BLACKSEA
          if ( sqrt(((lat1(j)-40.00)**2 + &
                     (lon1(i)-26.16)**2)) < 0.5*fscal ) then
            chym_dis(i,j) = mval(bs_fresh_flux,imon,iday)
          end if
#endif
        end do
      end do

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
!!
      end subroutine chymmodel
!!
      end module mod_chym_model

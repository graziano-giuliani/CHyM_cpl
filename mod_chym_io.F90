!
!-----------------------------------------------------------------------
!     Module for CHyM model I/O subroutines
!-----------------------------------------------------------------------
!
      module mod_chym_io
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
      public :: read_config
      public :: read_init
      public :: chym_out_init
      public :: chym_out
      public :: chym_rst_init
      public :: chym_rst
      public :: chym_ini
!
      contains
!
      subroutine read_config(ifile)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      character(len=*) :: ifile
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      real :: fdum
      character(len=100) :: pname
!
!-----------------------------------------------------------------------
!     Read parameters
!-----------------------------------------------------------------------
!
      call read_rec(ifile, 'IOUT', fdum, pname)
      iout = int(fdum+0.001)
      call read_rec(ifile, 'ISREAD', fdum, pname)
      isread = int(fdum+0.001)
      call read_rec(ifile, 'ISWRIT', fdum, pname)
      iswrit = int(fdum+0.001)
      call read_rec(ifile, 'TIDATE', fdum, tdate)
      read(tdate(1:4),'(i)') jahr1
      read(tdate(5:6),'(i)') jahr2
      read(tdate(7:8),'(i)') jahr3
      print*,"year",jahr1,"month",jahr2,"day",jahr3
      call read_rec(ifile, 'NSTEP', fdum, pname)
      nstep = int(fdum+0.001)
      fdum = 0.0
      call read_rec(ifile, 'UFAKRU', fdum, pname)
      ufakru = fdum
      call read_rec(ifile, 'UTHR', fdum, pname)
      thrriv = fdum
!
      call read_rec(ifile, 'TDNRES', fdum, dnres)
      call read_rec(ifile, 'TDNINI', fdum, dnini)
      call read_rec(ifile, 'TDNOUT', fdum, dnout)
      call read_rec(ifile, 'TDNSTK', fdum, chym_statikin)
!
      end subroutine read_config
!
      subroutine read_init()

      use netcdf
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: i,j,idir,ilnd
      integer :: ncid, varid
      integer :: contat
!
!-----------------------------------------------------------------------
!     Allocate variables
!-----------------------------------------------------------------------
!
      if (.not. allocated(chym_area)) allocate(chym_area(nlc,nbc))
      if (.not. allocated(chym_drai)) allocate(chym_drai(nlc,nbc))
      if (.not. allocated(chym_lsm)) allocate(chym_lsm(nlc,nbc))
      if (.not. allocated(chym_dx)) allocate(chym_dx(nlc,nbc))
      if (.not. allocated(chym_lat)) allocate(chym_lat(nlc,nbc))
      if (.not. allocated(chym_lon)) allocate(chym_lon(nlc,nbc))
      if (.not. allocated(lat1)) allocate(lat1(nbc))
      if (.not. allocated(lon1)) allocate(lon1(nlc))
      if (.not. allocated(alfa)) allocate(alfa(nlc,nbc))
      if (.not. allocated(fmap)) allocate(fmap(nlc,nbc))
      if (.not. allocated(accl)) allocate(accl(nlc,nbc))
      if (.not. allocated(luse)) allocate(luse(nlc,nbc))
      if (.not. allocated(port)) allocate(port(nlc,nbc))
      if (.not. allocated(wkm1)) allocate(wkm1(nlc,nbc))
      if (.not. allocated(bwet)) allocate(bwet(nlc,nbc))
      if (.not. allocated(h2o)) allocate(h2o(nlc,nbc))
      if (.not. allocated(manning)) allocate(manning(lntypes))

      port=0; wkm1=0; bwet=0; h2o=0
      chym_lsm(:,:) = 1.0
!
!-----------------------------------------------------------------------
!     Set coefficients
!-----------------------------------------------------------------------
!
      manning = 0.043
!
!-----------------------------------------------------------------------
!     Open, read and close restart file
!-----------------------------------------------------------------------
!
      pstep = 0
      if (isread /= 0) then
        print*, "read chym restart data"
        call chym_ini()
      end if
!
!-----------------------------------------------------------------------
!     Open, read and close ocean-atmosphere land-sea mask file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_open(trim(chym_statikin),                     &
         nf90_nowrite, ncid))

      call nio_check(nf90_inq_varid(ncid, 'lon', varid))
      call nio_check(nf90_get_var(ncid, varid, lon1))

      call nio_check(nf90_inq_varid(ncid, 'lat', varid))
      call nio_check(nf90_get_var(ncid, varid, lat1))

      call nio_check(nf90_inq_varid(ncid, 'fdm', varid))
      call nio_check(nf90_get_var(ncid, varid, fmap))

      call nio_check(nf90_inq_varid(ncid, 'acc', varid))
      call nio_check(nf90_get_var(ncid, varid, accl))

      call nio_check(nf90_inq_varid(ncid, 'lus', varid))
      call nio_check(nf90_get_var(ncid, varid, luse))

      call nio_check(nf90_inq_varid(ncid, 'aer', varid))
      call nio_check(nf90_get_var(ncid, varid, chym_area))

      call nio_check(nf90_inq_varid(ncid, 'dra', varid))
      call nio_check(nf90_get_var(ncid, varid, chym_drai))

      call nio_check(nf90_close(ncid))

      do i=1,nlc
        chym_lat(i,:) = lat1
      enddo
      do i=1,nbc
        chym_lon(:,i) = lon1
      enddo

      contat = 0
      do j = 2, chym_nlat-1
        do i = 2, chym_nlon-1
          idir = int(fmap(i,j))
          if ( idir >= 1 .and. idir <= 8 ) then   !5400 km^2 ~ 6 grid
            ilnd = int(luse(i+ir(idir),j+jr(idir)))
            if ( chym_drai(i,j) > thrriv .and. ilnd == mare ) then
              chym_lsm(i,j) = 0.0
              contat = contat + 1
            end if
          end if
#ifdef NILE
          ! Booca 1 : Lat 31.52, Lon 31.84 (Damietta) 50%
          ! Bocca 2 : Lat 31.47, Lon 30.36 (Rosetta)  50%
          if ( sqrt(((lat1(j)-31.525)**2 + &
                     (lon1(i)-31.843)**2)) < 0.03 .or. &
               sqrt(((lat1(j)-31.467)**2 + &
                     (lon1(i)-30.366)**2)) < 0.03) then
            chym_lsm(i,j) = 0.0
            contat = contat + 1
          end if
#endif
#ifdef BLACKSEA
          if ( sqrt(((lat1(j)-40.00)**2 + &
                     (lon1(i)-26.16)**2)) < 0.5*fscal ) then
            chym_lsm(i,j) = 0.0
            contat = contat + 1
          end if
#endif
        end do
      end do
      call runoffspeed
!
      print*,"Number of river mouths found : ",contat
      end subroutine read_init
!
      subroutine read_rec(ifile, key, value, pname)
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      character(len=*), intent(in) :: ifile, key
      character(len=*), intent(inout) :: pname
      real, intent(inout) :: value
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: ios
!
!-----------------------------------------------------------------------
!     Open configuration file
!-----------------------------------------------------------------------
!
      open(lu, file=trim(ifile), access='sequential', form='formatted', &
           status='old', iostat=ios)
      if (ios /= 0) then
         write(*,*) '[error] -- file '//trim(ifile)//' not found!'
         stop
      endif
!
!-----------------------------------------------------------------------
!     Read and find the parameter
!-----------------------------------------------------------------------
!
      ios = 0
      do while (ios == 0)
        read(lu, fmt='(A80)', iostat=ios) pname
        if (index(pname, key) /= 0) then
          if (key(1:1) == 't' .or. key(1:1) == 'T') then
            read(lu, '(A80)') pname
            write(*,*) 'Parameter: '//trim(key)//' = ', trim(pname)
          else
            read(lu, *) value
            write(*,*) 'Parameter: '//trim(key)//' = ', value
          end if
          exit
        end if
      end do
!
!-----------------------------------------------------------------------
!     Close configuration file
!-----------------------------------------------------------------------
!
      close(lu, status='keep', iostat=ios)
      if (ios /= 0) then
        write(*,*) '[error] -- '//trim(ifile)//'is not closed!'
      endif
      end subroutine read_rec
!
      subroutine chym_out_init()
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      character(len=100) :: str
!
!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
      if (.not. allocated(chymout%dimid)) allocate(chymout%dimid(3))
      if (.not. allocated(chymout%varid)) allocate(chymout%varid(4))
!
      call nio_check(nf90_create(trim(dnout), nf90_clobber,             &
                                  chymout%ncid))
!
!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_dim(chymout%ncid, 'lon',                  &
                                  nlc, chymout%dimid(1)))
      call nio_check(nf90_def_dim(chymout%ncid, 'lat',                  &
                                  nbc, chymout%dimid(2)))
      call nio_check(nf90_def_dim(chymout%ncid, 'time',                 &
                                  nf90_unlimited, chymout%dimid(3)))
!
!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymout%ncid, 'lon', nf90_real,       &
                     chymout%dimid(1), chymout%varid(1)))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymout%ncid, 'lat', nf90_real,       &
                     chymout%dimid(2), chymout%varid(2)))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                     'long_name', 'Latitude'))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                     'units', 'degrees_north'))
!
      call nio_check(nf90_def_var(chymout%ncid, 'time', nf90_int,       &
                     chymout%dimid(3), chymout%varid(3)))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                     'long_name', 'Time'))
      write(str,fmt='("days since ",I4,"-",I2.2,"-",I2.2," 00:00:00")') &
           jahr1, jahr2, jahr3
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                     'units', str))
!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymout%ncid, 'dis', nf90_real,       &
                     chymout%dimid, chymout%varid(4)))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'long_name', 'River Discharge'))
      call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'missing_value', 1.0e20))
!
!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------
!
      call nio_check(nf90_enddef(chymout%ncid))
!
!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(1), lon1))
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(2), lat1))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymout%ncid))
!
      end subroutine chym_out_init
!
      subroutine chym_rst_init()
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      character(len=100) :: str
!
!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------
!
      if (.not. allocated(chymrst%dimid)) allocate(chymrst%dimid(3))
      if (.not. allocated(chymrst%varid)) allocate(chymrst%varid(5))
!
      call nio_check(nf90_create(trim(dnres), nf90_clobber, chymrst%ncid))
!
!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_dim(chymrst%ncid, 'lon',                  &
                                  nlc, chymrst%dimid(1)))
      call nio_check(nf90_def_dim(chymrst%ncid, 'lat',                  &
                                  nbc, chymrst%dimid(2)))
      call nio_check(nf90_def_dim(chymrst%ncid, 'time',                 &
                                  nf90_unlimited, chymrst%dimid(3)))
!
!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymrst%ncid, 'lon', nf90_real,       &
                     chymrst%dimid(1), chymrst%varid(1)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1),       &
                     'long_name', 'Longitude'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1),       &
                     'units', 'degrees_east'))
!
      call nio_check(nf90_def_var(chymrst%ncid, 'lat', nf90_real,       &
                     chymrst%dimid(2), chymrst%varid(2)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(2),       &
                     'long_name', 'Latitude'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(2),       &
                     'units', 'degrees_north'))
!
      call nio_check(nf90_def_var(chymrst%ncid, 'time', nf90_int,       &
                     chymrst%dimid(3), chymrst%varid(3)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(3),       &
                     'long_name', 'Time'))
      write(str,fmt='("days since ",I4,"-",I2.2,"-",I2.2," 00:00:00")') &
           jahr1, jahr2, jahr3
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(3),       &
                     'units', str))

!
!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_def_var(chymrst%ncid, 'dis', nf90_real,       &
                     chymrst%dimid, chymrst%varid(4)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4),       &
                     'long_name', 'River Discharge'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4),       &
                     'missing_value', 1.0e20))
      call nio_check(nf90_def_var(chymrst%ncid, 'h2o', nf90_real,       &
                     chymrst%dimid, chymrst%varid(5)))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5),       &
                     'long_name', 'Total water'))
      call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5),       &
                     'missing_value', 1.0e20))
!
!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------
!
      call nio_check(nf90_enddef(chymrst%ncid))
!
!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(1), lon1))
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(2), lat1))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymrst%ncid))
!
      end subroutine chym_rst_init
!
      subroutine chym_out(istep)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: istep
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: len, start(4), count(4)
!
!-----------------------------------------------------------------------
!     Write data to file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_inquire_dimension(chymout%ncid,               &
                     chymout%dimid(3), len=len))
!
      start = (/ len+1, 1, 1, 1 /)
      count = (/ 1, 1, 1, 1 /)
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(3),       &
                    (/ istep-1 /), start, count))
!
      start = (/ 1, 1, len+1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_put_var(chymout%ncid, chymout%varid(4),       &
                     port, start, count))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymout%ncid))
!
!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------
!
      if (istep == nstep) call nio_check(nf90_close(chymout%ncid))
!
      end subroutine chym_out
!
      subroutine chym_rst(istep)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: istep
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: start(4), count(4), len
!
!-----------------------------------------------------------------------
!     Write data to file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_inquire_dimension(chymrst%ncid,               &
                     chymrst%dimid(3), len=len))
      start = (/ len+1, 1, 1, 1 /)
      count = (/ 1, 1, 1, 1 /)
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(3),       &
                    (/ istep-1 /), start, count))
!
      start = (/ 1, 1, len+1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(4),       &
                       port, start, count))
!
      start = (/ 1, 1, len+1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(5),       &
                       h2o, start, count))
!
!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_sync(chymrst%ncid))
!
!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------
!
      if (istep == nstep) call nio_check(nf90_close(chymrst%ncid))
!
      end subroutine chym_rst
!
      subroutine chym_ini()
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!
      integer :: ncid, varid, start(4), count(4)
!
!-----------------------------------------------------------------------
!     Open netCDF file
!-----------------------------------------------------------------------
!
      call nio_check(nf90_open(trim(dnini), nf90_nowrite, ncid))
!
!-----------------------------------------------------------------------
!     Read variables
!-----------------------------------------------------------------------
!
      call nio_check(nf90_inq_varid(ncid, 'dis', varid))
      start = (/ 1, 1, 1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_get_var(ncid, varid, port,                    &
                       start=start, count=count))
!
      call nio_check(nf90_inq_varid(ncid, 'h2o', varid))
      start = (/ 1, 1, 1, 1 /)
      count = (/ nlc, nbc, 1, 1 /)
      call nio_check(nf90_get_var(ncid, varid, h2o,                     &
                       start=start, count=count))
!
      call nio_check(nf90_inq_varid(ncid, 'time', varid))
      call nio_check(nf90_get_var(ncid, varid, pstep))
      pstep = pstep+1
!
!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------
!

      end subroutine chym_ini
!
      subroutine nio_check(status)
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
      use netcdf
!
      implicit none
!
!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------
!
      integer, intent(in) :: status

      if (status /= nf90_noerr) then
        print*, trim(nf90_strerror(status))
        stop 2
      end if
!
      end subroutine nio_check
!
      subroutine runoffspeed
      implicit none
      integer i,j,idir,iland
      real mann
      real alfamin,enne,gamma,delta,tresh,hrad
      alfa=0.0
      gamma=0.33
      delta=4.5                                       !cpar(8) in CHyM
      tresh=100.0                                     !cpar(6) in CHyM
      alfamin=0.1
      do j=2,chym_nlat-1
        do i=2,chym_nlon-1
          idir = int(fmap(i,j))
          iland = int(luse(i,j))
          mann = manning(iland)
          if ( idir >= 1 .and. idir <= 8 .and. &
              iland /= mare .and. iland > 0 ) then
            if (iland.gt.lntypes.or.iland.le.0) then
              print*,"Error in line: 845   in file: mod_chym_io"
              call exit(0)
            end if
            chym_dx(i,j)=geodistance(chym_lat(i,j),chym_lon(i,j),       &
                 chym_lat(i+ir(idir),j+jr(idir)),                       &
                 chym_lon(i+ir(idir),j+jr(idir)))
            if (chym_drai(i,j).gt.tresh) then
               enne=mann/delta
            else
               enne=mann/(1+(delta-1)*(1+(chym_drai(i,j)-tresh)/tresh))
            endif
            hrad=0.0015+0.050*((chym_drai(i,j)*1.e00)**gamma)
            !In CHyM 0.0015 = cpar( 2) ---> Alpha coefficients for
            !hydraulic radius (0.0015)
            !In CHyM 0.050 = cpar( 3) ---> Beta coefficients for
            !hydraulic radius (0.050)
            alfa(i,j)=((hrad**0.6666*accl(i,j)**0.5)/(enne))
            if (alfa(i,j).lt.alfamin) alfa(i,j)=alfamin
          endif
        enddo
      enddo
      end subroutine runoffspeed
!
      real function geodistance(latt1,lonn1,latt2,lonn2)
      implicit none
      real rad,dpi ; parameter(rad=6371000.0,dpi=6.2831855)
      real latt1,lonn1,latt2,lonn2,lt1,lt2,ln1,ln2,x,y
      lt1=latt1*dpi/360. ; lt2=latt2*dpi/360.
      ln1=lonn1*dpi/360. ; ln2=lonn2*dpi/360.
      if (abs(latt1-latt2).lt.0.2.and.abs(lonn1-lonn2).lt.0.2) then
          x=(rad*cos(lt1)*(ln1-ln2))*(rad*cos(lt2)*(ln1-ln2))
          y=(rad*(lt1-lt2))**2
          geodistance=sqrt(x+y)
      else
         x=sin(lt1)*sin(lt2)+cos(lt1)*cos(lt2)*cos((ln1)-(ln2))
         if (x.gt.1) x=1.0
         geodistance=acos(x)*rad
      endif
      if (geodistance.lt.0.1) geodistance=0.1
      return
      end function geodistance

      end module mod_chym_io

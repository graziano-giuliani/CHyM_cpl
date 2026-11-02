
!-----------------------------------------------------------------------
!     Module for CHyM model I/O subroutines
!-----------------------------------------------------------------------

module mod_chym_io

!-----------------------------------------------------------------------
!  Used module declarations
!-----------------------------------------------------------------------

   use, intrinsic :: iso_fortran_env
   use mpi
   use mod_chym_param

   implicit none
   private

   integer, public :: mycomm

   public :: read_config
   public :: read_init
   public :: chym_out_init
   public :: chym_out
   public :: chym_rst_init
   public :: chym_rst
   public :: chym_ini
   public :: chym_dispose

   contains

   subroutine chym_dispose
     use netcdf
     implicit none
     call nio_check(nf90_close(chymout%ncid),__LINE__)
     call nio_check(nf90_close(chymrst%ncid),__LINE__)
   end subroutine chym_dispose

   subroutine read_config(ifile)
    implicit none

!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------

    character(len=*) :: ifile

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

    integer :: iretval
    integer :: lun

    namelist /iniparam/ model_nsteps, efficiency, &
            nfarm, ifarm, irloss, irmonfac, riverrad, mouthsep
    namelist /inputparam/ isread, iswrit, nstep, &
            dnres, dnini, dnout, dnstt
    namelist /timeparam/ sdate, edate, calendar

!-----------------------------------------------------------------------
!     Read parameters
!-----------------------------------------------------------------------
    efficiency = 1.0
    nfarm = 0
    ifarm(:) = 0
    irloss = 0.0
    irmonfac(:) = 0.0
    mouthsep = 25.0

    open(newunit=lun, file=ifile, status='old', &
         action='read', iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'CHYM - Error opening namelist file'//trim(ifile)
      call mpi_abort(mycomm,1,iretval)
    end if
    model_nsteps = -1
    read(lun, nml=iniparam, iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'CHYM - Error reading iniparam namelist'
      call mpi_abort(mycomm,1,iretval)
    end if
    if ( model_nsteps < 0 ) then
      write(error_unit,*) 'CHYM - ERROR !'
      write(error_unit,*) 'CHYM - Model steps per day : ',model_nsteps
      write(error_unit,*) 'CHYM - ERROR !'
      call mpi_abort(mycomm,1,iretval)
    end if
    rewind(lun)
    read(lun, nml=inputparam, iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'CHYM - Error reading inputparam namelist'
      call mpi_abort(mycomm,1,iretval)
    end if
    rewind(lun)
    read(lun, nml=timeparam, iostat=iretval)
    if ( iretval /= 0 ) then
      write(error_unit,*) 'CHYM - Error reading timeparam namelist'
      call mpi_abort(mycomm,1,iretval)
    end if
    jahr1 = sdate/1000000
    jahr2 = (sdate-jahr1*1000000)/10000
    jahr3 = (sdate-(jahr1*1000000+jahr2*10000))/100
    jahr4 = sdate-(jahr1*1000000+jahr2*10000+jahr3*100)
    write(output_unit, *) "CHYM - Year ",jahr1,",month ",jahr2,",day ",jahr3
    close(lun)
  end subroutine read_config

  subroutine read_init()

    use netcdf
    implicit none

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

    integer :: i,j,ii,jj,iii,jjj,idir,ilnd
    integer :: ncid, varid, dimid, istat
    integer , dimension(2) :: idimid
    integer , dimension(3) :: ivarid
    integer :: contat

    call nio_check(nf90_open(trim(dnstt), nf90_nowrite, ncid),__LINE__)

    call nio_check(nf90_inq_dimid(ncid,'lon',dimid),__LINE__)
    call nio_check(nf90_inquire_dimension(ncid,dimid,len=chym_nlon),__LINE__)
    call nio_check(nf90_inq_dimid(ncid,'lat',dimid),__LINE__)
    call nio_check(nf90_inquire_dimension(ncid,dimid,len=chym_nlat),__LINE__)

    nbc = chym_nlat
    nlc = chym_nlon

!-----------------------------------------------------------------------
!     Allocate variables
!-----------------------------------------------------------------------

    if (.not. allocated(chym_area)) allocate(chym_area(nlc,nbc))
    if (.not. allocated(chym_drai)) allocate(chym_drai(nlc,nbc))
    if (.not. allocated(chym_lsm)) allocate(chym_lsm(nlc,nbc))
    if (.not. allocated(colx)) allocate(colx(nlc,nbc))
    if (.not. allocated(coly)) allocate(coly(nlc,nbc))
    if (.not. allocated(chym_dx)) allocate(chym_dx(nlc,nbc))
    if (.not. allocated(chym_lat)) allocate(chym_lat(nlc,nbc))
    if (.not. allocated(chym_mask)) allocate(chym_mask(nlc,nbc))
    if (.not. allocated(chym_lon)) allocate(chym_lon(nlc,nbc))
    if (.not. allocated(corner_lat)) allocate(corner_lat(4,nlc,nbc))
    if (.not. allocated(corner_lon)) allocate(corner_lon(4,nlc,nbc))
    if (.not. allocated(alfa)) allocate(alfa(nlc,nbc))
    if (.not. allocated(omdis)) allocate(omdis(nlc,nbc))
    if (.not. allocated(fmap)) allocate(fmap(nlc,nbc))
    if (.not. allocated(accl)) allocate(accl(nlc,nbc))
    if (.not. allocated(luse)) allocate(luse(nlc,nbc))
    if (.not. allocated(farm)) allocate(farm(nlc,nbc))
    if (.not. allocated(port)) allocate(port(nlc,nbc))
    if (.not. allocated(wkm1)) allocate(wkm1(nlc,nbc))
    if (.not. allocated(bwet)) allocate(bwet(nlc,nbc))
    if (.not. allocated(h2o)) allocate(h2o(nlc,nbc))
    if (.not. allocated(chym_runoff)) allocate(chym_runoff(nlc,nbc))
    if (.not. allocated(chym_dis)) allocate(chym_dis(nlc,nbc))
    if (.not. allocated(manning)) allocate(manning(lntypes))

    port = 0
    wkm1 = 0
    bwet = 0
    h2o = 0
    chym_lsm = 0.0
    chym_runoff = 0.0
    chym_dis = 0.0

    call nio_check(nf90_inq_varid(ncid, 'manning', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, manning),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'lon', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, chym_lon),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'corner_lon', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, corner_lon),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'lat', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, chym_lat),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'corner_lat', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, corner_lat),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'msk', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, chym_mask),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'fdm', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, fmap),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'acc', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, accl),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'lus', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, luse),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'aer', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, chym_area),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'dra', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, chym_drai),__LINE__)
    call nio_check(nf90_inq_varid(ncid, 'alf', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, alfa),__LINE__)
    istat = nf90_inq_varid(ncid, 'dis', varid)
    if ( istat == nf90_noerr ) then 
      call nio_check(nf90_get_var(ncid, varid, omdis),__LINE__)
    else
      omdis = riverrad
    end if
    call nio_check(nf90_close(ncid),__LINE__)

!-----------------------------------------------------------------------
!     Open, read and close restart file
!-----------------------------------------------------------------------

    pstep = 0
    if (isread /= 0) then
      call chym_ini()
      write(output_unit, *) "CHYM - Successfully read chym restart data"
    end if

    ! Every land cell draining into the coupled sea (chym_mask == 1) is a
    ! river mouth, whatever its drainage area: all of them deliver water.
    ! fdm = 9 marks an endorheic sink (no outflow).
    contat = 0
    chym_dx(:,:) = 0.1
    do j = 2, chym_nlat-1
      do i = 2, chym_nlon-1
        idir = fmap(i,j)
        if ( idir >= 1 .and. idir <= 8 ) then
          ii = i + ir(idir)
          jj = j + jr(idir)
          chym_dx(i,j) = geodistance(chym_lat(i,j),chym_lon(i,j),     &
                                     chym_lat(ii,jj),chym_lon(ii,jj))
          ilnd = luse(ii,jj)
          if ( luse(i,j) /= ocean .and. ilnd == ocean .and. &
               chym_mask(ii,jj) == 1 ) then
            if ( chym_lsm(ii,jj) <= 0.0 ) contat = contat + 1
            chym_lsm(ii,jj) = max(chym_lsm(ii,jj),omdis(i,j),1.0)
          end if
        end if
#ifdef NILE
        ! Booca 1 : Lat 31.52, Lon 31.84 (Damietta) 50%
        ! Bocca 2 : Lat 31.47, Lon 30.36 (Rosetta)  50%
        if ( is_inbox(lat_damietta,lon_damietta, &
               corner_lat(:,i,j),corner_lon(:,i,j)) ) then
          call find_nearest_land(i,j,iii,jjj)
          call find_nearest_ocean(iii,jjj,ii,jj)
          idamietta = ii
          jdamietta = jj
          if ( chym_lsm(ii,jj) < 0.5 ) then
            chym_lsm(ii,jj) = 200.0
            contat = contat + 1
          end if
          write(output_unit, *) 'CHYM - Damietta is at   ',ii,jj
        end if
        if ( is_inbox(lat_rosetta,lon_rosetta, &
               corner_lat(:,i,j),corner_lon(:,i,j)) ) then
          call find_nearest_land(i,j,iii,jjj)
          call find_nearest_ocean(iii,jjj,ii,jj)
          irosetta = ii
          jrosetta = jj
          if ( chym_lsm(ii,jj) < 0.5 ) then
            chym_lsm(ii,jj) = 200.0
            contat = contat + 1
          end if
          write(output_unit, *) 'CHYM - Rosetta is at    ',ii,jj
        end if
#endif
#ifdef BLACKSEA
        if (is_inbox(lat_dardanelli,lon_dardanelli, &
              corner_lat(:,i,j),corner_lon(:,i,j))) then
          call find_nearest_land(i,j,iii,jjj)
          call find_nearest_ocean(iii,jjj,ii,jj)
          idardanelli = ii
          jdardanelli = jj
          if ( chym_lsm(ii,jj) < 0.5 ) then
            chym_lsm(ii,jj) = 4000.0
            contat = contat + 1
          end if
          write(output_unit, *) 'CHYM - Dardanelli is at ',ii,jj
        end if
#endif
#ifdef AZOV
        if (is_inbox(lat_kerch,lon_kerch, &
              corner_lat(:,i,j),corner_lon(:,i,j))) then
          call find_nearest_land(i,j,iii,jjj)
          call find_nearest_ocean(iii,jjj,ii,jj)
          ikerch = ii
          jkerch = jj
          if ( chym_lsm(ii,jj) < 0.5 ) then
            chym_lsm(ii,jj) = max(omdis(i,j),1.0)
            contat = contat + 1
          end if
          write(output_unit, *) 'CHYM - Kerch is at      ',ii,jj
        end if
#endif
      end do
    end do

    call aggregate_mouths(contat)

    do j = 1, chym_nlat
      do i = 1, chym_nlon
        if ( any(ifarm(1:nfarm) == luse(i,j) ) ) then
          farm(i,j) = .true.
        else
          farm(i,j) = .false.
        end if
      end do
    end do

    call nio_check(nf90_create('rivermouth.nc', nf90_clobber,ncid),__LINE__)
    call nio_check(nf90_def_dim(ncid,'lon',nlc,idimid(1)),__LINE__)
    call nio_check(nf90_def_dim(ncid,'lat',nbc,idimid(2)),__LINE__)
    call nio_check(nf90_def_var(ncid,'lon',nf90_real,idimid,&
                          ivarid(1)),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(1), &
                          'standard_name','longitude'),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(1), &
                          'long_name','Longitude'),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(1), &
                          'units','degrees_east'),__LINE__)
    call nio_check(nf90_def_var(ncid,'lat', &
                          nf90_real,idimid,ivarid(2)),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(2), &
                          'standard_name','latitude'),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(2), &
                          'long_name','Latitude'),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(2), &
                          'units','degrees_north'),__LINE__)
    call nio_check(nf90_def_var(ncid,'lsm', &
                          nf90_real,idimid,ivarid(3)),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(3), &
                          'standard_name','binary_mask'),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(3), &
                          'long_name','River Mouths'),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(3),'units','1'),__LINE__)
    call nio_check(nf90_put_att(ncid,ivarid(3), &
                          'coordinates','lat lon'),__LINE__)
    call nio_check(nf90_enddef(ncid),__LINE__)
    call nio_check(nf90_put_var(ncid,ivarid(1), chym_lon),__LINE__)
    call nio_check(nf90_put_var(ncid,ivarid(2), chym_lat),__LINE__)
    call nio_check(nf90_put_var(ncid,ivarid(3), chym_lsm),__LINE__)
    call nio_check(nf90_close(ncid),__LINE__)
    write(output_unit, *) "CHYM - Diagnostic mouth position file created"
    write(output_unit, *) "CHYM - Total number of river mouths found : ",contat
  end subroutine read_init

  subroutine chym_out_init()

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

    character(len=100) :: str

!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------

    if (.not. allocated(chymout%dimid)) allocate(chymout%dimid(3))
    if (.not. allocated(chymout%varid)) allocate(chymout%varid(5))

    call nio_check(nf90_create(trim(dnout), nf90_clobber, &
                               chymout%ncid),__LINE__)

!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------

    call nio_check(nf90_def_dim(chymout%ncid, 'lon',                  &
                                nlc, chymout%dimid(1)),__LINE__)
    call nio_check(nf90_def_dim(chymout%ncid, 'lat',                  &
                                nbc, chymout%dimid(2)),__LINE__)
    call nio_check(nf90_def_dim(chymout%ncid, 'time',                 &
                                nf90_unlimited, chymout%dimid(3)),__LINE__)

!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------

    call nio_check(nf90_def_var(chymout%ncid, 'lon', nf90_real,       &
                     chymout%dimid(1:2), chymout%varid(1)),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                     'long_name', 'Longitude'),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(1),       &
                     'units', 'degrees_east'),__LINE__)

    call nio_check(nf90_def_var(chymout%ncid, 'lat', nf90_real,       &
                     chymout%dimid(1:2), chymout%varid(2)),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                     'long_name', 'Latitude'),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(2),       &
                     'units', 'degrees_north'),__LINE__)

    call nio_check(nf90_def_var(chymout%ncid, 'time', nf90_int,       &
                     chymout%dimid(3), chymout%varid(3)),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                     'long_name', 'Time'),__LINE__)
    write(str,fmt='("days since ",I4,"-",I2.2,"-",I2.2," 00:00:00")') &
           jahr1, jahr2, jahr3
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(3),       &
                     'units', str),__LINE__)

!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------

    call nio_check(nf90_def_var(chymout%ncid, 'dis', nf90_real,       &
                     chymout%dimid, chymout%varid(4)),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'long_name', 'River Discharge'),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'missing_value', 1.0e20),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(4),       &
                     'coordinates', 'lat lon'),__LINE__)
    call nio_check(nf90_def_var(chymout%ncid, 'ro', nf90_real,        &
                     chymout%dimid, chymout%varid(5)),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(5),       &
                     'long_name', 'Atmosphere Runoff'),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(5),       &
                     'missing_value', 1.0e20),__LINE__)
    call nio_check(nf90_put_att(chymout%ncid, chymout%varid(5),       &
                     'coordinates', 'lat lon'),__LINE__)

!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------

    call nio_check(nf90_enddef(chymout%ncid),__LINE__)

!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------

    call nio_check(nf90_put_var(chymout%ncid, &
                           chymout%varid(1), chym_lon),__LINE__)
    call nio_check(nf90_put_var(chymout%ncid, &
                           chymout%varid(2), chym_lat),__LINE__)

!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------

    call nio_check(nf90_sync(chymout%ncid),__LINE__)
    chymout%nrec = 0

  end subroutine chym_out_init

  subroutine chym_rst_init()

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

    character(len=100) :: str

!-----------------------------------------------------------------------
!     Create netCDF file
!-----------------------------------------------------------------------

    if (.not. allocated(chymrst%dimid)) allocate(chymrst%dimid(3))
    if (.not. allocated(chymrst%varid)) allocate(chymrst%varid(6))

    call nio_check(nf90_create(trim(dnres), nf90_clobber, &
                               chymrst%ncid),__LINE__)

!-----------------------------------------------------------------------
!     Define dimensions
!-----------------------------------------------------------------------

    call nio_check(nf90_def_dim(chymrst%ncid, 'lon',                  &
                                nlc, chymrst%dimid(1)),__LINE__)
    call nio_check(nf90_def_dim(chymrst%ncid, 'lat',                  &
                                nbc, chymrst%dimid(2)),__LINE__)
    call nio_check(nf90_def_dim(chymrst%ncid, 'time',                 &
                                nf90_unlimited, chymrst%dimid(3)),    &
                                __LINE__)

!-----------------------------------------------------------------------
!     Define dimension variables
!-----------------------------------------------------------------------

    call nio_check(nf90_def_var(chymrst%ncid, 'lon', nf90_real,       &
                     chymrst%dimid(1:2), chymrst%varid(1)),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1),       &
                     'long_name', 'Longitude'),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(1),       &
                     'units', 'degrees_east'),__LINE__)

    call nio_check(nf90_def_var(chymrst%ncid, 'lat', nf90_real,       &
                     chymrst%dimid(1:2), chymrst%varid(2)),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(2),       &
                     'long_name', 'Latitude'),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(2),       &
                     'units', 'degrees_north'),__LINE__)

    call nio_check(nf90_def_var(chymrst%ncid, 'time', nf90_int,       &
                     chymrst%dimid(3), chymrst%varid(3)),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(3),       &
                     'long_name', 'Time'),__LINE__)
    write(str,fmt='("days since ",I4,"-",I2.2,"-",I2.2," 00:00:00")') &
           jahr1, jahr2, jahr3
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(3),       &
                     'units', str),__LINE__)

!-----------------------------------------------------------------------
!     Define variables
!-----------------------------------------------------------------------

    call nio_check(nf90_def_var(chymrst%ncid, 'dis', nf90_real,       &
                     chymrst%dimid, chymrst%varid(4)),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4),       &
                     'long_name', 'River Discharge'),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4),       &
                     'missing_value', 1.0e20),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(4),       &
                     'coordinates', 'lat lon'),__LINE__)
    call nio_check(nf90_def_var(chymrst%ncid, 'h2o', nf90_real,       &
                     chymrst%dimid, chymrst%varid(5)),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5),       &
                     'long_name', 'Total water'),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5),       &
                     'missing_value', 1.0e20),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(5),       &
                     'coordinates', 'lat lon'),__LINE__)
    call nio_check(nf90_def_var(chymrst%ncid, 'rno', nf90_real,       &
                     chymrst%dimid, chymrst%varid(6)),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(6),       &
                     'long_name', 'Total runoff'),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(6),       &
                     'missing_value', 1.0e20),__LINE__)
    call nio_check(nf90_put_att(chymrst%ncid, chymrst%varid(6),       &
                     'coordinates', 'lat lon'),__LINE__)

!-----------------------------------------------------------------------
!     Exit define mode
!-----------------------------------------------------------------------

    call nio_check(nf90_enddef(chymrst%ncid),__LINE__)

!-----------------------------------------------------------------------
!     Fill coordinate variables
!-----------------------------------------------------------------------

    call nio_check(nf90_put_var(chymrst%ncid, &
                             chymrst%varid(1), chym_lon),__LINE__)
    call nio_check(nf90_put_var(chymrst%ncid, &
                             chymrst%varid(2), chym_lat),__LINE__)

!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------

    call nio_check(nf90_sync(chymrst%ncid),__LINE__)
    chymrst%nrec = 0

  end subroutine chym_rst_init

  subroutine chym_out(istep)

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------

    integer, intent(in) :: istep

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

    integer :: ilen, start(4), count(4)

!-----------------------------------------------------------------------
!     Write data to file
!-----------------------------------------------------------------------

    ilen = chymout%nrec + 1
    start = (/ ilen, 1, 1, 1 /)
    count = (/ 1, 1, 1, 1 /)
    call nio_check(nf90_put_var(chymout%ncid, chymout%varid(3),       &
                    (/ istep-1 /), start, count),__LINE__)

    start = (/ 1, 1, ilen, 1 /)
    count = (/ nlc, nbc, 1, 1 /)
    call nio_check(nf90_put_var(chymout%ncid, chymout%varid(4),       &
                     port, start, count),__LINE__)
    call nio_check(nf90_put_var(chymout%ncid, chymout%varid(5),       &
                     chym_runoff, start, count),__LINE__)

!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------

    call nio_check(nf90_sync(chymout%ncid),__LINE__)
    chymout%nrec = ilen
  end subroutine chym_out

  subroutine chym_rst(istep)

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------

    integer, intent(in) :: istep

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

    integer :: start(4), count(4), ilen

!-----------------------------------------------------------------------
!     Write data to file
!-----------------------------------------------------------------------

    ilen = chymrst%nrec + 1
    start = (/ ilen, 1, 1, 1 /)
    count = (/ 1, 1, 1, 1 /)
    call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(3),  &
                    (/ istep-1 /), start, count),__LINE__)

    start = (/ 1, 1, ilen, 1 /)
    count = (/ nlc, nbc, 1, 1 /)
    call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(4),  &
                       port, start, count),__LINE__)
    call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(5),  &
                       h2o, start, count),__LINE__)
    call nio_check(nf90_put_var(chymrst%ncid, chymrst%varid(6),  &
                       chym_runoff, start, count),__LINE__)

!-----------------------------------------------------------------------
!     Sync file
!-----------------------------------------------------------------------

    call nio_check(nf90_sync(chymrst%ncid),__LINE__)
    chymrst%nrec = ilen
  end subroutine chym_rst

  subroutine chym_ini()

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------

    integer :: ncid, ncstat, dimid, varid, nt, start(3), count(3)
    integer :: last_time(1)

!-----------------------------------------------------------------------
!     Open netCDF file
!-----------------------------------------------------------------------

    call nio_check(nf90_open(trim(dnini), nf90_nowrite, ncid),__LINE__)

!-----------------------------------------------------------------------
!   Select last timestep
!-----------------------------------------------------------------------

    call nio_check(nf90_inq_dimid(ncid,'time',dimid),__LINE__)
    call nio_check(nf90_inquire_dimension(ncid,dimid,len=nt),__LINE__)

!-----------------------------------------------------------------------
!     Read variables
!-----------------------------------------------------------------------

    write(output_unit, *) "CHYM - Reading timestep ",nt
    start = (/   1,   1, nt /)
    count = (/ nlc, nbc,  1 /)

    call nio_check(nf90_inq_varid(ncid, 'dis', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, port,                    &
                       start=start, count=count),__LINE__)

    call nio_check(nf90_inq_varid(ncid, 'h2o', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, h2o,                     &
                       start=start, count=count),__LINE__)

    write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Discharge max value: ", &
                             maxval(port)
    write(output_unit,fmt='(1x,A,F16.2)')"CHYM - H2o max value      : ", &
                             maxval(h2o)

    ncstat = nf90_inq_varid(ncid, 'rno', varid)
    if ( ncstat == nf90_noerr ) then
      write(output_unit, *) "CHYM - Runoff present, reading it"
      call nio_check(nf90_get_var(ncid, varid, chym_runoff,             &
                         start=start, count=count),__LINE__)
      write(output_unit,fmt='(1x,A,F16.2)')"CHYM - Runoff max value   : ", &
                          maxval(chym_runoff)*86400*1000.0, ' mm/day'
    end if

    start(1) = nt
    count(1) = 1
    write(output_unit, *) "CHYM - Reading time information..."
    call nio_check(nf90_inq_varid(ncid, 'time', varid),__LINE__)
    call nio_check(nf90_get_var(ncid, varid, last_time, &
                         start=start(1:1), count=count(1:1)),__LINE__)
    pstep = last_time(1)+1
    write(output_unit, *) "CHYM - last timestep is ",last_time(1)

!-----------------------------------------------------------------------
!     Close file
!-----------------------------------------------------------------------

    call nio_check(nf90_close(ncid),__LINE__)

  end subroutine chym_ini

  subroutine nio_check(istatus,line)

!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------

    use netcdf

    implicit none

!-----------------------------------------------------------------------
!     Imported variable declarations
!-----------------------------------------------------------------------

    integer, intent(in) :: istatus
    integer, intent(in) :: line
    integer :: ierr

    if (istatus /= nf90_noerr) then
      write(error_unit, *) 'CHYM - At line ',line
      write(error_unit, *) trim(nf90_strerror(istatus))
      call mpi_abort(mycomm,2,ierr)
    end if

  end subroutine nio_check

  ! Merge river mouths closer than mouthsep (km) into outlets (collectors).
  ! The largest mouths become collectors first; prescribed outlets (Nile,
  ! Dardanelles, Kerch) always do. Every other mouth sends its water to the
  ! nearest collector, which is less than mouthsep away by construction.
  ! The ocean coupler sees only the collectors, all at least mouthsep apart,
  ! so RegESM's isolated-point river detection keeps every one of them.
  subroutine aggregate_mouths(nout)
    implicit none
    integer, intent(out) :: nout
    integer :: nm, n, k, l, i, j, nbest
    integer, allocatable :: mi(:), mj(:), ord(:), col(:), tgt(:)
    real, allocatable :: w(:), key(:)
    logical, allocatable :: forced(:), done(:)
    real :: d, dbest

    nm = count(chym_lsm > 0.0)
    colx = 0
    coly = 0
    nout = 0
    if ( nm == 0 ) return
    allocate(mi(nm), mj(nm), w(nm), key(nm), forced(nm), done(nm))
    allocate(ord(nm), col(nm), tgt(nm))
    n = 0
    do j = 1, chym_nlat
      do i = 1, chym_nlon
        if ( chym_lsm(i,j) > 0.0 ) then
          n = n + 1
          mi(n) = i
          mj(n) = j
          w(n) = chym_lsm(i,j)
          forced(n) = is_prescribed(i,j)
        end if
      end do
    end do
    ! Processing order: prescribed outlets, then decreasing weight
    key = w
    where ( forced ) key = huge(1.0)
    done = .false.
    do k = 1, nm
      n = maxloc(key, dim=1, mask=.not. done)
      ord(k) = n
      done(n) = .true.
    end do
    ! Greedy selection of the collectors
    nout = 0
    do k = 1, nm
      n = ord(k)
      if ( .not. forced(n) ) then
        dbest = huge(1.0)
        do l = 1, nout
          d = geodistance(chym_lat(mi(n),mj(n)),chym_lon(mi(n),mj(n)), &
                chym_lat(mi(col(l)),mj(col(l))),chym_lon(mi(col(l)),mj(col(l))))
          dbest = min(dbest, d)
          if ( dbest <= mouthsep*1000.0 ) exit
        end do
        if ( dbest <= mouthsep*1000.0 ) cycle
      end if
      nout = nout + 1
      col(nout) = n
    end do
    ! Assign every mouth to its nearest collector
    do n = 1, nm
      dbest = huge(1.0)
      nbest = 1
      do l = 1, nout
        d = geodistance(chym_lat(mi(n),mj(n)),chym_lon(mi(n),mj(n)), &
              chym_lat(mi(col(l)),mj(col(l))),chym_lon(mi(col(l)),mj(col(l))))
        if ( d < dbest ) then
          dbest = d
          nbest = l
        end if
      end do
      tgt(n) = col(nbest)
    end do
    ! Mouths that are collectors map onto themselves (also the forced ones)
    do l = 1, nout
      tgt(col(l)) = col(l)
    end do
    ! New outlet map: weight of a collector = sum of its members' weights
    chym_lsm = 0.0
    do n = 1, nm
      colx(mi(n),mj(n)) = mi(tgt(n))
      coly(mi(n),mj(n)) = mj(tgt(n))
      chym_lsm(mi(tgt(n)),mj(tgt(n))) = chym_lsm(mi(tgt(n)),mj(tgt(n))) + w(n)
    end do
    write(output_unit,'(1x,a,i6,a,i6,a,f6.1,a)') 'CHYM - ', nm, &
      ' river mouth cells aggregated into ', nout, &
      ' outlets (minimum separation ', mouthsep, ' km)'
    deallocate(mi, mj, w, key, forced, done, ord, col, tgt)
  end subroutine aggregate_mouths

  logical function is_prescribed(i,j)
    implicit none
    integer, intent(in) :: i, j
    is_prescribed = .false.
#ifdef NILE
    if ( i == idamietta .and. j == jdamietta ) is_prescribed = .true.
    if ( i == irosetta .and. j == jrosetta ) is_prescribed = .true.
#endif
#ifdef BLACKSEA
    if ( i == idardanelli .and. j == jdardanelli ) is_prescribed = .true.
#endif
#ifdef AZOV
    if ( i == ikerch .and. j == jkerch ) is_prescribed = .true.
#endif
  end function is_prescribed

  real function geodistance(latt1,lonn1,latt2,lonn2)
    implicit none
    real, parameter :: rad = 6371000.0
    real, parameter :: dpi = 6.2831855
    real :: latt1,lonn1,latt2,lonn2,lt1,lt2,ln1,ln2,x,y
    lt1 = latt1*dpi/360.
    lt2 = latt2*dpi/360.
    ln1 = lonn1*dpi/360.
    ln2 = lonn2*dpi/360.
    if ( abs(latt1-latt2) < 0.2 .and. &
         abs(lonn1-lonn2) < 0.2 ) then
      x = (rad*cos(lt1)*(ln1-ln2))*(rad*cos(lt2)*(ln1-ln2))
      y = (rad*(lt1-lt2))**2
      geodistance = sqrt(x+y)
    else
      x = sin(lt1)*sin(lt2)+cos(lt1)*cos(lt2)*cos((ln1)-(ln2))
      if ( x > 1.0 ) x = 1.0
      geodistance = acos(x)*rad
    endif
    if ( geodistance < 0.1 ) geodistance = 0.1
  end function geodistance

end module mod_chym_io

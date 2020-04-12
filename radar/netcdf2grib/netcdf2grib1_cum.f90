program netcdf2grib1_cum

!================================================================
!
! Trasforma i file netdcf contenenti la precipitazione cumulata 
! da radar (cosi' come codificati ad ARPA-SIMC) in file GRIB1.
!
!================================================================

use netcdf
use char_utilities
use datetime_class
use missing_values

use log4fortran
use gridinfo_class
use grid_class
use grid_id_class

use vol7d_level_class
use vol7d_timerange_class
use volgrid6d_var_class
use err_handling

implicit none


character(LEN=100):: name_nc, delta

integer :: il, iun, ier, iargc, n, i, j, ii, jj

integer :: ncid, ndim, nvar, natt, nunlimdimid
integer :: id_lon, id_lat, id_geo, id_mesh, id_time
integer :: dim_lon, dim_lat, dim_geo, dim_mesh
integer :: new_lat, new_lon
integer :: varid, varid_time, date(3), hour(2), delta_t, dim_t
integer :: varid_pr, varid_geo, varid_mesh 

character(len=100) :: time, dum
character(len=100) :: data_char,ora_char
character(len=100) :: varname
character(len=10)  :: pr_units
!character(:),allocatable :: unpr !funziona solo nuove versioni di gfortran
real :: acc_t

real, allocatable :: cum_pr_mm(:,:,:)
doubleprecision, allocatable :: geo_lim(:), mesh_xy(:)

! Var per il logging
integer :: category 
character(len=512):: a_name

! Var per scrivere il grib
character(len=80)::fileout
type(gridinfo_def) :: gridinfo
type(arrayof_gridinfo) :: gridinfov
type(griddim_def) :: griddim
integer :: nx, ny 
integer,parameter:: component_flag=0  
type(grid_id) :: gaid_template
type(vol7d_level) :: level
type(vol7d_timerange) :: timerange
type(volgrid6d_var) :: var
type(datetime) :: date_time
doubleprecision :: xmin, xmax, ymin, ymax
character(len=80) :: tipo='regular_ll'

!-------------------------------
! get launcher name
call l4f_launcher(a_name,a_name_force="netcdf2grib1_cum")
! log4fortran init
ier=l4f_init()
! set a_name
category=l4f_category_get(a_name//".main")
call l4f_category_log(category,L4F_INFO,"Inizio")

! Lettura del numero degli argomenti esterni
n=iargc() 

!-------------------------------
! Parametri in input: nome file
if (iargc() < 1) then
  call l4f_category_log(category,L4F_ERROR,&
   "Mancano gli argomenti. Uso: netcdf2grib1_cum filencdf")
  call raise_error() 
endif

! L'argomento passato e' il nome del file
call getarg(1,name_nc)

!------------------------------------------------------------------
! Lettura dal file NETCDF delle variabili necessarie per la 
! scrittura del file GRIB.

! Apertura del file netcdf 
call check( nf90_open(name_nc,NF90_NOWRITE,ncid),"Apertura" )

! Estraggo le dimensioni delle variabili immagazzinate
call check( nf90_inq_dimid(ncid,"lon" ,id_lon ) )
call check( nf90_inq_dimid(ncid,"lat" ,id_lat ) )
call check( nf90_inq_dimid(ncid,"geo_dim" ,id_geo  ) )
call check( nf90_inq_dimid(ncid,"mesh_dim",id_mesh ) )
call check( nf90_inquire_dimension(ncid,id_lon ,len=dim_lon  ) )
call check( nf90_inquire_dimension(ncid,id_lat ,len=dim_lat  ) )
call check( nf90_inquire_dimension(ncid,id_geo ,len=dim_geo  ) )
call check( nf90_inquire_dimension(ncid,id_mesh,len=dim_mesh ) )

! Estraggo l'istante di emissione del dato
call check( nf90_inq_dimid(ncid,"time",id_time) )
call check( nf90_inq_varid(ncid,"time",varid_time) )
call check( nf90_get_att(ncid, varid_time,"units",time) )
call check( nf90_inquire_dimension(ncid,id_time,len=dim_t) )

! Time e' definito come "hour before AAAA-MM-GG hh:mm:0"
! Nelle vecchie cumulate potrebbe esserci scritto "hours"
! invece di "hour", "since" invece "before". Per ovviare
! al problema divido la stringa in 4 sottostringhe.
! In questo modo la data e' nella sottostringa 3 e l'ora
! nella sottostringa 4

! Controllo che l'unita' di cumulazione non siano i minuti.
if (time(1:7) == 'minutes') then 
  call l4f_category_log(category,L4F_ERROR, &
  "Le cumulate in minuti non sono ancora gestite, esco")
  call raise_error()   
else if (time(1:4) == 'hour') then
  dum=time
  do i=2,4
    j=index(dum,' ')
    dum=dum(j+1:len(dum))
    j=j+1
    if (i == 3) then
      read(dum(  1: 4 ),'(I4)') date(1)
      read(dum(  6: 7 ),'(I2)') date(2)
      read(dum(  9:10 ),'(I2)') date(3)
    else if (i == 4) then
      read(dum( 1:2 ),'(I2)') hour(1)
      read(dum( 4:6 ),'(I2)') hour(2)
    endif
  enddo
endif

! Trasformo data e ora in stringa
data_char=trim(to_char(date(1),form='(I4)'))//   &
 trim(to_char(date(2),form='(I2.2)'))// &
 trim(to_char(date(3),form='(I2.2)'))
ora_char =trim(to_char(hour(1),form='(I2.2)'))// &
 trim(to_char(hour(2),form='(I2.2)'))

! Estraggo le variabili necessarie al grib
call check( nf90_inquire(ncid, ndim, nvar, natt, nunlimdimid) )

do varid=1,nvar
  call check( nf90_inquire_variable(ncid, varid, varname) )
  if ( varname == 'cum_pr_mm') then
    varid_pr=varid
  else if ( varname == 'geo_dim') then
    varid_geo=varid
  else if ( varname == 'mesh_dim') then
    varid_mesh=varid
  endif  
enddo

allocate ( geo_lim(dim_geo ) )
allocate ( mesh_xy(dim_mesh) )
allocate ( cum_pr_mm(dim_lon,dim_lat,dim_t) )

call check( nf90_get_var(ncid,varid_geo ,geo_lim) )
call check( nf90_get_var(ncid,varid_mesh,mesh_xy) )
call check( nf90_get_var(ncid,varid_pr  ,cum_pr_mm) )
! Estraggo gli attributi del campo di pioggia 
call check( nf90_get_att(ncid,varid_pr,"units"       ,pr_units) )
call check( nf90_get_att(ncid,varid_pr,"accum_time_h",acc_t   ) )

! Verifico che l'unita' di misura sia 'mm':
! Per una qualche ragione non chiara c'è un carattere nullo nella lettura
! di pr_units che non viene eliminato dalla funzione TRIM. Elimino il
! problema con un trucco
!allocate(character(len(trim(pr_units))-1) :: unpr) 
!unpr=pr_units(1:len(unpr))
if (pr_units(1:len(trim(pr_units))-1) /= 'mm') then
  call l4f_category_log(category,L4F_ERROR, &
  "L'unita' di misura non e' mm, esco")
  call raise_error()
endif  

if (acc_t == 0.) then
  call l4f_category_log(category,L4F_WARN, &
   "Accumulation time (acc_t) not defined! Default= 1.0 hour")
  acc_t= 1.0
endif

! Sostituisco il valor mancante con rmiss
where ( cum_pr_mm < 0. )
  cum_pr_mm=rmiss
end where

!==============================================================================
!  SCRITTURA DEL GRIB
!==============================================================================
fileout="radar_SRT_"//trim(data_char)//trim(ora_char)//"_"//&
 trim(to_char(int(acc_t)))//"h.grib1"

call l4f_category_log(category,L4F_INFO,"Output file= "//trim(fileout))

! Definizione della griglia
call init(griddim,proj_type=tipo,nx=dim_lon,ny=dim_lat, &
  xmin=geo_lim(2),xmax=geo_lim(4),ymin=geo_lim(1),ymax=geo_lim(3), &
  component_flag=component_flag,categoryappend="generated")

gaid_template=grid_id_new(grib_api_template="regular_ll_sfc_grib1")

! Data di emissione
call init(date_time,year=date(1),month=date(2),day=date(3),&
  hour=hour(1),minute=hour(2),msec=00)
! Timerange locale SIMC di analisi di precipitazione
call init(timerange,timerange=1,p1=0,p2=int(acc_t*60))
! Livelli 
call init(level, level1=1, l1=0, level2=imiss, l2=imiss)
! Variabile precipitazione 
call init(var,centre=200, number=61,category=2)

call init (gridinfo,gaid_template,griddim,date_time, &
 timerange,level,var,clone=.false.)

call grib_set(grid_id_get_gaid(gridinfo%gaid),"centre",80)
call grib_set(grid_id_get_gaid(gridinfo%gaid),"generatingProcessIdentifier",1)
call grib_set(grid_id_get_gaid(gridinfo%gaid),'missingValue',rmiss)
call grib_set(grid_id_get_gaid(gridinfo%gaid),'packingType','grid_simple')
call grib_set(grid_id_get_gaid(gridinfo%gaid),"bitmapPresent",1)

call encode_gridinfo (gridinfo,cum_pr_mm(:,:,1))
!call display(gridinfo)

call l4f_category_log(category,L4F_INFO,"export to GRIB")

call insert(gridinfov, gridinfo)
call export(gridinfov,filename=fileout,categoryappend="gridinfo scritto")

call l4f_category_log(category,L4F_INFO,"end")

call delete (gridinfov)
! close logger
call l4f_category_delete(category)
ier=l4f_fini()

deallocate (cum_pr_mm)
deallocate (geo_lim)
deallocate (mesh_xy)

!=========================
! CHECK 
!=========================
contains
subroutine check(status,message)
integer, intent(in) :: status
character(len=*),optional :: message    

if (status /= nf90_noerr) then
  call l4f_category_log(category,L4F_ERROR,&
   nf90_strerror(status)//optio_c(message,40))
  call raise_error()
end if
end subroutine check


end program netcdf2grib1_cum

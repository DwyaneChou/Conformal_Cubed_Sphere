module output_mod
  use netcdf
  use constants_mod
  use parameters_mod
  use mesh_mod
  implicit none
    
    character(13) :: ncFile = 'ccs_output.nc'
    
    contains
    subroutine history_init
      
      integer status
      integer ncid
      integer lon_dim_id,lat_dim_id
      integer patch_dim_id
      integer x_id,y_id
      integer lon_id,lat_id
      
      status = nf90_create(ncFile, NF90_CLOBBER + NF90_NETCDF4 , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'lon'   ,Nx,lon_dim_id )
      status = nf90_def_dim(ncid,'lat'   ,Ny,lat_dim_id )
      status = nf90_def_dim(ncid,'nPatch',Nf,patch_dim_id)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'x'  ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),x_id       )
      status = nf90_def_var(ncid,'y'  ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),y_id       )
      status = nf90_def_var(ncid,'lon',NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),lon_id     )
      status = nf90_def_var(ncid,'lat',NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),lat_id     )
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_att'
      status = nf90_put_att(ncid,nf90_global       ,'dx'       ,dx)
      status = nf90_put_att(ncid,nf90_global       ,'dy'       ,dy)
      status = nf90_put_att(ncid,nf90_global       ,'ids'      ,ids)
      status = nf90_put_att(ncid,nf90_global       ,'ide'      ,ide)
      status = nf90_put_att(ncid,nf90_global       ,'jds'      ,jds)
      status = nf90_put_att(ncid,nf90_global       ,'jde'      ,jde)
      status = nf90_put_att(ncid,nf90_global       ,'ifs'      ,ifs)
      status = nf90_put_att(ncid,nf90_global       ,'ife'      ,ife)
      
      status = nf90_put_att(ncid,lon_id            ,'units'    ,'degree_east' )
      status = nf90_put_att(ncid,lat_id            ,'units'    ,'degree_north')
      
      status = nf90_put_att(ncid,lon_id            ,'long_name','longitude on sphere coordinate for Cells' )
      status = nf90_put_att(ncid,lat_id            ,'long_name','latitude on sphere coordinate for Cells'  )
      
      ! Define coordinates
      status = nf90_put_att(ncid, x_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, y_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lon_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lat_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,x_id       , mesh%x        )
      status = nf90_put_var(ncid,y_id       , mesh%y        )
      status = nf90_put_var(ncid,lon_id     , mesh%lon * R2D)
      status = nf90_put_var(ncid,lat_id     , mesh%lat * R2D)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_init
    
    subroutine handle_err(status)
      implicit none
      integer,intent(in)::status
            
      if(status/=nf90_noerr)then
          print*, trim(nf90_strerror(status))
          stop "Stopped by netCDF"
      endif  
    endsubroutine handle_err
end module output_mod
    
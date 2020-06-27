
    MODULE mesh_mod
      use constants_mod
      use parameters_mod
      use cuco
      use ghost_mod
      use math_mod
      use jimcc_m
      
      implicit none
      
      ! coordinate
      type mesh_info
        real, dimension(:,:,:    ), allocatable :: xi       ! x-dir local coordinate on map panel
        real, dimension(:,:,:    ), allocatable :: eta      ! y-dir local coordinate on map panel
        real, dimension(:,:,:    ), allocatable :: x        ! cartesian coordinate x
        real, dimension(:,:,:    ), allocatable :: y        ! cartesian coordinate y
        real, dimension(:,:,:    ), allocatable :: z        ! cartesian coordinate z
        real, dimension(:,:,:    ), allocatable :: lon      ! longitude on cells
        real, dimension(:,:,:    ), allocatable :: lat      ! latitude on cells
        real, dimension(:,:,:,:,:), allocatable :: jab      ! jacobian matrix
        
      end type mesh_info
      
      type(mesh_info), target :: mesh
      
      contains
      
      subroutine initMesh
        integer :: iPV, jPV, iCell, jCell, iPatch, iVar, iDOF, jDOF
        integer :: iPVs, iPVe, jPVs, jPVe
        
        ! Allocate arrays in structures
        allocate( mesh%xi (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%eta(      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%x  (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%y  (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%z  (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lon(      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lat(      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%jab(3, 2, ids:ide, ids:ide, ifs:ife) )
        
        ! Calculate mesh infomation on VIA
        do iPatch = ifs, ife
          do jCell = jds, jde
            do iCell = ids, ide
              mesh%xi (iCell, jCell, iPatch) = ( (iCell - 1) * dx / 2. + x_min ) * D2M
              mesh%eta(iCell, jCell, iPatch) = ( (jCell - 1) * dy / 2. + y_min ) * D2M
            end do
          end do
        end do
        
      end subroutine initMesh
      
      subroutine generateConformalCubedSphere
      
        real xx,yy
        real lon,lat
        real dlondx,dlatdx
        real dlondy,dlatdy
        real dlondz,dlatdz
        real sph2cart_matrix(2,3)
        
        real dxdxi ,dydxi ,dzdxi
        real dxdeta,dydeta,dzdeta
        real dxdr  ,dydr  ,dzdr
      
        real,dimension(2):: xm
        real,dimension(3):: xc
        
        integer, parameter :: stencil_width = 2 * xhalo+1
        real fc(stencil_width)
        real dh
        
        real, dimension(:,:,:), allocatable :: x_tmp
        real, dimension(:,:,:), allocatable :: y_tmp
        real, dimension(:,:,:), allocatable :: z_tmp
        real, dimension(:,:,:), allocatable :: lon_tmp
        real, dimension(:,:,:), allocatable :: lat_tmp
        
        integer i,j,k
        integer iPatch
        
        allocate(x_tmp  (ips:ipe,jps:jpe,ifs:ife))
        allocate(y_tmp  (ips:ipe,jps:jpe,ifs:ife))
        allocate(z_tmp  (ips:ipe,jps:jpe,ifs:ife))
        allocate(lon_tmp(ips:ipe,jps:jpe,ifs:ife))
        allocate(lat_tmp(ips:ipe,jps:jpe,ifs:ife))
        
        ! Generating grids
        if(trim(adjustl(gen_method))=='CCAM')then
          call inrot()
          do k = ifs,ife
            !$OMP PARALLEL DO PRIVATE(i,xx,yy,xc,lon,lat,iPatch)
            do j = jds, jde
              do i = ids, ide
                xx = mesh%xi (i,j,k) / 2. + 0.5
                yy = mesh%eta(i,j,k) / 2. + 0.5
                
                call mtoc(xx, yy, k, xc)
                call cart2sph(lon,lat,xc(1),xc(2),xc(3))
                
                ! reset longitude to [0,2*pi)
                if(lon<0    ) lon = 2. * pi + lon
                if(lon>2.*pi) lon = lon - 2. * pi
                
                if(k == 1) iPatch = 4
                if(k == 2) iPatch = 1
                if(k == 3) iPatch = 2
                if(k == 4) iPatch = 3
                if(k == 5) iPatch = 5
                if(k == 6) iPatch = 6
                !iPatch = k
                
                mesh%x  (i,j,iPatch) = xc(1)
                mesh%y  (i,j,iPatch) = xc(2)
                mesh%z  (i,j,iPatch) = xc(3)
                mesh%lon(i,j,iPatch) = lon
                mesh%lat(i,j,iPatch) = lat
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
          
          ! rotate relative position in panel 5 and 6
          x_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%x
          y_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%y
          z_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%z
          lon_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lon
          lat_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lat
          
          k = 5
          do i = ids,ide
            mesh%x  (ids:ide,i,k) = x_tmp  (ide-i+1,jds:jde,k)
            mesh%y  (ids:ide,i,k) = y_tmp  (ide-i+1,jds:jde,k)
            mesh%z  (ids:ide,i,k) = z_tmp  (ide-i+1,jds:jde,k)
            mesh%lon(ids:ide,i,k) = lon_tmp(ide-i+1,jds:jde,k)
            mesh%lat(ids:ide,i,k) = lat_tmp(ide-i+1,jds:jde,k)
          enddo
          
          k = 6
          do i = ids,ide
            mesh%x  (ids:ide,ide-i+1,k) = x_tmp  (ide-i+1,ide:ids:-1,k)
            mesh%y  (ids:ide,ide-i+1,k) = y_tmp  (ide-i+1,ide:ids:-1,k)
            mesh%z  (ids:ide,ide-i+1,k) = z_tmp  (ide-i+1,ide:ids:-1,k)
            mesh%lon(ids:ide,ide-i+1,k) = lon_tmp(ide-i+1,ide:ids:-1,k)
            mesh%lat(ids:ide,ide-i+1,k) = lat_tmp(ide-i+1,ide:ids:-1,k)
          enddo
        elseif(trim(adjustl(gen_method))=='Rancic')then
          call inicuco
          do k = 1,Nf
            !$OMP PARALLEL DO PRIVATE(xm,xc,i)
             do j = 1, ny
                do i = 1, nx
                  xm(1) = mesh%xi (i,j,k)
                  xm(2) = mesh%eta(i,j,k)
                  call xmtoxc(xm,xc,k)
                  mesh%x(i,j,k) = xc(1)
                  mesh%y(i,j,k) = xc(2)
                  mesh%z(i,j,k) = xc(3)
                  call cart2sph(mesh%lon(i,j,k),mesh%lat(i,j,k),mesh%x(i,j,k),mesh%y(i,j,k),mesh%z(i,j,k))
                end do! i loop 
             end do! j loop
             !$OMP END PARALLEL DO
          end do! k loop
        endif
        
        x_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%x
        y_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%y
        z_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%z
        lon_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lon
        lat_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lat
        
        call CubedSphereFillHalo(x_tmp)
        call CubedSphereFillHalo(y_tmp)
        call CubedSphereFillHalo(z_tmp)
        
        ! Calculate matrices
        dh = dx * D2R / 2.
        do k = ifs,ife
          !$OMP PARALLEL DO PRIVATE(i,lon,lat,dlondx,dlondy,dlondz,dlatdx,dlatdy,dlatdz,sph2cart_matrix,fc)
          do j = jds, jde
            do i = ids, ide
              lon = mesh%lon(i,j,k)
              lat = mesh%lat(i,j,k)
              
              dlondx = -sin(lon); dlatdx = -cos(lon)*sin(lat)
              dlondy =  cos(lon); dlatdy = -sin(lon)*sin(lat)
              dlondz = 0.       ; dlatdz =  cos(lat)
              
              sph2cart_matrix(1,1) = dlondx
              sph2cart_matrix(1,2) = dlondy
              sph2cart_matrix(1,3) = dlondz
              sph2cart_matrix(2,1) = dlatdx
              sph2cart_matrix(2,2) = dlatdy
              sph2cart_matrix(2,3) = dlatdz
              
              !dxdxi
              fc(1:stencil_width) = x_tmp(i-xhalo:i+xhalo,j,k)
              mesh%jab(1,1,i,j,k) = center_difference(fc,dh)
              !dydxi
              fc(1:stencil_width) = y_tmp(i-xhalo:i+xhalo,j,k)
              mesh%jab(2,1,i,j,k) = center_difference(fc,dh)
              !dzdxi
              fc(1:stencil_width) = z_tmp(i-xhalo:i+xhalo,j,k)
              mesh%jab(3,1,i,j,k) = center_difference(fc,dh)
              !dxdeta
              fc(1:stencil_width) = x_tmp(i,j-xhalo:j+xhalo,k)
              mesh%jab(1,2,i,j,k) = center_difference(fc,dh)
              !dydeta
              fc(1:stencil_width) = y_tmp(i,j-xhalo:j+xhalo,k)
              mesh%jab(2,2,i,j,k) = center_difference(fc,dh)
              !dzdeta
              fc(1:stencil_width) = z_tmp(i,j-xhalo:j+xhalo,k)
              mesh%jab(3,2,i,j,k) = center_difference(fc,dh)
            enddo
          enddo
          !$OMP END PARALLEL DO
        enddo
        
      end subroutine generateConformalCubedSphere
      
    END MODULE mesh_mod


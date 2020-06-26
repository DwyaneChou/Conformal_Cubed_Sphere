
    MODULE mesh_mod
      use constants_mod
      use parameters_mod
      
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
        
      end type mesh_info
      
      type(mesh_info), target :: mesh
      
      contains
      
      subroutine initMesh
        integer :: iPV, jPV, iCell, jCell, iPatch, iVar, iDOF, jDOF
        integer :: iPVs, iPVe, jPVs, jPVe
        
        ! Allocate arrays in structures
        allocate( mesh%xi (ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%eta(ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%x  (ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%y  (ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%z  (ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lon(ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lat(ids:ide, jds:jde, ifs:ife) )
        
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
      
    END MODULE mesh_mod


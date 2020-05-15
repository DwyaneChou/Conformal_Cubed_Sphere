
    MODULE mesh_mod
      use constants_mod
      use parameters_mod
      
      implicit none
      
      ! coordinate
      type mesh_info
        real, dimension(:,:,:    ), allocatable :: x        ! central angle on x direction for cells on each patch, unit: radian
        real, dimension(:,:,:    ), allocatable :: y        ! central angle on y direction for cells on each patch, unit: radian
        real, dimension(:,:,:    ), allocatable :: lon      ! longitude on cells
        real, dimension(:,:,:    ), allocatable :: lat      ! latitude on cells
        real, dimension(:,:,:    ), allocatable :: sqrtG    ! jacobian of Transformation, sqrt(G)
        real, dimension(:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
        real, dimension(:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
        real, dimension(:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
        real, dimension(:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform
        
        real, dimension(:,:,:    ), allocatable :: f       ! Coriolis parameter
        
        real, dimension(:,:,:    ), allocatable :: sinlon  ! sin(longitude)
        real, dimension(:,:,:    ), allocatable :: coslon  ! cos(longitude)
        real, dimension(:,:,:    ), allocatable :: sinlat  ! sin(latitude)
        real, dimension(:,:,:    ), allocatable :: coslat  ! cos(latitude)
        
      end type mesh_info
      
      type(mesh_info), target :: mesh
      
      contains
      
      subroutine initMesh
        integer :: iPV, jPV, iCell, jCell, iPatch, iVar, iDOF, jDOF
        integer :: iPVs, iPVe, jPVs, jPVe
        
        ! Allocate arrays in structures
        allocate( mesh%x        (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%y        (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lon      (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lat      (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%sqrtG    (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixG  (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixIG (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixA  (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixIA (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%f        (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%sinlon   (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%coslon   (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%sinlat   (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%coslat   (      ids:ide, ids:ide, ifs:ife) )
        
        ! Calculate mesh infomation on VIA
        do iPatch = ifs, ife
          do jCell = jds, jde
            do iCell = ids, ide
              mesh%x(iCell, jCell, iPatch) = (iCell - 0.5) * dx + x_min * R2M
              mesh%y(iCell, jCell, iPatch) = (jCell - 0.5) * dy + y_min * R2M
              
              !call pointProjPlane2Sphere(mesh%lon(iCell, jCell, iPatch), mesh%lat(iCell, jCell, iPatch), &
              !                           mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch), iPatch)
              !
              !mesh%sinlon(iCell, jCell, iPatch) = sin(mesh%lon(iCell, jCell, iPatch))
              !mesh%coslon(iCell, jCell, iPatch) = cos(mesh%lon(iCell, jCell, iPatch))
              !
              !mesh%sinlat(iCell, jCell, iPatch) = sin(mesh%lat(iCell, jCell, iPatch))
              !mesh%coslat(iCell, jCell, iPatch) = cos(mesh%lat(iCell, jCell, iPatch))
              !
              !call calc_matrixG (mesh%matrixG (:, :, iCell, jCell, iPatch), mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch))
              !call calc_matrixIG(mesh%matrixIG(:, :, iCell, jCell, iPatch), mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch))
              !call calc_matrixA (mesh%matrixA (:, :, iCell, jCell, iPatch), mesh%lon(iCell, jCell, iPatch), mesh%lat(iCell, jCell, iPatch), iPatch)
              !call calc_matrixIA(mesh%matrixIA(:, :, iCell, jCell, iPatch), mesh%lon(iCell, jCell, iPatch), mesh%lat(iCell, jCell, iPatch), iPatch)
              !call calc_Jacobian(mesh%sqrtG   (      iCell, jCell, iPatch), mesh%x  (iCell, jCell, iPatch), mesh%y  (iCell, jCell, iPatch))
              !
              !mesh%f(iCell, jCell, iPatch) = 2. * Omega * mesh%sinlat(iCell, jCell, iPatch)
            end do
          end do
        end do
      end subroutine initMesh
      
    END MODULE mesh_mod


    ! Conformal Cubed Sphere mesh generator
    program CCS
    use parameters_mod
    use mesh_mod
    use math_mod
    use jimcc_m
    use output_mod
    implicit none
    
    real lon,lat
    real xx,yy,xc(3)
    
    integer i,j,iPatch
    
    call initParameters
    
    call initMesh
    
    call inrot()
    
    do iPatch = 1,6
      do j = jds, jde
        do i = ids, ide
          xx = mesh%x(i,j,iPatch)
          yy = mesh%y(i,j,iPatch)
          
          call mtoc(xx, yy, iPatch, xc)
          call cart2sph(lon,lat,xc(1),xc(2),xc(3))
          
          mesh%lon(i,j,iPatch) = lon
          mesh%lat(i,j,iPatch) = lat
        enddo
      enddo
    enddo
          
    call history_init
    
    end program CCS

MODULE spatial_operators_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  use tend_mod
  use projection_mod
  
  implicit none
  
    contains
    
    subroutine spatial_operator(stat,tend)
      type(stat_field), target, intent(in   ) :: stat
      type(tend_field), target, intent(inout) :: tend
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: E
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: K
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: flux_x
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: flux_y
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: vorticity
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: phit
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dEdy
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dvdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dudy
      
      integer i,j,iPatch
      integer ip1,jp1,ip2,jp2
      integer im1,jm1,im2,jm2
      
      phit  = stat%phi + mesh%phis
      
      K = 0.5 * ( stat%u * stat%uc + stat%v * stat%vc )
      E = phit + K
      
      flux_x = stat%phiG * stat%uc
      flux_y = stat%phiG * stat%vc
      
      do iPatch = ifs, ife
        do j = jts, jte
          do i = its, ite
            ip1 = i + 1
            im1 = i - 1
            jp1 = j + 1
            jm1 = j - 1
            ip2 = i + 2
            im2 = i - 2
            jp2 = j + 2
            jm2 = j - 2
            
            dfluxdx  (i,j,iPatch) = ( flux_x(im2,j,iPatch) - 8. * flux_x(im1,j,iPatch) + 8. * flux_x(ip1,j,iPatch) - flux_x(ip2,j,iPatch) ) / ( 12. * dx )
            dEdx     (i,j,iPatch) = ( E     (im2,j,iPatch) - 8. * E     (im1,j,iPatch) + 8. * E     (ip1,j,iPatch) - E     (ip2,j,iPatch) ) / ( 12. * dx )
            dvdx     (i,j,iPatch) = ( stat%v(im2,j,iPatch) - 8. * stat%v(im1,j,iPatch) + 8. * stat%v(ip1,j,iPatch) - stat%v(ip2,j,iPatch) ) / ( 12. * dx )
            dfluxdy  (i,j,iPatch) = ( flux_y(i,jm2,iPatch) - 8. * flux_y(i,jm1,iPatch) + 8. * flux_y(i,jp1,iPatch) - flux_y(i,jp2,iPatch) ) / ( 12. * dy )
            dEdy     (i,j,iPatch) = ( E     (i,jm2,iPatch) - 8. * E     (i,jm1,iPatch) + 8. * E     (i,jp1,iPatch) - E     (i,jp2,iPatch) ) / ( 12. * dy )
            dudy     (i,j,iPatch) = ( stat%u(i,jm2,iPatch) - 8. * stat%u(i,jm1,iPatch) + 8. * stat%u(i,jp1,iPatch) - stat%u(i,jp2,iPatch) ) / ( 12. * dy )
            
            vorticity(i,j,iPatch) = dvdx(i,j,iPatch) - dudy(i,j,iPatch) + mesh%sqrtG(i,j,iPatch) * mesh%f(i,j,iPatch)
          enddo
        enddo
      enddo
      
      tend%phiG = - ( dfluxdx + dfluxdy )
      tend%u    = - dEdx + vorticity * stat%vc
      tend%v    = - dEdy - vorticity * stat%uc
      
    end subroutine spatial_operator
    
    ! convert vector from patch to sphere
    subroutine convert_wind_P2SP(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i)
        do j = jts, jte
          do i = its, ite
            call covProjPlane2Sphere(stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), stat%u(i,j,iPatch), stat%v(i,j,iPatch), mesh%matrixA(:,:,i,j,iPatch), mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
    end subroutine convert_wind_P2SP
    
    subroutine convert_wind_SP2P(stat)
      type(stat_field), intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i)
        do j = jcs, jce
          do i = ics, ice
            call covProjSphere2Plane(stat%u(i,j,iPatch), stat%v(i,j,iPatch), stat%zonal_wind(i,j,iPatch), stat%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
    end subroutine convert_wind_SP2P
    
    subroutine convert_wind_cov2contrav(stat)
      type(stat_field), target, intent(inout) :: stat
      
      integer i,j,iPatch
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i)
        do j = jcs, jce
          do i = ics, ice
            call cov2contrav(stat%uc(i,j,iPatch),stat%vc(i,j,iPatch),stat%u(i,j,iPatch),stat%v(i,j,iPatch),mesh%matrixIG(:,:,i,j,iPatch))
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
    end subroutine convert_wind_cov2contrav
    
END MODULE spatial_operators_mod


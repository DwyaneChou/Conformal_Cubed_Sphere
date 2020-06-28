module test_case_mod
  use parameters_mod
  use mesh_mod
  use ghost_mod
  use stat_mod
  use projection_mod
  use spatial_operators_mod
  use diag_mod
  implicit none
  
  contains
    
  subroutine initTestCase
    integer i,j,iPatch

    real, dimension(:,:,:), pointer :: phiG
    real, dimension(:,:,:), pointer :: u
    real, dimension(:,:,:), pointer :: v
    
    phiG => stat(0)%phiG
    u    => stat(0)%u   
    v    => stat(0)%v   
    
    if(case_num == 1)then
      print*,''
      print*,' gaussian hill is selected'
      call gaussian_hill(stat(0))
    else
      print*, 'Unknown test case, choose from case_num=1'
      stop
    endif
    
    ! initialize phiG
    stat(0)%phiG = stat(0)%phi * mesh%sqrtG
    
    do iPatch = ifs, ife
      do j = jcs, jce
        do i = ics, ice
          call covProjSphere2Plane    (stat(0)%u (i,j,iPatch), stat(0)%v (i,j,iPatch), stat(0)%zonal_wind(i,j,iPatch), stat(0)%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch), mesh%matrixG(:,:,i,j,iPatch))
          call contravProjSphere2Plane(stat(0)%uC(i,j,iPatch), stat(0)%vC(i,j,iPatch), stat(0)%zonal_wind(i,j,iPatch), stat(0)%meridional_wind(i,j,iPatch), mesh%matrixIA(:,:,i,j,iPatch))
        enddo
      enddo
    enddo
      
    print*,''
    print*,'max/min value of phi : ',maxval(stat(0)%phi),minval(stat(0)%phi)
    print*,'max/min value of u   : ',maxval(stat(0)%u),minval(stat(0)%u)
    print*,'max/min value of v   : ',maxval(stat(0)%v),minval(stat(0)%v)
    print*,'max/min value of uc  : ',maxval(stat(0)%uc),minval(stat(0)%uc)
    print*,'max/min value of vc  : ',maxval(stat(0)%vc),minval(stat(0)%vc)
    print*,'max/min value of us  : ',maxval(stat(0)%zonal_wind),minval(stat(0)%zonal_wind)
    print*,'max/min value of vs  : ',maxval(stat(0)%meridional_wind),minval(stat(0)%meridional_wind)
  end subroutine initTestCase

  ! Global steady flow
  subroutine gaussian_hill(stat)
    type(stat_field), intent(inout) :: stat
    real,dimension(ics:ice,jcs:jce,ifs:ife) :: u,v,phi                 ! working array
    
    real    :: u0
    real    :: alpha0
    real    :: b0
    real    :: lambda_c
    real    :: theta_c
    real    :: gh0
    
    real    :: X,Y,Z
    real    :: Xc,Yc,Zc
    
    integer :: i,j,iPatch
    
    gh0      = 100. * gravity
    alpha0   = 0.
    b0       = 0.5 * 10.e-14
    lambda_c = 0. !3. * pi / 2.
    theta_c  = 0.
    u0       = 2. * pi * radius / (12. * 86400.)
    Xc       = radius * cos(lambda_c) * cos(theta_c)
    Yc       = radius * sin(lambda_c) * cos(theta_c)
    Zc       = radius * sin(theta_c)
    
    do iPatch = ifs, ife
      do j = jcs, jce
        do i = ics, ice
          X = radius * cos(mesh%lon(i,j,iPatch)) * cos(mesh%lat(i,j,iPatch))
          Y = radius * sin(mesh%lon(i,j,iPatch)) * cos(mesh%lat(i,j,iPatch))
          Z = radius * sin(mesh%lat(i,j,iPatch))
          
          phi(i,j,iPatch) = gh0 * exp(-b0 * ((X-Xc)**2 + (Y-Yc)**2 + (Z-Zc)**2))
          u  (i,j,iPatch) = u0 * (cos(mesh%lat(i,j,iPatch)) * cos(alpha0) + sin(alpha0) * cos(mesh%lon(i,j,iPatch)) * sin(mesh%lat(i,j,iPatch)))
          v  (i,j,iPatch) = -u0 * sin(alpha0) * sin(mesh%lon(i,j,iPatch))
        enddo
      enddo
    enddo
    
    stat%zonal_wind      = u
    stat%meridional_wind = v
    stat%phi             = phi
    
    mesh%phis = 0.
    
  end subroutine gaussian_hill
  
end module test_case_mod
    
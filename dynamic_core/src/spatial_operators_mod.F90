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
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: flux_x
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: flux_y
      
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdx
      real, dimension(ics:ice,jcs:jce,ifs:ife) :: dfluxdy
      
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext_n
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext_p
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext_n
      
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_x_ext
      real, dimension(ids:ide,jds:jde,ifs:ife) :: flux_y_ext
      
      real maxeigen_x
      real maxeigen_y
      
      integer i,j,iPatch
      integer ip1,jp1,ip2,jp2,ip3,jp3
      integer im1,jm1,im2,jm2,im3,jm3
      integer ic,jc
      
      flux_x = stat%phiG * stat%uc
      flux_y = stat%phiG * stat%vc
      
      call weno(flux_x_ext_p, flux_x   , 1, 1)
      call weno(flux_x_ext_n, flux_x   ,-1, 1)
      call weno(flux_y_ext_p, flux_y   , 1, 2)
      call weno(flux_y_ext_n, flux_y   ,-1, 2)
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,ip1,jc,maxeigen_x)
        do j = jts, jte
          do i = its-1, ite
            ip1 = 2 * i + 1
            jc  = 2 * j
            
            maxeigen_x = 0.5 * ( stat%uc(i,j,iPatch) + stat%uc(i+1,j,iPatch) )
            
            flux_x_ext(ip1,jc,iPatch) = 0.5 * ( flux_x_ext_p(ip1,jc,iPatch) + flux_x_ext_n(ip1,jc,iPatch) - sign(1.,maxeigen_x) * ( flux_x_ext_n(ip1,jc,iPatch) - flux_x_ext_p(ip1,jc,iPatch) ) )
          enddo
        enddo
        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(i,ic,jp1,maxeigen_y)
        do j = jts-1, jte
          do i = its, ite
            ic  = 2 * i
            jp1 = 2 * j + 1
            
            maxeigen_y = 0.5 * ( stat%vc(i,j,iPatch) + stat%vc(i,j+1,iPatch) )
            
            flux_y_ext(ic,jp1,iPatch) = 0.5 * ( flux_y_ext_p(ic,jp1,iPatch) + flux_y_ext_n(ic,jp1,iPatch) - sign(1.,maxeigen_y) * ( flux_y_ext_n(ic,jp1,iPatch) - flux_y_ext_p(ic,jp1,iPatch) ) )
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      do iPatch = ifs, ife
        !$OMP PARALLEL DO PRIVATE(i,ic,jc,ip1,im1,jp1,jm1)
        do j = jts, jte
          do i = its, ite
            ic  = 2 * i
            jc  = 2 * j
            ip1 = 2 * i + 1
            im1 = 2 * i - 1
            jp1 = 2 * j + 1
            jm1 = 2 * j - 1
            
            dfluxdx  (i,j,iPatch) = ( flux_x_ext(ip1,jc,iPatch) - flux_x_ext(im1,jc,iPatch) ) / dx
            dfluxdy  (i,j,iPatch) = ( flux_y_ext(ic,jp1,iPatch) - flux_y_ext(ic,jm1,iPatch) ) / dy
          enddo
        enddo
        !$OMP END PARALLEL DO
      enddo
      
      tend%phiG = - ( dfluxdx + dfluxdy )
      tend%u    = 0.
      tend%v    = 0.
      
    end subroutine spatial_operator
    
    subroutine weno(field_ext,field,flux_dir,axis_dir)
      real   , dimension(ids:ide,jds:jde,ifs:ife), intent(out) :: field_ext
      real   , dimension(ics:ice,jcs:jce,ifs:ife), intent(in ) :: field
      integer,                                     intent(in ) :: flux_dir ! 1 for positive flux, -1 for negative flux
      integer,                                     intent(in ) :: axis_dir ! 1 for x axis, 2 for y axis
      
      real Qx(5)
      real Qy(5)
      
      integer i,j,iPatch
      integer ic,jc
      integer iext
      integer jext
      
      if(axis_dir==1)then
        if(flux_dir>0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(i,Qx,jc,iext)
            do j = jts, jte
              do i = its-1, ite
                Qx(1) = field(i-2,j,iPatch)
                Qx(2) = field(i-1,j,iPatch)
                Qx(3) = field(i  ,j,iPatch)
                Qx(4) = field(i+1,j,iPatch)
                Qx(5) = field(i+2,j,iPatch)
                
                jc   = 2 * j
                iext = 2 * i + 1
                call WENO_limiter(field_ext(iext,jc,iPatch),Qx,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        elseif(flux_dir<0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(i,Qx,jc,iext)
            do j = jts, jte
              do i = its, ite+1
                Qx(1) = field(i-2,j,iPatch)
                Qx(2) = field(i-1,j,iPatch)
                Qx(3) = field(i  ,j,iPatch)
                Qx(4) = field(i+1,j,iPatch)
                Qx(5) = field(i+2,j,iPatch)
                
                jc   = 2 * j
                iext = 2 * i - 1
                call WENO_limiter(field_ext(iext,jc,iPatch),Qx,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        endif
      endif
        
      if(axis_dir==2)then
        if(flux_dir>0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(j,Qy,ic,jext)
            do i = its, ite
              do j = jts-1, jte
                Qy(1) = field(i,j-2,iPatch)
                Qy(2) = field(i,j-1,iPatch)
                Qy(3) = field(i,j  ,iPatch)
                Qy(4) = field(i,j+1,iPatch)
                Qy(5) = field(i,j+2,iPatch)
                
                ic   = 2 * i
                jext = 2 * j + 1
                call WENO_limiter(field_ext(ic,jext,iPatch),Qy,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        elseif(flux_dir<0)then
          do iPatch = ifs, ife
            !$OMP PARALLEL DO PRIVATE(j,Qy,ic,jext)
            do i = its, ite
              do j = jts, jte+1
                Qy(1) = field(i,j-2,iPatch)
                Qy(2) = field(i,j-1,iPatch)
                Qy(3) = field(i,j  ,iPatch)
                Qy(4) = field(i,j+1,iPatch)
                Qy(5) = field(i,j+2,iPatch)
                
                ic   = 2 * i
                jext = 2 * j - 1
                call WENO_limiter(field_ext(ic,jext,iPatch),Qy,flux_dir)
              enddo
            enddo
            !$OMP END PARALLEL DO
          enddo
        endif
      endif
        
    end subroutine weno
    
    ! 1D WENO slope limiter, according to Sun,2015
    ! "A Slope Constrained 4th OrderMulti-Moment Finite Volume Method with WENO Limiter"
    ! and Jiang and Shu, 1996
    subroutine WENO_limiter(Qrec,Q,dir)
      real              , intent(out) :: Qrec
      real, dimension(5), intent(in ) :: Q
      integer           , intent(in ) :: dir
      
      integer, parameter :: nStencil = 3
      real   , parameter :: weno_coef(3)  = [0.1, 0.6, 0.3]
      real   , parameter :: eps           = 1.E-2
      
      real Qim(nStencil-1)
      real Qip(nStencil-1)
      real Qi
      
      real, dimension(nStencil) :: stencil
      real, dimension(nStencil) :: coefA
      real, dimension(nStencil) :: coefB
      real, dimension(nStencil) :: alpha
      real, dimension(nStencil) :: beta
      real, dimension(nStencil) :: omega
      
      real tau40
      real tau41
      real tau5
      
      integer iStencil
      
      Qim(2) = Q(1)
      Qim(1) = Q(2)
      Qi     = Q(3)
      Qip(1) = Q(4)
      Qip(2) = Q(5)
      
      if(dir>0)then
        stencil (1) =  Qim(2)/3. - 7./6. * Qim(1) + 11./6. * Qi     
        stencil (2) = -Qim(1)/6. + 5./6. * Qi     +  1./3. * Qip(1) 
        stencil (3) =  Qi    /3. + 5./6. * Qip(1) -  1./6. * Qip(2)
        
        coefA(1) = Qim(2) - 2. * Qim(1) + Qi
        coefA(2) = Qim(1) - 2. * Qi     + Qip(1)
        coefA(3) = Qi     - 2. * Qip(1) + Qip(2)
        
        coefB(1) =      Qim(2) - 4. * Qim(1) + 3. * Qi
        coefB(2) =      Qim(1) -      Qip(1)
        coefB(3) = 3. * Qi     - 4. * Qip(1) +      Qip(2)
      elseif(dir<0)then
        stencil (1) =  Qip(2)/3. - 7./6. * Qip(1) + 11./6. * Qi     
        stencil (2) = -Qip(1)/6. + 5./6. * Qi     +  1./3. * Qim(1) 
        stencil (3) =  Qi    /3. + 5./6. * Qim(1) -  1./6. * Qim(2)
        
        coefA(1) = Qip(2) - 2. * Qip(1) + Qi
        coefA(2) = Qip(1) - 2. * Qi     + Qim(1)
        coefA(3) = Qi     - 2. * Qim(1) + Qim(2)
        
        coefB(1) =      Qip(2) - 4. * Qip(1) + 3. * Qi
        coefB(2) =      Qip(1) -      Qim(1)
        coefB(3) = 3. * Qi     - 4. * Qim(1) +      Qim(2)
      endif
      
      beta = coefA**2 * 13. / 12. + coefB**2 * 0.25
      
      ! WENO-Z
      tau40 = abs( beta(1) - beta(2) )
      tau41 = abs( beta(2) - beta(3) )
      tau5  = abs( beta(3) - beta(1) )
      
      ! xi
      if( tau40<=minval(beta) .and. tau41>minval(beta) )then
        omega = [1./3.,2./3.,0.]
      elseif( tau40>minval(beta) .and. tau41<minval(beta) )then
        omega = [0.,2./3.,1./3.]
      else
        alpha = weno_coef * ( 1. + tau5 / ( eps + beta ) )
        omega = alpha / sum(alpha)
      endif
      
      !omega = weno_coef
      
      Qrec = dot_product( stencil, omega )
      
    end subroutine WENO_limiter
    
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


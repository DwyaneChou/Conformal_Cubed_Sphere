    ! Conformal Cubed Sphere mesh generator
    program CCS
      use parameters_mod
      use mesh_mod
      use output_mod
      use cuco
      use math_mod
      implicit none
      
      real,dimension(2):: xm
      real,dimension(3):: xc

      integer i,j,n
      
      call initParameters
      
      call initMesh

      do n = 1,Nf
        !$OMP PARALLEL DO PRIVATE(xm,xc)
         do j = 1, ny
            do i = 1, nx
              xm(1) = mesh%xi (i,j,n)
              xm(2) = mesh%eta(i,j,n)
              call xmtoxc(xm,xc,n)
              mesh%x(i,j,n) = xc(1)
              mesh%y(i,j,n) = xc(2)
              mesh%z(i,j,n) = xc(3)
              call cart2sph(mesh%lon(i,j,n),mesh%lat(i,j,n),mesh%x(i,j,n),mesh%y(i,j,n),mesh%z(i,j,n))
            end do! i loop 
         end do! j loop
         !$OMP END PARALLEL DO
      end do! n loop 
      
      call history_init
    
    end program CCS

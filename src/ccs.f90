    ! Conformal Cubed Sphere mesh generator
    program CCS
      use parameters_mod
      use mesh_mod
      use setxyz_m
      use xyzinfo_m
      use newmpar_m, only : il,jl,kl,npanels,ifull,iquad
      use ind_m
      use output_mod
      implicit none
      
      integer             :: diag
      integer             :: id,jd
      real                :: rlong0, rlat0
      !real                :: schmidt
      real   , parameter  :: schm13 = 0.1
      !integer             :: ntang
      real                :: rearth
      
      integer iq
      integer i,j,n
      
      call initParameters
      
      call initMesh
      
      il      = nx
      jl      = il * 6
      kl      = 1
      npanels = jl/il - 1
      ifull   = il * jl
      iquad   = 1+il*((8*npanels)/(npanels+4))
      diag    = 3
      id      = 1
      jd      = 1
      rlong0  = 0.
      rlat0   = 0.
      !schmidt = 1.
      !ntang   = 0
      rearth  = radius
      
      call setxyz( il, jl, kl, npanels, ifull, iquad, diag, id, jd, &
                   rlong0, rlat0, schmidt, schm13, ntang, rearth )
      call setaxu(ifull)
      
      do n = 0, npanels 
         do j = 1, ny
            do i = 1, nx
              iq = ind(i,j,n)
              mesh%x  (i,j,n+1) = x    (iq)
              mesh%y  (i,j,n+1) = y    (iq)
              mesh%lon(i,j,n+1) = rlong(iq)
              mesh%lat(i,j,n+1) = rlat (iq)
            end do! i loop 
         end do! j loop 
      end do! n loop 
      
      call history_init
    
    end program CCS

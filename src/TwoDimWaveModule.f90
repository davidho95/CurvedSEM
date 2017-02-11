module TwoDimWaveModule

  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: PI = 3.141592653589793

  integer, dimension(2), parameter :: NUM_SPEC_EL = (/1, 1/)
  integer, parameter :: NUM_GLL = 3
  integer, dimension(2), parameter :: NUM_GLOBAL_POINTS = NUM_SPEC_EL * (NUM_GLL - 1) + (/1, 1/)
  contains

    subroutine Setup_mesh_2d()
      real(dp), dimension(0:NUM_SPEC_EL(1), 0:NUM_SPEC_EL(2), 2) :: anchor_points
      real(dp), dimension(NUM_GLL) :: gll_points, gll_weights
      real(dp), dimension(NUM_GLL, NUM_GLL) :: hprime

      integer, parameter :: TOTAL_NUM_SPEC = NUM_SPEC_EL(1) * NUM_SPEC_EL(2)
      real(dp), dimension(NUM_GLL, NUM_GLL, TOTAL_NUM_SPEC, 2) :: global_points

      integer i_gll, j_gll, i_spec, i_spec_x, i_spec_y

      ! get the GLL points and weights
      call zwgljd(gll_points,gll_weights,NUM_GLL,0.0_dp,0.0_dp)
      if(mod(NUM_GLL,2) /= 0) gll_points((NUM_GLL-1)/2+1) = 0.0_dp
      
      ! get the derivatives of the Lagrange polynomials at 
      ! the GLL points; recall that  hprime(i,j)=h'_{j}(gll_points_{i}) 
      do j_gll=1,NUM_GLL
         do i_gll=1,NUM_GLL
            hprime(i_gll,j_gll) = lagrange_deriv_GLL(j_gll-1,i_gll-1,gll_points,NUM_GLL)
         end do
      end do

      ! Setup anchor points
      do i_spec_x = 0, NUM_SPEC_EL(1)
        do i_spec_y = 0, NUM_SPEC_EL(2)
          anchor_points(i_spec_x, i_spec_y, 1) = 2 * PI &
           * dble(i_spec_x) / dble(NUM_SPEC_EL(1))
          anchor_points(i_spec_x, i_spec_y, 2) = 2 * PI &
           * dble(i_spec_y) / dble(NUM_SPEC_EL(2)) 
        enddo
      enddo

      ! Setup global points
      i_spec = 0
      do i_spec_x = 1, NUM_SPEC_EL(1)
       do i_spec_y = 1, NUM_SPEC_EL(2)
          i_spec = i_spec + 1
          do i_gll = 1, NUM_GLL
             do j_gll = 1, NUM_GLL
                global_points(i_gll,j_gll,i_spec, 1) =&
                  anchor_points(i_spec_x-1, i_spec_y-1, 1) +&
                  (gll_points(i_gll)+1.0_dp)*0.5_dp*&
                  (anchor_points(i_spec_x, i_spec_y, 1) -&
                  anchor_points(i_spec_x - 1, i_spec_y - 1, 1))
                global_points(i_gll,j_gll,i_spec, 2) =&
                  anchor_points(i_spec_x-1, i_spec_y-1, 2) +&
                  (gll_points(i_gll)+1.0_dp)*0.5_dp*&
                  (anchor_points(i_spec_x, i_spec_y, 2) -&
                  anchor_points(i_spec_x - 1, i_spec_y - 1, 2))
             enddo
          enddo
       enddo
    enddo

    end subroutine Setup_mesh_2d

end module TwoDimWaveModule
# if defined BIO_HADV_WENO5 || \
     defined TS_HADV_WENO5  || defined TS_VADV_WENO5 || \
     defined UV_HADV_WENO5  || defined UV_VADV_WENO5 || \
     defined W_HADV_WENO5   || defined W_VADV_WENO5  || \
     defined BEDLOAD_WENO5  || \
     defined BIO_HADV_WENO7 || \
     defined TS_HADV_WENO7  || defined TS_VADV_WENO7 || \
     defined UV_HADV_C8     || defined UV_HADV_UP7   || \
     defined UV_HADV_WENO7  || \
     defined BIO_HADV_WENO9 || \
     defined TS_HADV_WENO9  || defined TS_VADV_WENO9 || \
     defined UV_HADV_C10     || defined UV_HADV_UP9   || \
     defined UV_HADV_WENO9

      INTERFACE

      ! flux3_weno
      function flux3_weno( q_im2, q_im1, q_i, q_ip1, ua)
!$acc routine(flux3_weno) seq
      REAL    :: flux3_weno
      REAL    :: q_im2, q_im1, q_i, q_ip1, ua
      end function flux3_weno 

      ! flux5_weno
      function flux5_weno(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
!$acc routine(flux5_weno) seq
      implicit none
      REAL :: flux5_weno
      REAL, INTENT(IN) :: q_im3, q_im2, q_im1,
     &                    q_i,   q_ip1, q_ip2, ua
      end function flux5_weno

      ! flux7_weno
      function flux7_weno(q_im4, q_im3, q_im2, q_im1,
     &                    q_i,   q_ip1, q_ip2, q_ip3, ua)
!$acc routine(flux7_weno) seq
      implicit none
      REAL :: flux7_weno
      REAL, INTENT(IN) :: q_im4, q_im3, q_im2, q_im1,
     &                    q_i,   q_ip1, q_ip2, q_ip3, ua
      end function flux7_weno

      ! flux9_weno
      function flux9_weno(q_im5, q_im4, q_im3, q_im2, q_im1,
     &                    q_i,   q_ip1, q_ip2, q_ip3, q_ip4, ua)
!$acc routine(flux9_weno) seq
      implicit none
      REAL :: flux9_weno
      REAL, INTENT(IN) :: q_im5, q_im4, q_im3, q_im2, q_im1,
     &                    q_i,   q_ip1, q_ip2, q_ip3, q_ip4, ua
      end function flux9_weno

      END INTERFACE
#endif

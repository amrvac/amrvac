!> User can use subroutine get_atm_para to generate 1D solar stmosphere.
!> User should provide heights (h), number density at the bottom,
!> number of points (nh), and the product of avaerage ion mass and 
!> gravity (mg) at each point. 
!> User can select temperature profile.
!> This subroutine will return density and pressure at each point.
!> please use cgs unit.
module mod_solar_atmosphere
  use mod_global_parameters
  use mod_constants
  implicit none

  integer :: n_valc,n_hong
  double precision :: h_valc(1:52),T_valc(1:52)
  double precision :: h_hong(1:300),T_hong(1:300)

  data n_valc / 52 /

  data    h_valc /    -75,    -50,     25,      0,     50, & 
                      100,    150,    250,    350,    450, &
                      515,    555,    605,    655,    705, &
                      755,    855,    905,    980,   1065, &
                     1180,   1280,   1380,   1515,   1605, &
                     1785,   1925,   1990,   2016,   2050, &
                     2070,   2080,   2090,   2104,   2107, &
                     2109,   2113,   2115,   2120,   2129, &
                     2160,   2200,   2230,   2255,   2263, &
                     2267,   2271,   2274,   2280,   2290, &
                     2298,   2543 /

  data    T_valc /   8320,   7610,   6910,   6420,   5840, &  
                     5455,   5180,   4780,   4465,   4220, &
                     4170,   4230,   4420,   4730,   5030, &
                     5280,   5650,   5755,   5925,   6040, &
                     6150,   6220,   6280,   6370,   6440, &
                     6630,   6940,   7160,   7360,   7660, &
                     7940,   8180,   8440,   9500,  10700, &
                    12300,  18500,  21000,  22500,  23000, &
                    23500,  24000,  24200,  24500,  25500, &
                    28000,  32000,  37000,  50000,  89100, &
                   141000, 447000 /


  data n_hong / 300 /

  data h_hong / -3.7616e+07, -3.7384e+07, -3.7045e+07, -3.6561e+07, -3.5880e+07,& 
                -3.4950e+07, -3.3718e+07, -3.2153e+07, -3.0250e+07, -2.8050e+07,&
                -2.5627e+07, -2.3088e+07, -2.0556e+07, -1.8151e+07, -1.5972e+07,&
                -1.4004e+07, -1.2234e+07, -1.0625e+07, -9.1510e+06, -7.7773e+06,&
                -6.4599e+06, -5.1584e+06, -3.8455e+06, -2.4946e+06, -1.0577e+06,&
                 4.7767e+05,  2.1041e+06,  3.8237e+06,  5.6198e+06,  7.4817e+06,&
                 9.3930e+06,  1.1339e+07,  1.3304e+07,  1.5278e+07,  1.7250e+07,&
                 1.9213e+07,  2.1163e+07,  2.3099e+07,  2.5019e+07,  2.6922e+07,&
                 2.8809e+07,  3.0678e+07,  3.2532e+07,  3.4371e+07,  3.6196e+07,&
                 3.8009e+07,  3.9811e+07,  4.1605e+07,  4.3396e+07,  4.5186e+07,&
                 4.6977e+07,  4.8773e+07,  5.0576e+07,  5.2390e+07,  5.4218e+07,&
                 5.6065e+07,  5.7934e+07,  5.9833e+07,  6.1766e+07,  6.3736e+07,&
                 6.5744e+07,  6.7791e+07,  6.9880e+07,  7.2011e+07,  7.4185e+07,&
                 7.6403e+07,  7.8660e+07,  8.0954e+07,  8.3284e+07,  8.5647e+07,&
                 8.8038e+07,  9.0456e+07,  9.2896e+07,  9.5357e+07,  9.7838e+07,&
                 1.0033e+08,  1.0285e+08,  1.0537e+08,  1.0791e+08,  1.1046e+08,&
                 1.1302e+08,  1.1559e+08,  1.1818e+08,  1.2077e+08,  1.2337e+08,&
                 1.2599e+08,  1.2862e+08,  1.3126e+08,  1.3392e+08,  1.3659e+08,&
                 1.3927e+08,  1.4198e+08,  1.4470e+08,  1.4744e+08,  1.5019e+08,&
                 1.5295e+08,  1.5571e+08,  1.5847e+08,  1.6121e+08,  1.6389e+08,&
                 1.6648e+08,  1.6895e+08,  1.7122e+08,  1.7325e+08,  1.7499e+08,&
                 1.7642e+08,  1.7754e+08,  1.7838e+08,  1.7900e+08,  1.7944e+08,&
                 1.7974e+08,  1.7996e+08,  1.8010e+08,  1.8020e+08,  1.8027e+08,&
                 1.8032e+08,  1.8035e+08,  1.8037e+08,  1.8039e+08,  1.8041e+08,&
                 1.8042e+08,  1.8044e+08,  1.8045e+08,  1.8046e+08,  1.8047e+08,&
                 1.8049e+08,  1.8051e+08,  1.8053e+08,  1.8055e+08,  1.8059e+08,&
                 1.8063e+08,  1.8068e+08,  1.8074e+08,  1.8082e+08,  1.8092e+08,&
                 1.8103e+08,  1.8118e+08,  1.8135e+08,  1.8157e+08,  1.8182e+08,&
                 1.8214e+08,  1.8252e+08,  1.8299e+08,  1.8357e+08,  1.8428e+08,&
                 1.8515e+08,  1.8621e+08,  1.8749e+08,  1.8904e+08,  1.9089e+08,&
                 1.9307e+08,  1.9561e+08,  1.9852e+08,  2.0180e+08,  2.0544e+08,&
                 2.0941e+08,  2.1368e+08,  2.1819e+08,  2.2291e+08,  2.2781e+08,&
                 2.3284e+08,  2.3798e+08,  2.4320e+08,  2.4849e+08,  2.5383e+08,&
                 2.5922e+08,  2.6463e+08,  2.7008e+08,  2.7554e+08,  2.8102e+08,&
                 2.8651e+08,  2.9202e+08,  2.9753e+08,  3.0305e+08,  3.0858e+08,&
                 3.1412e+08,  3.1966e+08,  3.2520e+08,  3.3075e+08,  3.3630e+08,&
                 3.4186e+08,  3.4741e+08,  3.5297e+08,  3.5854e+08,  3.6410e+08,&
                 3.6967e+08,  3.7523e+08,  3.8080e+08,  3.8637e+08,  3.9195e+08,&
                 3.9752e+08,  4.0309e+08,  4.0867e+08,  4.1424e+08,  4.1982e+08,&
                 4.2540e+08,  4.3098e+08,  4.3656e+08,  4.4214e+08,  4.4772e+08,&
                 4.5330e+08,  4.5888e+08,  4.6446e+08,  4.7004e+08,  4.7563e+08,&
                 4.8121e+08,  4.8679e+08,  4.9238e+08,  4.9796e+08,  5.0355e+08,&
                 5.0913e+08,  5.1472e+08,  5.2030e+08,  5.2589e+08,  5.3148e+08,&
                 5.3706e+08,  5.4265e+08,  5.4824e+08,  5.5382e+08,  5.5941e+08,&
                 5.6500e+08,  5.7059e+08,  5.7617e+08,  5.8176e+08,  5.8735e+08,&
                 5.9294e+08,  5.9853e+08,  6.0412e+08,  6.0971e+08,  6.1530e+08,&
                 6.2089e+08,  6.2648e+08,  6.3207e+08,  6.3766e+08,  6.4325e+08,&
                 6.4884e+08,  6.5443e+08,  6.6002e+08,  6.6561e+08,  6.7120e+08,&
                 6.7679e+08,  6.8238e+08,  6.8798e+08,  6.9357e+08,  6.9916e+08,&
                 7.0475e+08,  7.1034e+08,  7.1593e+08,  7.2152e+08,  7.2712e+08,&
                 7.3271e+08,  7.3830e+08,  7.4389e+08,  7.4948e+08,  7.5508e+08,&
                 7.6067e+08,  7.6626e+08,  7.7185e+08,  7.7744e+08,  7.8304e+08,&
                 7.8863e+08,  7.9422e+08,  7.9981e+08,  8.0540e+08,  8.1100e+08,&
                 8.1659e+08,  8.2218e+08,  8.2777e+08,  8.3337e+08,  8.3896e+08,&
                 8.4455e+08,  8.5014e+08,  8.5574e+08,  8.6133e+08,  8.6692e+08,&
                 8.7252e+08,  8.7811e+08,  8.8370e+08,  8.8929e+08,  8.9489e+08,&
                 9.0048e+08,  9.0607e+08,  9.1166e+08,  9.1726e+08,  9.2285e+08,&
                 9.2844e+08,  9.3404e+08,  9.3963e+08,  9.4522e+08,  9.5082e+08,&
                 9.5641e+08,  9.6200e+08,  9.6759e+08,  9.7319e+08,  9.7878e+08,& 
                 9.8437e+08,  9.8997e+08,  9.9556e+08,  1.0012e+09,  1.0067e+09 /

  data T_hong / 1.0000e+04, 1.0000e+04, 1.0001e+04, 1.0002e+04, 1.0003e+04,&  
                1.0005e+04, 1.0008e+04, 1.0012e+04, 1.0018e+04, 1.0026e+04,&
                1.0036e+04, 1.0047e+04, 1.0058e+04, 1.0011e+04, 9.6700e+03,&
                9.3468e+03, 8.9643e+03, 8.6116e+03, 8.2357e+03, 7.8466e+03,&
                7.4734e+03, 7.1441e+03, 6.8397e+03, 6.5263e+03, 6.2978e+03,&
                6.1165e+03, 5.9274e+03, 5.7787e+03, 5.6233e+03, 5.4920e+03,&
                5.3684e+03, 5.2602e+03, 5.1621e+03, 5.0768e+03, 4.9988e+03,&
                4.9202e+03, 4.8435e+03, 4.7743e+03, 4.7107e+03, 4.6492e+03,&
                4.5903e+03, 4.5333e+03, 4.4783e+03, 4.4272e+03, 4.3832e+03,&
                4.3391e+03, 4.2948e+03, 4.2540e+03, 4.2277e+03, 4.2067e+03,&
                4.1940e+03, 4.1851e+03, 4.1956e+03, 4.2262e+03, 4.2723e+03,&
                4.3381e+03, 4.4204e+03, 4.5278e+03, 4.6445e+03, 4.7643e+03,&
                4.8856e+03, 5.0037e+03, 5.1176e+03, 5.2250e+03, 5.3211e+03,&
                5.4063e+03, 5.4888e+03, 5.5778e+03, 5.6443e+03, 5.7069e+03,&
                5.7603e+03, 5.8132e+03, 5.8700e+03, 5.9167e+03, 5.9592e+03,&
                5.9938e+03, 6.0256e+03, 6.0561e+03, 6.0805e+03, 6.1051e+03,&
                6.1298e+03, 6.1515e+03, 6.1720e+03, 6.1904e+03, 6.2096e+03,&
                6.2281e+03, 6.2458e+03, 6.2654e+03, 6.2877e+03, 6.3134e+03,&
                6.3429e+03, 6.3754e+03, 6.4108e+03, 6.4497e+03, 6.4949e+03,&
                6.5456e+03, 6.6066e+03, 6.6837e+03, 6.7794e+03, 6.8983e+03,&
                7.0383e+03, 7.1943e+03, 7.3537e+03, 7.4943e+03, 7.6133e+03,&
                7.7131e+03, 7.8049e+03, 7.8909e+03, 7.9677e+03, 8.0333e+03,&
                8.0851e+03, 8.1219e+03, 8.1534e+03, 8.2225e+03, 8.4087e+03,&
                8.7744e+03, 9.3187e+03, 1.0097e+04, 1.1128e+04, 1.2398e+04,&
                1.3941e+04, 1.5857e+04, 1.8286e+04, 2.1320e+04, 2.4676e+04,&
                2.7878e+04, 3.0947e+04, 3.4019e+04, 3.7182e+04, 4.0498e+04,&
                4.4021e+04, 4.7801e+04, 5.1890e+04, 5.6344e+04, 6.1230e+04,&
                6.6619e+04, 7.2587e+04, 7.9209e+04, 8.6564e+04, 9.4709e+04,&
                1.0364e+05, 1.1329e+05, 1.2363e+05, 1.3469e+05, 1.4653e+05,&
                1.5920e+05, 1.7275e+05, 1.8719e+05, 2.0250e+05, 2.1859e+05,&
                2.3535e+05, 2.5260e+05, 2.7013e+05, 2.8767e+05, 3.0498e+05,&
                3.2184e+05, 3.3811e+05, 3.5372e+05, 3.6862e+05, 3.8283e+05,&
                3.9636e+05, 4.0926e+05, 4.2157e+05, 4.3333e+05, 4.4459e+05,&
                4.5538e+05, 4.6575e+05, 4.7573e+05, 4.8535e+05, 4.9464e+05,&
                5.0362e+05, 5.1232e+05, 5.2076e+05, 5.2895e+05, 5.3691e+05,&
                5.4467e+05, 5.5222e+05, 5.5958e+05, 5.6677e+05, 5.7379e+05,&
                5.8066e+05, 5.8738e+05, 5.9395e+05, 6.0040e+05, 6.0672e+05,&
                6.1291e+05, 6.1900e+05, 6.2497e+05, 6.3084e+05, 6.3661e+05,&
                6.4228e+05, 6.4786e+05, 6.5335e+05, 6.5876e+05, 6.6408e+05,&
                6.6933e+05, 6.7450e+05, 6.7960e+05, 6.8463e+05, 6.8959e+05,&
                6.9449e+05, 6.9932e+05, 7.0409e+05, 7.0880e+05, 7.1346e+05,&
                7.1806e+05, 7.2260e+05, 7.2709e+05, 7.3154e+05, 7.3593e+05,&
                7.4027e+05, 7.4457e+05, 7.4882e+05, 7.5303e+05, 7.5719e+05,&
                7.6131e+05, 7.6539e+05, 7.6943e+05, 7.7343e+05, 7.7739e+05,&
                7.8131e+05, 7.8518e+05, 7.8898e+05, 7.9272e+05, 7.9641e+05,&
                8.0004e+05, 8.0362e+05, 8.0715e+05, 8.1063e+05, 8.1406e+05,&
                8.1745e+05, 8.2079e+05, 8.2409e+05, 8.2734e+05, 8.3056e+05,&
                8.3373e+05, 8.3687e+05, 8.3997e+05, 8.4303e+05, 8.4606e+05,&
                8.4906e+05, 8.5205e+05, 8.5502e+05, 8.5796e+05, 8.6090e+05,&
                8.6381e+05, 8.6670e+05, 8.6958e+05, 8.7244e+05, 8.7528e+05,&
                8.7810e+05, 8.8091e+05, 8.8371e+05, 8.8648e+05, 8.8924e+05,&
                8.9199e+05, 8.9472e+05, 8.9743e+05, 9.0013e+05, 9.0281e+05,&
                9.0548e+05, 9.0814e+05, 9.1078e+05, 9.1341e+05, 9.1602e+05,&
                9.1862e+05, 9.2121e+05, 9.2378e+05, 9.2634e+05, 9.2889e+05,&
                9.3142e+05, 9.3394e+05, 9.3645e+05, 9.3895e+05, 9.4143e+05,&
                9.4391e+05, 9.4637e+05, 9.4882e+05, 9.5126e+05, 9.5368e+05,&
                9.5610e+05, 9.5850e+05, 9.6089e+05, 9.6328e+05, 9.6565e+05,&
                9.6801e+05, 9.7036e+05, 9.7270e+05, 9.7503e+05, 9.7735e+05,&
                9.7965e+05, 9.8195e+05, 9.8424e+05, 9.8652e+05, 9.8879e+05,&
                9.9105e+05, 9.9330e+05, 9.9554e+05, 9.9778e+05, 9.9778e+05 /


contains

  subroutine get_atm_para(h,rho,pth,mg,nh,rho0,Tcurve)
    ! h [cm]
    ! rho [cm^-3]
    ! pth [dyne cm^-2]
    ! mg [g cm s^-2]
    ! nh -- number of points
    ! rho0 -- number density at h(1)
    ! Tcurve -- 'VAL-C' | 'Hong2017' | 'YS2001'

    integer :: nh
    double precision :: rho0
    double precision :: h(nh),rho(nh),pth(nh),mg(nh)
    character(20) :: Tcurve

    double precision :: Te(nh)
    integer :: j
    double precision :: invT,dh

    select case(Tcurve)
      case('VAL-C')
        call get_Te_VALC(h,Te,nh)
        if (mype==0) print *, 'VAL-C temperature curve'

      case('Hong2017')
        call get_Te_Hong(h,Te,nh)
        if (mype==0) print *, 'Temperature curve in Hong et al. (2017) is employed'

      case('YS2001')
        call get_NT_YS(h,rho,pth,nh,rho0)
        if (mype==0) print *, 'Temperature and density curves in Yokoyama & Shibata (2001) are employed'

      case default
        call mpistop("This temperature curve is unknown")

    end select


    ! density and pressure profiles
    if (Tcurve=='VAL-C' .or. Tcurve=='Hong2017') then
      rho(1)=rho0
      pth(1)=rho(1)*const_kb*Te(1)

      invT=0.d0
      do j=2,nh
        dh=h(j)-h(j-1)
        invT=invT+dh*(mg(j)/(const_kb*Te(j))+mg(j-1)/(const_kb*Te(j-1)))*0.5d0
        pth(j)=pth(1)*dexp(invT)
        rho(j)=pth(j)/(const_kb*Te(j))
      end do
    endif

  end subroutine get_atm_para

  subroutine get_Te_VALC(h,Te,nh)
    use mod_interpolation
    use mod_constants

    integer :: nh
    double precision :: h(nh),Te(nh)

    integer :: ih,j,imin,imax,n_table
    double precision :: h_table(n_valc),T_table(n_valc)
    double precision :: unit_h,unit_T,htra,Ttra,Fc,kappa,dTdh

    ! temperature profile
    unit_h=1.d5 !  km -> cm
    unit_T=1.0  !  K -> K
    Fc=6*2.d5
    kappa=8.d-7

    n_table=n_valc
    h_table=h_valc*unit_h
    T_table=T_valc*unit_T
    htra=h_table(n_table)
    Ttra=T_table(n_table)
    dTdh=(T_table(2)-T_table(1))/(h_table(2)-h_table(1))

    do ih=1,nh
      if (h(ih)>h_table(n_table)) then
      ! above transition region
        Te(ih)=(3.5d0*Fc*(h(ih)-htra)/kappa+Ttra**3.5d0)**(2.d0/7.d0)
      endif

      if (h(ih)<=h_table(1)) then
      ! below photosphere
        Te(ih)=(h(ih)-h_table(1))*dTdh+T_table(1)
      endif
    enddo


    ! inside the table
    imin=nh
    imax=nh-1
    if (h(1)>=h_table(1) .and. h(1)<=h_table(n_table)) then
      imin=1
    else
      do ih=2,nh
        if (h(ih-1)<h_table(1) .and. h(ih)>=h_table(1) .and. h(ih)<=h_table(n_table)) imin=ih
      enddo
    endif
    if (h(nh)>=h_table(1) .and. h(nh)<=h_table(n_table)) then
      imax=nh
    else
      do ih=1,nh-1
        if (h(ih)<=h_table(n_table) .and. h(ih+1)>h_table(n_table) .and. h(ih)>=h_table(1)) imax=ih
      enddo
    endif

    if (imin<=imax) then
      call interp_linear(h_table,T_table,n_valc,h(imin:imax),Te(imin:imax),imax-imin+1)
    endif

  end subroutine get_Te_VALC

  subroutine get_Te_Hong(h,Te,nh)
    use mod_interpolation
    use mod_constants

    integer :: nh
    double precision :: h(nh),Te(nh)

    integer :: ih,j,imin,imax,n_table
    double precision :: h_table(n_hong),T_table(n_hong)
    double precision :: dTdh

    n_table=n_hong
    h_table=h_hong
    T_table=T_hong

    do ih=1,nh
      if (h(ih)>=h_table(n_table)) then
      ! above max height of the table
        dTdh=(T_table(n_table)-T_table(n_table-1))/(h_table(n_table)-h_table(n_table-1))
        Te(ih)=(h(ih)-h_table(n_table))*dTdh+T_table(n_table)
      endif

      if (h(ih)<=h_table(1)) then
      ! below min height of the table
        dTdh=(T_table(2)-T_table(1))/(h_table(2)-h_table(1))
        Te(ih)=(h(ih)-h_table(1))*dTdh+T_table(1)
      endif
    enddo


    ! inside the table
    imin=nh
    imax=nh-1
    if (h(1)>=h_table(1) .and. h(1)<=h_table(n_table)) then
      imin=1
    else
      do ih=2,nh
        if (h(ih-1)<h_table(1) .and. h(ih)>=h_table(1) .and. h(ih)<=h_table(n_table)) imin=ih
      enddo
    endif
    if (h(nh)>=h_table(1) .and. h(nh)<=h_table(n_table)) then
      imax=nh
    else
      do ih=1,nh-1
        if (h(ih)<=h_table(n_table) .and. h(ih+1)>h_table(n_table) .and. h(ih)>=h_table(1)) imax=ih
      enddo
    endif

    if (imin<=imax) then
      call interp_linear(h_table,T_table,n_table,h(imin:imax),Te(imin:imax),imax-imin+1)
    endif

  end subroutine get_Te_Hong

  subroutine get_NT_YS(h,rho,pth,nh,rho0)

    integer :: nh
    double precision :: rho0
    double precision :: h(nh),rho(nh),pth(nh)

    integer :: ih
    double precision :: htra,wtra,rhoc

    htra=0.3d9 ! height of initial transition region
    wtra=0.06d9 ! width of initial transition region 
    rhoc=rho0/1d5
    pth=0.47/1.67

    do ih=1,nh
      rho(ih)=rhoc+(rho0-rhoc)*(1.d0-tanh((h(ih)-htra)/wtra))
    enddo

  end subroutine get_NT_YS

end module mod_solar_atmosphere

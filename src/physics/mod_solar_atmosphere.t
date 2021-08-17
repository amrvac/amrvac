!> User can use subroutine get_atm_para to generate 1D solar stmosphere.
!> User should provide heights (h), number density at h=0,
!> number of points (nh), and the gravity (grav) at each point. 
!> User can select temperature profile.
!> This subroutine will return density and pressure at each point.
module mod_solar_atmosphere
  use mod_global_parameters
  use mod_constants
  implicit none

  integer :: n_valc,n_hong,n_fontenla,n_alc7
  double precision :: h_valc(1:52),T_valc(1:52)
  double precision :: h_hong(1:300),T_hong(1:300)
  double precision :: h_fontenla(1:56),T_fontenla(1:56)
  double precision :: h_alc7(1:140),T_alc7(1:140)

  data n_alc7 / 140 /

  data h_alc7 /  -100.0,   -90.0,   -80.0,   -70.0,   -60.0, & 
                  -50.0,   -40.0,   -30.0,   -20.0,   -10.0, &
                    0.0,    10.0,    20.0,    35.0,    50.0, &
                   75.0,   100.0,   125.0,   150.0,   175.0, &
                  200.0,   250.0,   300.0,   350.0,   400.0, &
                  450.0,   490.0,   525.0,   560.0,   615.0, &
                  660.0,   700.0,   750.0,   800.0,   854.0, &
                  900.0,   946.0,   971.0,  1003.0,  1032.0, &
                 1065.0,  1101.0,  1143.0,  1214.0,  1299.0, &
                 1398.0,  1520.0,  1617.0,  1722.0,  1820.0, &
                 1894.0,  1946.0,  1989.0,  2024.0,  2055.0, &
                 2083.0,  2098.0,  2110.0,  2120.0,  2126.0, &
                 2130.0,  2132.0,  2134.0,  2136.0,  2138.0, &
                 2141.9,  2145.2,  2147.1,  2148.7,  2150.1, &
                 2151.2,  2152.0,  2152.9,  2153.5,  2154.1, &
                 2154.7,  2155.5,  2156.2,  2156.8,  2157.6, &
                 2158.2,  2158.8,  2159.4,  2160.0,  2160.6, &
                 2161.3,  2162.3,  2163.3,  2164.6,  2166.0, &
                 2167.7,  2170.0,  2172.3,  2175.3,  2178.7, &
                 2182.4,  2186.1,  2191.5,  2196.5,  2201.9, &
                 2207.7,  2213.7,  2220.0,  2227.7,  2235.3, &
                 2243.3,  2251.6,  2259.4,  2269.7,  2282.4, &
                 2294.4,  2308.6,  2325.8,  2346.1,  2378.2, &
                 2418.7,  2457.9,  2495.6,  2535.5,  2578.2, &
                 2628.3,  2688.1,  2761.5,  2844.4,  2970.8, &
                 3133.8,  3425.0,  3752.7,  4295.9,  4968.6, &
                 5762.8,  7360.8,  8974.2, 11596.2, 15392.0, &
                21133.1, 26676.6, 36079.5, 47009.3, 68084.4 /

  data T_alc7 /    9380,    9120,    8850,    8540,    8220, &
                   7900,    7590,    7280,    7020,    6780, &
                   6583,    6397,    6231,    6006,    5826, &
                   5607,    5431,    5288,    5165,    5080, &
                   5010,    4907,    4805,    4700,    4590, &
                   4485,    4435,    4410,    4400,    4435, &
                   4510,    4640,    4840,    5090,    5430, &
                   5720,    5969,    6100,    6225,    6315, &
                   6400,    6474,    6531,    6576,    6598, &
                   6610,    6623,    6633,    6643,    6652, &
                   6660,    6667,    6674,    6680,    6686, &
                   6694,    6700,    6706,    6718,    6740, &
                   6768,    6800,    6870,    6992,    7248, &
                   7950,    9115,   10980,   13200,   15760, &
                  18140,   20510,   23100,   25120,   27130, &
                  29500,   32260,   34580,   36870,   39400, &
                  41450,   43400,   45140,   46800,   48490, &
                  50140,   52690,   55020,   57790,   60790, &
                  63950,   68000,   71810,   76330,   81220, &
                  86120,   90640,   96860,  102300,  107800, &
                 113200,  118700,  124000,  130200,  135800, &
                 141600,  147200,  152300,  158600,  166000, &
                 172600,  180100,  188600,  198100,  212000, &
                 227800,  241900,  254400,  266700,  279000, &
                 292500,  307300,  324000,  341300,  365000, &
                 392000,  432900,  471600,  524300,  577400, &
                 629300,  711700,  778300,  865000,  996370, &
                1080000, 1170000, 1294000, 1410000, 1586000 / 



  data n_fontenla / 56 /

  data h_fontenla /  -100.0,   -90.0,   -80.0,   -70.0,   -60.0, &
                  -50.0,   -40.0,   -30.0,   -20.0,   -10.0, &
                   -0.7,     9.5,    20.0,    35.0,    50.0, &
                   75.0,   100.0,   125.0,   150.0,   175.0, &
                  200.0,   250.0,   300.0,   347.0,   398.0, &
                  447.5,   485.0,   519.0,   570.0,   625.0, &
                  696.0,   772.0,   812.0,   838.0,   874.5, &
                  899.5,   926.0,   947.0,   965.0,  1020.0, &
                 1060.0,  1120.0,  1180.0,  1245.0,  1319.0, &
                 1447.0,  1571.0,  1687.0,  1801.0,  1911.0, &
                 2013.0,  2142.0,  2263.0,  2357.0,  2425.0, &
                 2472.0 /

  data T_fontenla /    9460,    9220,    8940,    8620,    8280, & 
                   7940,    7600,    7270,    6960,    6690, &
                   6490,    6310,    6150,    5950,    5780, &
                   5550,    5380,    5245,    5130,    5035, &
                   4960,    4835,    4740,    4660,    4580, &
                   4510,    4450,    4390,    4300,    4200, &
                   4080,    3955,    3890,    3850,    3800, &
                   3810,    3950,    4150,    4400,    5350, &
                   5770,    6120,    6270,    6390,    6450, &
                   6490,    6530,    6570,    6620,    6670, &
                   6720,    6770,    6830,    6910,    7030, &
                   7210 /



  data n_valc / 52 /

  data    h_valc /    -75,    -50,    -25,      0,     50, & 
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

  data h_hong / -3.7616d+07, -3.7384d+07, -3.7045d+07, -3.6561d+07, -3.5880d+07,& 
                -3.4950d+07, -3.3718d+07, -3.2153d+07, -3.0250d+07, -2.8050d+07,&
                -2.5627d+07, -2.3088d+07, -2.0556d+07, -1.8151d+07, -1.5972d+07,&
                -1.4004d+07, -1.2234d+07, -1.0625d+07, -9.1510d+06, -7.7773d+06,&
                -6.4599d+06, -5.1584d+06, -3.8455d+06, -2.4946d+06, -1.0577d+06,&
                 4.7767d+05,  2.1041d+06,  3.8237d+06,  5.6198d+06,  7.4817d+06,&
                 9.3930d+06,  1.1339d+07,  1.3304d+07,  1.5278d+07,  1.7250d+07,&
                 1.9213d+07,  2.1163d+07,  2.3099d+07,  2.5019d+07,  2.6922d+07,&
                 2.8809d+07,  3.0678d+07,  3.2532d+07,  3.4371d+07,  3.6196d+07,&
                 3.8009d+07,  3.9811d+07,  4.1605d+07,  4.3396d+07,  4.5186d+07,&
                 4.6977d+07,  4.8773d+07,  5.0576d+07,  5.2390d+07,  5.4218d+07,&
                 5.6065d+07,  5.7934d+07,  5.9833d+07,  6.1766d+07,  6.3736d+07,&
                 6.5744d+07,  6.7791d+07,  6.9880d+07,  7.2011d+07,  7.4185d+07,&
                 7.6403d+07,  7.8660d+07,  8.0954d+07,  8.3284d+07,  8.5647d+07,&
                 8.8038d+07,  9.0456d+07,  9.2896d+07,  9.5357d+07,  9.7838d+07,&
                 1.0033d+08,  1.0285d+08,  1.0537d+08,  1.0791d+08,  1.1046d+08,&
                 1.1302d+08,  1.1559d+08,  1.1818d+08,  1.2077d+08,  1.2337d+08,&
                 1.2599d+08,  1.2862d+08,  1.3126d+08,  1.3392d+08,  1.3659d+08,&
                 1.3927d+08,  1.4198d+08,  1.4470d+08,  1.4744d+08,  1.5019d+08,&
                 1.5295d+08,  1.5571d+08,  1.5847d+08,  1.6121d+08,  1.6389d+08,&
                 1.6648d+08,  1.6895d+08,  1.7122d+08,  1.7325d+08,  1.7499d+08,&
                 1.7642d+08,  1.7754d+08,  1.7838d+08,  1.7900d+08,  1.7944d+08,&
                 1.7974d+08,  1.7996d+08,  1.8010d+08,  1.8020d+08,  1.8027d+08,&
                 1.8032d+08,  1.8035d+08,  1.8037d+08,  1.8039d+08,  1.8041d+08,&
                 1.8042d+08,  1.8044d+08,  1.8045d+08,  1.8046d+08,  1.8047d+08,&
                 1.8049d+08,  1.8051d+08,  1.8053d+08,  1.8055d+08,  1.8059d+08,&
                 1.8063d+08,  1.8068d+08,  1.8074d+08,  1.8082d+08,  1.8092d+08,&
                 1.8103d+08,  1.8118d+08,  1.8135d+08,  1.8157d+08,  1.8182d+08,&
                 1.8214d+08,  1.8252d+08,  1.8299d+08,  1.8357d+08,  1.8428d+08,&
                 1.8515d+08,  1.8621d+08,  1.8749d+08,  1.8904d+08,  1.9089d+08,&
                 1.9307d+08,  1.9561d+08,  1.9852d+08,  2.0180d+08,  2.0544d+08,&
                 2.0941d+08,  2.1368d+08,  2.1819d+08,  2.2291d+08,  2.2781d+08,&
                 2.3284d+08,  2.3798d+08,  2.4320d+08,  2.4849d+08,  2.5383d+08,&
                 2.5922d+08,  2.6463d+08,  2.7008d+08,  2.7554d+08,  2.8102d+08,&
                 2.8651d+08,  2.9202d+08,  2.9753d+08,  3.0305d+08,  3.0858d+08,&
                 3.1412d+08,  3.1966d+08,  3.2520d+08,  3.3075d+08,  3.3630d+08,&
                 3.4186d+08,  3.4741d+08,  3.5297d+08,  3.5854d+08,  3.6410d+08,&
                 3.6967d+08,  3.7523d+08,  3.8080d+08,  3.8637d+08,  3.9195d+08,&
                 3.9752d+08,  4.0309d+08,  4.0867d+08,  4.1424d+08,  4.1982d+08,&
                 4.2540d+08,  4.3098d+08,  4.3656d+08,  4.4214d+08,  4.4772d+08,&
                 4.5330d+08,  4.5888d+08,  4.6446d+08,  4.7004d+08,  4.7563d+08,&
                 4.8121d+08,  4.8679d+08,  4.9238d+08,  4.9796d+08,  5.0355d+08,&
                 5.0913d+08,  5.1472d+08,  5.2030d+08,  5.2589d+08,  5.3148d+08,&
                 5.3706d+08,  5.4265d+08,  5.4824d+08,  5.5382d+08,  5.5941d+08,&
                 5.6500d+08,  5.7059d+08,  5.7617d+08,  5.8176d+08,  5.8735d+08,&
                 5.9294d+08,  5.9853d+08,  6.0412d+08,  6.0971d+08,  6.1530d+08,&
                 6.2089d+08,  6.2648d+08,  6.3207d+08,  6.3766d+08,  6.4325d+08,&
                 6.4884d+08,  6.5443d+08,  6.6002d+08,  6.6561d+08,  6.7120d+08,&
                 6.7679d+08,  6.8238d+08,  6.8798d+08,  6.9357d+08,  6.9916d+08,&
                 7.0475d+08,  7.1034d+08,  7.1593d+08,  7.2152d+08,  7.2712d+08,&
                 7.3271d+08,  7.3830d+08,  7.4389d+08,  7.4948d+08,  7.5508d+08,&
                 7.6067d+08,  7.6626d+08,  7.7185d+08,  7.7744d+08,  7.8304d+08,&
                 7.8863d+08,  7.9422d+08,  7.9981d+08,  8.0540d+08,  8.1100d+08,&
                 8.1659d+08,  8.2218d+08,  8.2777d+08,  8.3337d+08,  8.3896d+08,&
                 8.4455d+08,  8.5014d+08,  8.5574d+08,  8.6133d+08,  8.6692d+08,&
                 8.7252d+08,  8.7811d+08,  8.8370d+08,  8.8929d+08,  8.9489d+08,&
                 9.0048d+08,  9.0607d+08,  9.1166d+08,  9.1726d+08,  9.2285d+08,&
                 9.2844d+08,  9.3404d+08,  9.3963d+08,  9.4522d+08,  9.5082d+08,&
                 9.5641d+08,  9.6200d+08,  9.6759d+08,  9.7319d+08,  9.7878d+08,& 
                 9.8437d+08,  9.8997d+08,  9.9556d+08,  1.0012d+09,  1.0067d+09 /

  data T_hong / 1.0000d+04, 1.0000d+04, 1.0001d+04, 1.0002d+04, 1.0003d+04,&  
                1.0005d+04, 1.0008d+04, 1.0012d+04, 1.0018d+04, 1.0026d+04,&
                1.0036d+04, 1.0047d+04, 1.0058d+04, 1.0011d+04, 9.6700d+03,&
                9.3468d+03, 8.9643d+03, 8.6116d+03, 8.2357d+03, 7.8466d+03,&
                7.4734d+03, 7.1441d+03, 6.8397d+03, 6.5263d+03, 6.2978d+03,&
                6.1165d+03, 5.9274d+03, 5.7787d+03, 5.6233d+03, 5.4920d+03,&
                5.3684d+03, 5.2602d+03, 5.1621d+03, 5.0768d+03, 4.9988d+03,&
                4.9202d+03, 4.8435d+03, 4.7743d+03, 4.7107d+03, 4.6492d+03,&
                4.5903d+03, 4.5333d+03, 4.4783d+03, 4.4272d+03, 4.3832d+03,&
                4.3391d+03, 4.2948d+03, 4.2540d+03, 4.2277d+03, 4.2067d+03,&
                4.1940d+03, 4.1851d+03, 4.1956d+03, 4.2262d+03, 4.2723d+03,&
                4.3381d+03, 4.4204d+03, 4.5278d+03, 4.6445d+03, 4.7643d+03,&
                4.8856d+03, 5.0037d+03, 5.1176d+03, 5.2250d+03, 5.3211d+03,&
                5.4063d+03, 5.4888d+03, 5.5778d+03, 5.6443d+03, 5.7069d+03,&
                5.7603d+03, 5.8132d+03, 5.8700d+03, 5.9167d+03, 5.9592d+03,&
                5.9938d+03, 6.0256d+03, 6.0561d+03, 6.0805d+03, 6.1051d+03,&
                6.1298d+03, 6.1515d+03, 6.1720d+03, 6.1904d+03, 6.2096d+03,&
                6.2281d+03, 6.2458d+03, 6.2654d+03, 6.2877d+03, 6.3134d+03,&
                6.3429d+03, 6.3754d+03, 6.4108d+03, 6.4497d+03, 6.4949d+03,&
                6.5456d+03, 6.6066d+03, 6.6837d+03, 6.7794d+03, 6.8983d+03,&
                7.0383d+03, 7.1943d+03, 7.3537d+03, 7.4943d+03, 7.6133d+03,&
                7.7131d+03, 7.8049d+03, 7.8909d+03, 7.9677d+03, 8.0333d+03,&
                8.0851d+03, 8.1219d+03, 8.1534d+03, 8.2225d+03, 8.4087d+03,&
                8.7744d+03, 9.3187d+03, 1.0097d+04, 1.1128d+04, 1.2398d+04,&
                1.3941d+04, 1.5857d+04, 1.8286d+04, 2.1320d+04, 2.4676d+04,&
                2.7878d+04, 3.0947d+04, 3.4019d+04, 3.7182d+04, 4.0498d+04,&
                4.4021d+04, 4.7801d+04, 5.1890d+04, 5.6344d+04, 6.1230d+04,&
                6.6619d+04, 7.2587d+04, 7.9209d+04, 8.6564d+04, 9.4709d+04,&
                1.0364d+05, 1.1329d+05, 1.2363d+05, 1.3469d+05, 1.4653d+05,&
                1.5920d+05, 1.7275d+05, 1.8719d+05, 2.0250d+05, 2.1859d+05,&
                2.3535d+05, 2.5260d+05, 2.7013d+05, 2.8767d+05, 3.0498d+05,&
                3.2184d+05, 3.3811d+05, 3.5372d+05, 3.6862d+05, 3.8283d+05,&
                3.9636d+05, 4.0926d+05, 4.2157d+05, 4.3333d+05, 4.4459d+05,&
                4.5538d+05, 4.6575d+05, 4.7573d+05, 4.8535d+05, 4.9464d+05,&
                5.0362d+05, 5.1232d+05, 5.2076d+05, 5.2895d+05, 5.3691d+05,&
                5.4467d+05, 5.5222d+05, 5.5958d+05, 5.6677d+05, 5.7379d+05,&
                5.8066d+05, 5.8738d+05, 5.9395d+05, 6.0040d+05, 6.0672d+05,&
                6.1291d+05, 6.1900d+05, 6.2497d+05, 6.3084d+05, 6.3661d+05,&
                6.4228d+05, 6.4786d+05, 6.5335d+05, 6.5876d+05, 6.6408d+05,&
                6.6933d+05, 6.7450d+05, 6.7960d+05, 6.8463d+05, 6.8959d+05,&
                6.9449d+05, 6.9932d+05, 7.0409d+05, 7.0880d+05, 7.1346d+05,&
                7.1806d+05, 7.2260d+05, 7.2709d+05, 7.3154d+05, 7.3593d+05,&
                7.4027d+05, 7.4457d+05, 7.4882d+05, 7.5303d+05, 7.5719d+05,&
                7.6131d+05, 7.6539d+05, 7.6943d+05, 7.7343d+05, 7.7739d+05,&
                7.8131d+05, 7.8518d+05, 7.8898d+05, 7.9272d+05, 7.9641d+05,&
                8.0004d+05, 8.0362d+05, 8.0715d+05, 8.1063d+05, 8.1406d+05,&
                8.1745d+05, 8.2079d+05, 8.2409d+05, 8.2734d+05, 8.3056d+05,&
                8.3373d+05, 8.3687d+05, 8.3997d+05, 8.4303d+05, 8.4606d+05,&
                8.4906d+05, 8.5205d+05, 8.5502d+05, 8.5796d+05, 8.6090d+05,&
                8.6381d+05, 8.6670d+05, 8.6958d+05, 8.7244d+05, 8.7528d+05,&
                8.7810d+05, 8.8091d+05, 8.8371d+05, 8.8648d+05, 8.8924d+05,&
                8.9199d+05, 8.9472d+05, 8.9743d+05, 9.0013d+05, 9.0281d+05,&
                9.0548d+05, 9.0814d+05, 9.1078d+05, 9.1341d+05, 9.1602d+05,&
                9.1862d+05, 9.2121d+05, 9.2378d+05, 9.2634d+05, 9.2889d+05,&
                9.3142d+05, 9.3394d+05, 9.3645d+05, 9.3895d+05, 9.4143d+05,&
                9.4391d+05, 9.4637d+05, 9.4882d+05, 9.5126d+05, 9.5368d+05,&
                9.5610d+05, 9.5850d+05, 9.6089d+05, 9.6328d+05, 9.6565d+05,&
                9.6801d+05, 9.7036d+05, 9.7270d+05, 9.7503d+05, 9.7735d+05,&
                9.7965d+05, 9.8195d+05, 9.8424d+05, 9.8652d+05, 9.8879d+05,&
                9.9105d+05, 9.9330d+05, 9.9554d+05, 9.9778d+05, 9.9778d+05 /


contains

  subroutine get_atm_para(h,rho,pth,grav,nh,Tcurve,hc,rhohc)
    use mod_global_parameters
    ! input:h,grav,nh,rho0,Tcurve; output:rho,pth (dimensionless units)
    ! nh -- number of points
    ! rho0 -- number density at h=0
    ! Tcurve -- 'VAL-C' | 'Hong2017' | 'SPRM305' | 'AL-C7'

    integer, intent(in) :: nh
    double precision :: h(nh),rho(nh),pth(nh),grav(nh)
    double precision :: rhohc,hc
    character(*) :: Tcurve

    double precision :: h_cgs(nh),Te_cgs(nh),Te(nh)
    integer :: j
    double precision :: invT,dh,rhot,dht,ratio

    h_cgs=h*unit_length

    select case(Tcurve)
      case('VAL-C')
        call get_Te_VALC(h_cgs,Te_cgs,nh)
        if (mype==0) print *, 'Temperature curve from Vernazza, Avrett & Loeser, 1981, ApJS, 45, 635'

      case('Hong2017')
        call get_Te_Hong(h_cgs,Te_cgs,nh)
        if (mype==0) print *, 'Temperature curve from Hong et al. 2017, ApJ, 845, 144'

      case('Fontenla')
        call get_Te_SPRM(h_cgs,Te_cgs,nh)
        if (mype==0) print *, 'Temperature curve from Fontenla et al. 2007, ApJ, 667, 1243'

      case('AL-C7')
        call get_Te_ALC7(h_cgs,Te_cgs,nh)
        if (mype==0) print *, 'Temperature curve from Avrett & Loeser 2008, ApJS, 175, 229'

      case default
        call mpistop("Unknown temperature curve")

    end select

    Te=Te_cgs/unit_temperature

    ! density and pressure profiles
    rho(1)=1.d5
    pth(1)=rho(1)*Te(1)

    rhot=1.d0
    invT=0.d0
    do j=2,nh
      dh=h(j)-h(j-1)
      invT=invT+dh*(grav(j)/Te(j)+grav(j-1)/Te(j-1))*0.5d0
      pth(j)=pth(1)*dexp(invT)
      rho(j)=pth(j)/Te(j)

      if (h(j-1)<=hc .and. h(j)>hc) then
        dht=hc-h(j-1)
        rhot=rho(j-1)+dht*(rho(j)-rho(j-1))/(h(j)-h(j-1))
      endif
    end do

    ratio=rhohc/rhot
    rho=rho*ratio
    pth=pth*ratio

  end subroutine get_atm_para

  subroutine get_Te_ALC7(h,Te,nh)
    use mod_interpolation
    use mod_constants

    integer :: nh
    double precision :: h(nh),Te(nh)

    integer :: ih,j,imin,imax,n_table
    double precision :: h_table(n_alc7),T_table(n_alc7)
    double precision :: unit_h,unit_T,htra,Ttra,Fc,kappa,dTdh

    ! temperature profile
    unit_h=1.d5 !  km -> cm
    unit_T=1.0  !  K -> K
    Fc=1.72d5
    kappa=8.d-7

    n_table=n_alc7
    h_table=h_alc7*unit_h
    T_table=T_alc7*unit_T
    htra=2196.52*unit_h
    Ttra=1.023d5
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
      call interp_linear(h_table,T_table,n_table,h(imin:imax),Te(imin:imax),imax-imin+1)
    endif

  end subroutine get_Te_ALC7

  subroutine get_Te_SPRM(h,Te,nh)
    use mod_interpolation
    use mod_constants

    integer :: nh
    double precision :: h(nh),Te(nh)

    integer :: ih,j,imin,imax,n_table
    double precision :: h_table(n_fontenla),T_table(n_fontenla)
    double precision :: unit_h,unit_T,dTdh
    double precision :: h1,h2,h3
    double precision :: Tpho,Ttop,htanh,wtra
    double precision :: htra,Ttra,Fc,kappa

    unit_h=1.d5 !  km -> cm
    unit_T=1.0  !  K -> K

    n_table=n_fontenla
    h_table=h_fontenla*unit_h
    T_table=T_fontenla*unit_T

    ! height for shift curve table/function
    h1=h_table(1)
    h2=h_table(n_table)
    h3=3271.d0*unit_h

    ! parameter for T curve in (h2,h3)
    htanh=3464.d0*unit_h
    Tpho=6726.d0
    Ttop=1.5d6
    wtra=246.d0*unit_h
    
    ! parameter for T curve above h3
    htra=3271.d0*unit_h
    Ttra=2.642d5
    kappa=8.d-7
    Fc=4.791d5

    ! parameter for T curve below h1
    dTdh=(T_table(2)-T_table(1))/(h_table(2)-h_table(1))


    do ih=1,nh
      ! below photosphere
      if (h(ih)<=h1) then
        Te(ih)=(h(ih)-h_table(1))*dTdh+T_table(1)
      endif

      ! high chromosphere and low transition region
      if (h(ih)>h2 .and. h(ih)<h3) then
        Te(ih)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((h(ih)-htanh)/wtra)+1.d0)
      endif

      ! high transition region and corona
      if (h(ih)>=h3) then
        Te(ih)=(3.5d0*Fc*(h(ih)-htra)/kappa+Ttra**3.5)**(2.d0/7.d0)
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
      call interp_linear(h_table,T_table,n_fontenla,h(imin:imax),Te(imin:imax),imax-imin+1)
    endif

  end subroutine 

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
    Fc=5.38d5
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
    double precision :: htra,Ttra,Fc,kappa

    n_table=n_hong
    h_table=h_hong
    T_table=T_hong

    htra=h_table(150)
    Ttra=T_table(150)   
    Fc=2.76d5 ! heat flux in high corona
    kappa=8.d-7

    do ih=1,nh
      if (h(ih)>=h_table(n_table)) then
      ! above max height of the table
        Te(ih)=(3.5d0*Fc*(h(ih)-htra)/kappa+Ttra**3.5d0)**(2.d0/7.d0)
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

end module mod_solar_atmosphere

! module mod_thermal_emission -- synthesize emission flux of some
! thermal lines
! EUV lines database: 
! 'He_II_304' 'Fe_IX_171' 'Fe_XXIV_193' 'Fe_XIV_211' 'Fe_XVI_335'
! 'Fe_XVIII_94' 'Fe_XXI_131'
! subroutines: 
! get_EUV: get local EUV emission intensity (for 1d, 2d and 3d)
! get_SXR: get local Soft X-ray emission intensity (for 1d, 2d and 3d)

module mod_thermal_emission
  use mod_global_parameters
  use mod_geometry
  use mod_physics
  use mod_comm_lib, only: mpistop

  implicit none

  integer :: n_aia
  double precision :: t_aia(1:101)
  double precision :: f_94(1:101),f_131(1:101),f_171(1:101)
  double precision :: f_193(1:101),f_211(1:101),f_304(1:101)
  double precision :: f_335(1:101)
  integer :: n_iris
  double precision :: t_iris(1:41)
  double precision :: f_1354(1:41)
  integer :: n_eis
  double precision :: t_eis1(1:60),t_eis2(1:60)
  double precision :: f_263(1:60),f_264(1:60),f_192(1:60),f_255(1:60)


  double precision :: vec_xI1(1:3),vec_xI2(1:3),vec_LOS(1:3)

  data n_aia / 101 /

  data t_aia / 4. ,  4.05, 4.1,  4.15, 4.2,  4.25, 4.3,  4.35, &
               4.4,  4.45, 4.5,  4.55, 4.6,  4.65, 4.7,  4.75, &
               4.8,  4.85, 4.9,  4.95, 5. ,  5.05, 5.1,  5.15, &
               5.2,  5.25, 5.3,  5.35, 5.4,  5.45, 5.5,  5.55, &
               5.6,  5.65, 5.7,  5.75, 5.8,  5.85, 5.9,  5.95, &
               6. ,  6.05, 6.1,  6.15, 6.2,  6.25, 6.3,  6.35, &
               6.4,  6.45, 6.5,  6.55, 6.6,  6.65, 6.7,  6.75, &
               6.8,  6.85, 6.9,  6.95, 7. ,  7.05, 7.1,  7.15, &
               7.2,  7.25, 7.3,  7.35, 7.4,  7.45, 7.5,  7.55, &
               7.6,  7.65, 7.7,  7.75, 7.8,  7.85, 7.9,  7.95, &
               8. ,  8.05, 8.1,  8.15, 8.2,  8.25, 8.3,  8.35, &
               8.4,  8.45, 8.5,  8.55, 8.6,  8.65, 8.7,  8.75, &
               8.8,  8.85, 8.9,  8.95, 9. /

  data f_94 / 4.25022959d-37, 4.35880298d-36, 3.57054296d-35, 2.18175426d-34, & 
              8.97592571d-34, 2.68512961d-33, 7.49559346d-33, 2.11603751d-32, &
              5.39752853d-32, 1.02935904d-31, 1.33822307d-31, 1.40884290d-31, &
              1.54933156d-31, 2.07543102d-31, 3.42026227d-31, 6.31171444d-31, &
              1.16559416d-30, 1.95360497d-30, 2.77818735d-30, 3.43552578d-30, &
              4.04061803d-30, 4.75470982d-30, 5.65553769d-30, 6.70595782d-30, &
              7.80680354d-30, 8.93247715d-30, 1.02618156d-29, 1.25979030d-29, &
              1.88526483d-29, 3.62448572d-29, 7.50553279d-29, 1.42337571d-28, &
              2.37912813d-28, 3.55232305d-28, 4.84985757d-28, 6.20662827d-28, &
              7.66193687d-28, 9.30403645d-28, 1.10519802d-27, 1.25786927d-27, &
              1.34362634d-27, 1.33185242d-27, 1.22302081d-27, 1.05677973d-27, &
              9.23064720d-28, 8.78570994d-28, 8.02397416d-28, 5.87681142d-28, &
              3.82272695d-28, 3.11492649d-28, 3.85736090d-28, 5.98893519d-28, &
              9.57553548d-28, 1.46650267d-27, 2.10365847d-27, 2.79406671d-27, &
              3.39420087d-27, 3.71077520d-27, 3.57296767d-27, 2.95114380d-27, &
              2.02913103d-27, 1.13361825d-27, 5.13405629d-28, 2.01305089d-28, &
              8.15781482d-29, 4.28366817d-29, 3.08701543d-29, 2.68693906d-29, &
              2.51764203d-29, 2.41773103d-29, 2.33996083d-29, 2.26997246d-29, &
              2.20316143d-29, 2.13810001d-29, 2.07424438d-29, 2.01149189d-29, &
              1.94980213d-29, 1.88917920d-29, 1.82963583d-29, 1.77116920d-29, &
              1.71374392d-29, 1.65740593d-29, 1.60214447d-29, 1.54803205d-29, &
              1.49510777d-29, 1.44346818d-29, 1.39322305d-29, 1.34441897d-29, &
              1.29713709d-29, 1.25132618d-29, 1.20686068d-29, 1.14226584d-29, &
              1.09866413d-29, 1.05635524d-29, 1.01532444d-29, 9.75577134d-30, &
              9.37102736d-30, 8.99873335d-30, 8.63860172d-30, 8.29051944d-30, &
              7.95414793d-30 /

  data f_131 / 3.18403601d-37,   3.22254703d-36,   2.61657920d-35, &
               1.59575286d-34,   6.65779556d-34,   2.07015132d-33, &
               6.05768615d-33,   1.76074833d-32,   4.52633001d-32, &
               8.57121883d-32,   1.09184271d-31,   1.10207963d-31, &
               1.11371658d-31,   1.29105226d-31,   1.80385897d-31, &
               3.27295431d-31,   8.92002136d-31,   3.15214579d-30, &
               9.73440787d-30,   2.22709702d-29,   4.01788984d-29, &
               6.27471832d-29,   8.91764995d-29,   1.18725647d-28, &
               1.52888040d-28,   2.05082946d-28,   3.47651873d-28, &
               8.80482184d-28,   2.66533063d-27,   7.05805149d-27, &
               1.46072515d-26,   2.45282476d-26,   3.55303726d-26, &
               4.59075911d-26,   5.36503515d-26,   5.68444094d-26, &
               5.47222296d-26,   4.81119761d-26,   3.85959059d-26, &
               2.80383406d-26,   1.83977650d-26,   1.11182849d-26, &
               6.50748885d-27,   3.96843481d-27,   2.61876319d-27, &
               1.85525324d-27,   1.39717024d-27,   1.11504283d-27, &
               9.38169611d-28,   8.24801234d-28,   7.43331919d-28, &
               6.74537063d-28,   6.14495760d-28,   5.70805277d-28, &
               5.61219786d-28,   6.31981777d-28,   9.19747307d-28, &
               1.76795732d-27,   3.77985446d-27,   7.43166191d-27, &
               1.19785603d-26,   1.48234676d-26,   1.36673114d-26, &
               9.61047146d-27,   5.61209353d-27,   3.04779780d-27, &
               1.69378976d-27,   1.02113491d-27,   6.82223774d-28, &
               5.02099099d-28,   3.99377760d-28,   3.36279037d-28, &
               2.94767378d-28,   2.65740865d-28,   2.44396277d-28, &
               2.28003967d-28,   2.14941419d-28,   2.04178995d-28, &
               1.95031045d-28,   1.87011994d-28,   1.79777869d-28, &
               1.73093957d-28,   1.66795789d-28,   1.60785455d-28, &
               1.55002399d-28,   1.49418229d-28,   1.44022426d-28, &
               1.38807103d-28,   1.33772767d-28,   1.28908404d-28, &
               1.24196208d-28,   1.17437501d-28,   1.12854330d-28, &
               1.08410498d-28,   1.04112003d-28,   9.99529904d-29, &
               9.59358806d-29,   9.20512291d-29,   8.83009123d-29, &
               8.46817043d-29,   8.11921928d-29 /

  data f_171 / 2.98015581d-42, 1.24696230d-40, 3.37614652d-39, 5.64103034d-38, &
               5.20550266d-37, 2.77785939d-36, 1.16283616d-35, 6.50007689d-35, &
               9.96177399d-34, 1.89586076d-32, 2.10982799d-31, 1.36946479d-30, &
               6.27396553d-30, 2.29955134d-29, 7.13430211d-29, 1.91024282d-28, &
               4.35358848d-28, 7.94807808d-28, 1.07431875d-27, 1.08399488d-27, &
               9.16212938d-28, 7.34715770d-28, 6.59246382d-28, 9.13541375d-28, &
               2.05939035d-27, 5.08206555d-27, 1.10148083d-26, 2.01884662d-26, &
               3.13578384d-26, 4.14367719d-26, 5.36067711d-26, 8.74170213d-26, &
               1.64161233d-25, 2.94587860d-25, 4.76298332d-25, 6.91765639d-25, &
               9.08825111d-25, 1.08496183d-24, 1.17440114d-24, 1.13943939d-24, &
               9.71696981d-25, 7.09593688d-25, 4.31376399d-25, 2.12708486d-25, &
               8.47429567d-26, 3.17608104d-26, 1.95898842d-26, 1.98064242d-26, &
               1.67706555d-26, 8.99126003d-27, 3.29773878d-27, 1.28896127d-27, &
               8.51169698d-28, 7.53520167d-28, 6.18268143d-28, 4.30034650d-28, &
               2.78152409d-28, 1.95437088d-28, 1.65896278d-28, 1.68740181d-28, &
               1.76054383d-28, 1.63978419d-28, 1.32880591d-28, 1.00833205d-28, &
               7.82252806d-29, 6.36181741d-29, 5.34633869d-29, 4.58013864d-29, &
               3.97833422d-29, 3.49414760d-29, 3.09790940d-29, 2.76786227d-29, &
               2.48806269d-29, 2.24823367d-29, 2.04016653d-29, 1.85977413d-29, &
               1.70367499d-29, 1.56966125d-29, 1.45570643d-29, 1.35964565d-29, &
               1.27879263d-29, 1.21016980d-29, 1.15132499d-29, 1.09959628d-29, &
               1.05307482d-29, 1.01040261d-29, 9.70657096d-30, 9.33214234d-30, &
               8.97689427d-30, 8.63761192d-30, 8.31149879d-30, 7.85162401d-30, &
               7.53828281d-30, 7.23559452d-30, 6.94341530d-30, 6.66137038d-30, &
               6.38929156d-30, 6.12669083d-30, 5.87346434d-30, 5.62943622d-30, & 
               5.39435202d-30 /

  data f_193 / 6.40066486d-32, 4.92737300d-31, 2.95342934d-30, 1.28061594d-29, & 
               3.47747667d-29, 5.88554792d-29, 7.72171179d-29, 9.75609282d-29, &
               1.34318963d-28, 1.96252638d-28, 2.70163878d-28, 3.63192965d-28, &
               5.28087341d-28, 8.37821446d-28, 1.39089159d-27, 2.31749718d-27, &
               3.77510689d-27, 5.85198594d-27, 8.26021568d-27, 1.04870405d-26, &
               1.25209374d-26, 1.47406787d-26, 1.77174067d-26, 2.24098537d-26, &
               3.05926105d-26, 4.50018853d-26, 6.84720216d-26, 1.00595861d-25, &
               1.30759222d-25, 1.36481773d-25, 1.15943558d-25, 1.01467304d-25, &
               1.04092532d-25, 1.15071251d-25, 1.27416033d-25, 1.38463476d-25, &
               1.47882726d-25, 1.57041238d-25, 1.69786224d-25, 1.94970397d-25, &
               2.50332918d-25, 3.58321431d-25, 5.18061550d-25, 6.60405549d-25, &
               6.64085365d-25, 4.83825816d-25, 2.40545020d-25, 8.59534098d-26, &
               2.90920638d-26, 1.33204845d-26, 9.03933926d-27, 7.78910836d-27, &
               7.29342321d-27, 7.40267022d-27, 8.05279981d-27, 8.13829291d-27, &
               6.92634262d-27, 5.12521880d-27, 3.59527615d-27, 2.69617560d-27, &
               2.84432713d-27, 5.06697306d-27, 1.01281903d-26, 1.63526978d-26, &
               2.06759342d-26, 2.19482312d-26, 2.10050611d-26, 1.89837248d-26, &
               1.66347131d-26, 1.43071097d-26, 1.21518419d-26, 1.02078343d-26, &
               8.46936184d-27, 6.93015742d-27, 5.56973237d-27, 4.38951754d-27, &
               3.38456457d-27, 2.55309556d-27, 1.88904224d-27, 1.38057546d-27, &
               1.00718330d-27, 7.43581116d-28, 5.63562931d-28, 4.43359435d-28, &
               3.63923535d-28, 3.11248143d-28, 2.75586846d-28, 2.50672237d-28, &
               2.32419348d-28, 2.18325682d-28, 2.06834486d-28, 1.93497044d-28, &
               1.84540751d-28, 1.76356504d-28, 1.68741425d-28, 1.61566157d-28, &
               1.54754523d-28, 1.48249410d-28, 1.42020176d-28, 1.36045230d-28, &
               1.30307965d-28 /

  data f_211 / 4.74439912d-42, 1.95251522d-40, 5.19700194d-39, 8.53120166d-38, &
               7.72745727d-37, 4.04158559d-36, 1.64853511d-35, 8.56295439d-35, &
               1.17529722d-33, 2.16867729d-32, 2.40472264d-31, 1.56418133d-30, &
               7.20032889d-30, 2.65838271d-29, 8.33196904d-29, 2.26128236d-28, &
               5.24295811d-28, 9.77791121d-28, 1.35913489d-27, 1.43957785d-27, &
               1.37591544d-27, 1.49029886d-27, 2.06183401d-27, 3.31440622d-27, &
               5.42497318d-27, 8.41100374d-27, 1.17941366d-26, 1.49269794d-26, &
               1.71506074d-26, 1.71266353d-26, 1.51434781d-26, 1.36766622d-26, &
               1.33483562d-26, 1.36834518d-26, 1.45829002d-26, 1.62575306d-26, &
               1.88773347d-26, 2.22026986d-26, 2.54930499d-26, 2.80758138d-26, &
               3.06176409d-26, 3.62799792d-26, 5.13226109d-26, 8.46260744d-26, &
               1.38486586d-25, 1.86192535d-25, 1.78007934d-25, 1.16548409d-25, &
               5.89293257d-26, 2.69952884d-26, 1.24891081d-26, 6.41273176d-27, &
               4.08282914d-27, 3.26463328d-27, 2.76230280d-27, 2.08986882d-27, &
               1.37658470d-27, 8.48489381d-28, 5.19304217d-28, 3.19312514d-28, &
               2.02968197d-28, 1.50171666d-28, 1.39164218d-28, 1.42448821d-28, &
               1.41714519d-28, 1.33341059d-28, 1.20759270d-28, 1.07259692d-28, &
               9.44895400d-29, 8.29030041d-29, 7.25440631d-29, 6.33479483d-29, &
               5.51563757d-29, 4.79002469d-29, 4.14990482d-29, 3.59384972d-29, &
               3.12010860d-29, 2.72624742d-29, 2.40734791d-29, 2.15543565d-29, &
               1.95921688d-29, 1.80682882d-29, 1.68695662d-29, 1.59020936d-29, &
               1.50940886d-29, 1.43956179d-29, 1.37731622d-29, 1.32049043d-29, &
               1.26771875d-29, 1.21803879d-29, 1.17074716d-29, 1.10507836d-29, &
               1.06022834d-29, 1.01703080d-29, 9.75436986d-30, 9.35349257d-30, &
               8.96744546d-30, 8.59527489d-30, 8.23678940d-30, 7.89144480d-30, & 
               7.55891138d-30 /

  data f_304 / 3.62695850d-32, 2.79969087d-31, 1.68340584d-30, 7.32681440d-30, &
               1.99967770d-29, 3.41296785d-29, 4.55409104d-29, 5.94994635d-29, &
               8.59864963d-29, 1.39787633d-28, 3.17701965d-28, 1.14474920d-27, &
               4.44845958d-27, 1.54785841d-26, 4.70265345d-26, 1.24524365d-25, &
               2.81535352d-25, 5.10093666d-25, 6.83545307d-25, 6.82110329d-25, &
               5.66886188d-25, 4.36205513d-25, 3.29265688d-25, 2.49802368d-25, &
               1.92527113d-25, 1.51058572d-25, 1.20596047d-25, 9.76884267d-26, &
               7.89979266d-26, 6.18224289d-26, 4.67298332d-26, 3.57934505d-26, &
               2.84535785d-26, 2.32853022d-26, 1.95228514d-26, 1.67880071d-26, &
               1.47608785d-26, 1.32199691d-26, 1.20070960d-26, 1.09378177d-26, &
               1.00031730d-26, 9.62434001d-27, 1.05063954d-26, 1.27267143d-26, &
               1.45923057d-26, 1.36746707d-26, 1.03466970d-26, 6.97647829d-27, &
               4.63141039d-27, 3.19031994d-27, 2.33373613d-27, 1.81589079d-27, &
               1.48446917d-27, 1.26611478d-27, 1.12617468d-27, 1.03625148d-27, &
               9.61400595d-28, 8.79016231d-28, 7.82612130d-28, 6.73762960d-28, &
               5.59717956d-28, 4.53010243d-28, 3.65712196d-28, 3.00958686d-28, &
               2.54011502d-28, 2.18102277d-28, 1.88736437d-28, 1.63817539d-28, &
               1.42283147d-28, 1.23631916d-28, 1.07526003d-28, 9.36797928d-29, &
               8.18565660d-29, 7.18152734d-29, 6.32523238d-29, 5.59513985d-29, &
               4.96614048d-29, 4.42518826d-29, 3.95487628d-29, 3.54690294d-29, &
               3.18953930d-29, 2.87720933d-29, 2.60186750d-29, 2.36011522d-29, &
               2.14717806d-29, 1.95905217d-29, 1.79287981d-29, 1.64562262d-29, &
               1.51489425d-29, 1.39876064d-29, 1.29496850d-29, 1.18665438d-29, &
               1.10240474d-29, 1.02643099d-29, 9.57780996d-30, 8.95465151d-30, &
               8.38950190d-30, 7.87283711d-30, 7.40136507d-30, 6.96804279d-30, & 
               6.56945323d-30 /

  data f_335 / 2.46882661d-32, 1.89476632d-31, 1.13216502d-30, 4.89532008d-30, & 
               1.32745970d-29, 2.25390335d-29, 3.00511672d-29, 3.96035934d-29, &
               5.77977656d-29, 8.58600736d-29, 1.14083000d-28, 1.48644411d-28, &
               2.15788823d-28, 3.51628877d-28, 6.12200698d-28, 1.08184987d-27, &
               1.85590697d-27, 2.91679107d-27, 3.94405396d-27, 4.63610680d-27, &
               5.13824456d-27, 5.66602209d-27, 6.30009232d-27, 7.03422868d-27, &
               7.77973918d-27, 8.32371831d-27, 8.56724316d-27, 8.62601374d-27, &
               8.13308844d-27, 6.53188216d-27, 4.55197029d-27, 3.57590087d-27, &
               3.59571707d-27, 4.03502770d-27, 4.54366411d-27, 4.96914990d-27, &
               5.24601170d-27, 5.39979250d-27, 5.43023669d-27, 5.26235042d-27, &
               4.91585495d-27, 4.52628362d-27, 4.13385020d-27, 3.67538967d-27, &
               3.39939742d-27, 3.81284533d-27, 5.02332701d-27, 6.19438602d-27, &
               6.49613071d-27, 6.04010475d-27, 5.24664275d-27, 4.37225997d-27, &
               3.52957182d-27, 2.76212276d-27, 2.08473158d-27, 1.50850518d-27, &
               1.04602472d-27, 7.13091243d-28, 5.34289645d-28, 5.21079581d-28, &
               6.22246365d-28, 6.99555864d-28, 6.29665489d-28, 4.45077026d-28, &
               2.67046793d-28, 1.52774686d-28, 9.18061770d-29, 6.09116074d-29, &
               4.48562572d-29, 3.59463696d-29, 3.05820218d-29, 2.70766652d-29, &
               2.46144034d-29, 2.27758450d-29, 2.13331183d-29, 2.01537836d-29, &
               1.91566180d-29, 1.82893912d-29, 1.75167748d-29, 1.68136168d-29, &
               1.61615595d-29, 1.55481846d-29, 1.49643236d-29, 1.44046656d-29, &
               1.38657085d-29, 1.33459068d-29, 1.28447380d-29, 1.23615682d-29, &
               1.18963296d-29, 1.14478976d-29, 1.10146637d-29, 1.04039479d-29, &
               9.98611410d-30, 9.58205147d-30, 9.19202009d-30, 8.81551313d-30, &
               8.45252127d-30, 8.10224764d-30, 7.76469090d-30, 7.43954323d-30, &
               7.12653873d-30 /


  data n_iris / 41 /

  data t_iris / 4.        , 4.1       , 4.2       , 4.3       , 4.40000001, &
                4.50000001, 4.60000001, 4.70000001, 4.80000001, 4.90000001, &
                5.00000001, 5.10000002, 5.20000002, 5.30000002, 5.40000002, &
                5.50000002, 5.60000002, 5.70000003, 5.80000003, 5.90000003, &
                6.00000003, 6.10000003, 6.20000003, 6.30000003, 6.40000004, &
                6.50000004, 6.60000004, 6.70000004, 6.80000004, 6.90000004, &
                7.00000004, 7.10000005, 7.20000005, 7.30000005, 7.40000005, &
                7.50000005, 7.60000005, 7.70000006, 7.80000006, 7.90000006, &
                8.00000006 /

  data f_1354 / 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 1.09503647d-39, &
                5.47214550d-36, 2.42433983d-33, 2.75295034d-31, 1.21929718d-29, &
                2.48392125d-28, 2.33268145d-27, 8.68623633d-27, 1.00166284d-26, &
                3.63126633d-27, 7.45174807d-28, 1.38224064d-28, 2.69270994d-29, &
                5.53314977d-30, 1.15313092d-30, 2.34195788d-31, 4.48242942d-32, &
                7.94976380d-33 /


  data n_eis  / 60 /

  data t_eis1 / 1.99526231d+05, 2.23872114d+05, 2.51188643d+05, 2.81838293d+05, & 
                3.16227766d+05, 3.54813389d+05, 3.98107171d+05, 4.46683592d+05, &
                5.01187234d+05, 5.62341325d+05, 6.30957344d+05, 7.07945784d+05, &
                7.94328235d+05, 8.91250938d+05, 1.00000000d+06, 1.12201845d+06, &
                1.25892541d+06, 1.41253754d+06, 1.58489319d+06, 1.77827941d+06, &
                1.99526231d+06, 2.23872114d+06, 2.51188643d+06, 2.81838293d+06, &
                3.16227766d+06, 3.54813389d+06, 3.98107171d+06, 4.46683592d+06, &
                5.01187234d+06, 5.62341325d+06, 6.30957344d+06, 7.07945784d+06, &
                7.94328235d+06, 8.91250938d+06, 1.00000000d+07, 1.12201845d+07, &
                1.25892541d+07, 1.41253754d+07, 1.58489319d+07, 1.77827941d+07, &
                1.99526231d+07, 2.23872114d+07, 2.51188643d+07, 2.81838293d+07, &
                3.16227766d+07, 3.54813389d+07, 3.98107171d+07, 4.46683592d+07, &
                5.01187234d+07, 5.62341325d+07, 6.30957344d+07, 7.07945784d+07, &
                7.94328235d+07, 8.91250938d+07, 1.00000000d+08, 1.12201845d+08, &
                1.25892541d+08, 1.41253754d+08, 1.58489319d+08, 1.77827941d+08 /

  data t_eis2 / 1.99526231d+06, 2.23872114d+06, 2.51188643d+06, 2.81838293d+06, & 
                3.16227766d+06, 3.54813389d+06, 3.98107171d+06, 4.46683592d+06, &
                5.01187234d+06, 5.62341325d+06, 6.30957344d+06, 7.07945784d+06, &
                7.94328235d+06, 8.91250938d+06, 1.00000000d+07, 1.12201845d+07, &
                1.25892541d+07, 1.41253754d+07, 1.58489319d+07, 1.77827941d+07, &
                1.99526231d+07, 2.23872114d+07, 2.51188643d+07, 2.81838293d+07, &
                3.16227766d+07, 3.54813389d+07, 3.98107171d+07, 4.46683592d+07, &
                5.01187234d+07, 5.62341325d+07, 6.30957344d+07, 7.07945784d+07, &
                7.94328235d+07, 8.91250938d+07, 1.00000000d+08, 1.12201845d+08, &
                1.25892541d+08, 1.41253754d+08, 1.58489319d+08, 1.77827941d+08, &
                1.99526231d+08, 2.23872114d+08, 2.51188643d+08, 2.81838293d+08, &
                3.16227766d+08, 3.54813389d+08, 3.98107171d+08, 4.46683592d+08, &
                5.01187234d+08, 5.62341325d+08, 6.30957344d+08, 7.07945784d+08, &
                7.94328235d+08, 8.91250938d+08, 1.00000000d+09, 1.12201845d+09, &
                1.25892541d+09, 1.41253754d+09, 1.58489319d+09, 1.77827941d+09 /

  data f_263 /  0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, & 
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, &
                0.00000000d+00, 4.46454917d-45, 3.26774829d-42, 1.25292566d-39, &
                2.66922338d-37, 3.28497742d-35, 2.38677554d-33, 1.03937729d-31, &
                2.75075687d-30, 4.47961733d-29, 4.46729177d-28, 2.64862689d-27, &
                8.90863800d-27, 1.72437548d-26, 2.22217752d-26, 2.27999477d-26, &
                2.08264363d-26, 1.78226687d-26, 1.45821699d-26, 1.14675379d-26, &
                8.63082492d-27, 6.15925429d-27, 4.11252514d-27, 2.51530564d-27, &
                1.37090986d-27, 6.42443134d-28, 2.48392636d-28, 7.59187874d-29, &
                1.77852938d-29, 3.23945221d-30, 4.90533903d-31, 6.75458158d-32, &
                9.06878868d-33, 1.23927474d-33, 1.75769395d-34, 2.60710914d-35, &
                4.04318030d-36, 6.53500581d-37, 1.09365022d-37, 1.88383322d-38, &
                3.31425233d-39, 5.90964084d-40, 1.06147549d-40, 1.90706170d-41, &
                3.41331584d-42, 6.07310718d-43, 1.07364738d-43, 1.89085498d-44, &
                3.32598922d-45, 5.87125640d-46, 0.00000000d+00, 0.00000000d+00 /

  data f_264 /  0.00000000d+00, 2.81670057d-46, 1.28007268d-43, 2.54586603d-41, & 
                2.67887256d-39, 1.68413285d-37, 6.85702304d-36, 1.91797284d-34, &
                3.84675839d-33, 5.69939170d-32, 6.36224608d-31, 5.39176489d-30, &
                3.45478458d-29, 1.64848693d-28, 5.71476364d-28, 1.39909997d-27, &
                2.37743056d-27, 2.86712530d-27, 2.65206348d-27, 2.07175767d-27, &
                1.47866767d-27, 1.01087374d-27, 6.79605811d-28, 4.54746770d-28, &
                3.04351751d-28, 2.03639149d-28, 1.35940991d-28, 9.01451939d-29, &
                5.91289972d-29, 3.81821178d-29, 2.41434696d-29, 1.48871220d-29, &
                8.93362094d-30, 5.21097445d-30, 2.95964719d-30, 1.64278748d-30, &
                8.95571660d-31, 4.82096011d-31, 2.57390991d-31, 1.36821781d-31, &
                7.27136350d-32, 3.87019426d-32, 2.06883430d-32, 1.11228884d-32, &
                6.01883313d-33, 3.27790676d-33, 1.79805012d-33, 9.93085346d-34, &
                5.52139556d-34, 3.08881387d-34, 1.73890315d-34, 9.84434964d-35, &
                5.60603378d-35, 3.20626492d-35, 1.84111068d-35, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00 /

  data f_192 /  0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 4.35772105d-44, & 
                1.26162319d-41, 1.97471205d-39, 1.83409019d-37, 1.08206288d-35, &
                4.27914363d-34, 1.17943846d-32, 2.32565755d-31, 3.33087991d-30, &
                3.47013260d-29, 2.60375866d-28, 1.37737127d-27, 5.01053913d-27, &
                1.23479810d-26, 2.11310542d-26, 2.71831513d-26, 2.89851163d-26, &
                2.77312376d-26, 2.50025229d-26, 2.18323661d-26, 1.86980322d-26, &
                1.58035034d-26, 1.31985651d-26, 1.08733133d-26, 8.81804906d-27, &
                7.00417973d-27, 5.43356567d-27, 4.09857884d-27, 2.99651764d-27, &
                2.11902962d-27, 1.45014127d-27, 9.62291023d-28, 6.21548647d-28, &
                3.92807578d-28, 2.44230375d-28, 1.50167782d-28, 9.17611405d-29, &
                5.58707641d-29, 3.40570915d-29, 2.08030862d-29, 1.27588676d-29, &
                7.86535588d-30, 4.87646151d-30, 3.03888897d-30, 1.90578649d-30, &
                1.20195947d-30, 7.61955060d-31, 4.85602199d-31, 3.11049969d-31, &
                2.00087065d-31, 1.29223740d-31, 8.37422008d-32, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00 /

  data f_255 /  0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 1.76014287d-44, & 
                5.07057938d-42, 7.90473970d-40, 7.31852999d-38, 4.30709255d-36, &
                1.70009061d-34, 4.67925160d-33, 9.21703546d-32, 1.31918676d-30, &
                1.37393161d-29, 1.03102379d-28, 5.45694018d-28, 1.98699648d-27, &
                4.90346776d-27, 8.40524725d-27, 1.08321456d-26, 1.15714525d-26, &
                1.10905152d-26, 1.00155023d-26, 8.75799161d-27, 7.50935839d-27, &
                6.35253533d-27, 5.30919268d-27, 4.37669455d-27, 3.55185164d-27, &
                2.82347055d-27, 2.19257595d-27, 1.65589541d-27, 1.21224987d-27, &
                8.58395132d-28, 5.88163935d-28, 3.90721447d-28, 2.52593407d-28, &
                1.59739995d-28, 9.93802874d-29, 6.11343388d-29, 3.73711135d-29, &
                2.27618743d-29, 1.38793199d-29, 8.48060787d-30, 5.20305940d-30, &
                3.20867365d-30, 1.99011277d-30, 1.24064551d-30, 7.78310544d-31, &
                4.91013681d-31, 3.11338381d-31, 1.98451675d-31, 1.27135460d-31, &
                8.17917486d-32, 5.28280497d-32, 3.42357159d-32, 0.00000000d+00, &
                0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00 /

  abstract interface
    subroutine get_subr1(w,x,ixI^L,ixO^L,res)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,nw)
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(out):: res(ixI^S)
    end subroutine get_subr1

  end interface

  abstract interface
    subroutine get_2var_subr_te(ixI^L, ixO^L, w, val1, val2)
      use mod_global_parameters
      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, nw)
      double precision, intent(out):: val1(ixI^S), val2(ixI^S)
    end subroutine get_2var_subr_te
  end interface

  type te_fluid

    procedure (get_subr1), pointer, nopass :: get_rho => null()
    procedure (get_subr1), pointer, nopass :: get_pthermal => null()
    procedure (get_subr1), pointer, nopass :: get_var_Rfactor => null()
    procedure (get_2var_subr_te), pointer, nopass :: get_ne_nH => null()

  end type te_fluid

  type radsyn_euv_cache
    integer :: igrid=0
    integer :: level=0
    integer :: rft=1
    integer :: los_min=0
    integer :: los_max=0
    double precision, allocatable :: source(:^D&)
    double precision, allocatable :: opacity(:^D&)
    double precision, allocatable :: sourcev(:^D&)
    double precision, allocatable :: xface1(:),xface2(:),xface3(:)
    double precision, allocatable :: rface(:),thetaface(:),phiface(:)
    double precision, allocatable :: rface2(:),theta_cos(:),phi_sin(:),phi_cos(:)
    double precision :: box_min(1:3)=0.d0
    double precision :: box_max(1:3)=0.d0
    integer :: ixPmin1=1
    integer :: ixPmax1=0
    integer :: ixPmin2=1
    integer :: ixPmax2=0
    logical :: has_pixels=.false.
  end type radsyn_euv_cache

  character(len=std_len) :: ray_method_active='legacy'
  logical :: sph_use_dda=.false.


  contains

    subroutine check_synthetic_emission_options(datatype)
      use mod_global_parameters

      character(len=*), intent(in) :: datatype

      if (trim(radiation_transfer) /= 'thin' .and. trim(radiation_transfer) /= 'thick') then
        call mpistop("bad radiation_transfer")
      endif

      sph_use_dda=.false.
      select case(trim(ray_method))
      case('auto','')
        if (datatype=='image_euv' .and. coordinate==spherical) then
          ray_method_active='spherical'
          sph_use_dda=.true.
        else if (datatype=='image_euv' .and. coordinate==cartesian .and. dat_resolution .and. slab) then
          ray_method_active='cart'
        else
          ray_method_active='legacy'
        endif
      case('legacy')
        ray_method_active='legacy'
      case('cart','cart_dda')
        ray_method_active='cart'
      case('spherical','sph_intersection')
        ray_method_active='spherical'
      case('sph_dda','spherical_dda')
        ray_method_active='spherical'
        sph_use_dda=.true.
      case default
        call mpistop("bad ray_method")
      end select

      if (trim(emission_model) /= 'auto' .and. trim(emission_model) /= 'euv_aia' .and. &
          trim(emission_model) /= 'white_light' .and. trim(emission_model) /= 'radio_ff' .and. &
          trim(emission_model) /= 'pseudo_current') then
        call mpistop("bad emission_model")
      endif

      if ((output_tau .or. output_absorption_fraction) .and. trim(radiation_transfer) /= 'thick') then
        call mpistop("tau and absorption-fraction output need thick transfer")
      endif

      if (radsyn_pixel_batch<1) then
        call mpistop("radsyn_pixel_batch must be positive")
      endif
      if (radsyn_segment_batch_factor<0) then
        call mpistop("radsyn_segment_batch_factor must be non-negative")
      endif
      if (radsyn_segment_memory_mb<=zero) then
        call mpistop("radsyn_segment_memory_mb must be positive")
      endif
      if (radsyn_segment_comm_factor<1) then
        call mpistop("radsyn_segment_comm_factor must be positive")
      endif

      if (instrument_postprocess) then
        if (datatype /= 'image_euv' .or. .not. dat_resolution) then
          call mpistop("instrument_postprocess currently needs dat-resolution EUV images")
        endif
        if (trim(ray_method_active) == 'spherical') then
          call mpistop("instrument_postprocess is not yet supported for spherical rays")
        endif
        if (trim(emission_model) == 'pseudo_current') then
          call mpistop("instrument_postprocess currently supports only EUV AIA or radio_ff images")
        endif
        if (trim(emission_model) == 'radio_ff' .and. radio_beam_fwhm<=zero) then
          call mpistop("radio_ff instrument_postprocess needs radio_beam_fwhm > 0 arcsec")
        endif
      endif

      select case(trim(emission_model))
      case('auto')
        continue
      case('euv_aia')
        if (datatype /= 'image_euv' .and. datatype /= 'spectrum_euv') then
          call mpistop("emission_model=euv_aia is only valid for EUV synthesis")
        endif
      case('white_light')
        if (datatype /= 'image_whitelight') then
          call mpistop("emission_model=white_light is only valid for white-light synthesis")
        endif
      case('radio_ff')
        if (datatype /= 'image_euv') then
          call mpistop("emission_model=radio_ff currently reuses EUV-image convert types")
        endif
        if (radio_frequency<=zero) then
          call mpistop("emission_model=radio_ff needs radio_frequency > 0")
        endif
      case('pseudo_current')
        if (datatype /= 'image_euv') then
          call mpistop("emission_model=pseudo_current is only valid for EUV-image convert types")
        endif
        if (trim(radiation_transfer) /= 'thin') then
          call mpistop("emission_model=pseudo_current currently supports only thin transfer")
        endif
      end select

      if (trim(ray_method_active) == 'cart') then
        if (datatype /= 'image_euv' .or. coordinate /= cartesian .or. .not. slab) then
          call mpistop("ray_method=cart needs Cartesian EUV slab images")
        endif
      endif
      if (trim(ray_method_active) == 'spherical') then
        {^IFONED
        call mpistop("ray_method=spherical currently needs 3D spherical grids")
        }
        {^IFTWOD
        call mpistop("ray_method=spherical currently needs 3D spherical grids")
        }
        {^IFTHREED
        if (datatype /= 'image_euv' .or. coordinate /= spherical .or. &
            (trim(radiation_transfer) /= 'thin' .and. trim(radiation_transfer) /= 'thick')) then
          call mpistop("bad ray_method=spherical mode")
        endif
        if (trim(emission_model) /= 'auto' .and. trim(emission_model) /= 'euv_aia') then
          call mpistop("ray_method=spherical currently supports only EUV AIA emission")
        endif
        if (xprobmin2<=1.d-10 .or. xprobmax2>=dpi-1.d-10) then
          call mpistop("ray_method=spherical does not support polar-axis crossing domains")
        endif
        if (xprobmax3<=xprobmin3 .or. xprobmax3-xprobmin3>=2.d0*dpi-1.d-10) then
          call mpistop("ray_method=spherical does not support phi-wrapping domains")
        endif
        }
      endif
      if (trim(radiation_transfer) == 'thick') then
        if (datatype /= 'image_euv') then
          call mpistop("thick transfer is only defined for EUV images")
        endif
        if (trim(ray_method_active) == 'spherical') then
          continue
        else if (trim(ray_method_active) == 'cart') then
          if (.not. slab) call mpistop("cartesian thick EUV currently needs slab output")
        else if (.not. slab .or. .not. dat_resolution) then
          call mpistop("thick EUV currently needs Cartesian dat_resolution output")
        endif
        if (trim(ray_method_active) /= 'cart' .and. &
            trim(ray_method_active) /= 'spherical' .and. &
            .not. ((LOS_phi==0 .and. LOS_theta==90) .or. &
                   (LOS_phi==90 .and. LOS_theta==90) .or. LOS_theta==0)) then
          call mpistop("thick EUV currently needs x/y/z-aligned LOS")
        endif
      endif
    end subroutine check_synthetic_emission_options

    subroutine integrate_transfer_step_first_order(emissivity,opacity,path_length,intensity,tau)
      ! First-order formal-solution step used by the planned ordered LOS transfer.
      double precision, intent(in) :: emissivity,opacity,path_length
      double precision, intent(inout) :: intensity,tau

      double precision :: dtau

      if (path_length<=zero) return
      intensity=intensity+transfer_attenuation(tau)*max(zero,emissivity)*path_length
      dtau=max(zero,opacity)*path_length
      tau=tau+dtau
    end subroutine integrate_transfer_step_first_order

    logical function radsyn_euv_has_doppler_output()
      radsyn_euv_has_doppler_output=trim(emission_model)/='pseudo_current' .and. &
           trim(emission_model)/='radio_ff' .and. &
           .not. (coordinate==spherical .and. trim(ray_method_active)=='spherical')
    end function radsyn_euv_has_doppler_output

    integer function radsyn_euv_num_outputs(has_doppler,has_thick) result(num_outputs)
      logical, intent(in) :: has_doppler,has_thick

      num_outputs=1
      if (has_doppler) num_outputs=num_outputs+1
      if (has_thick .and. output_tau) num_outputs=num_outputs+1
      if (has_thick .and. output_absorption_fraction) num_outputs=num_outputs+1
    end function radsyn_euv_num_outputs

    subroutine normalize_euv_doppler(nI1,nI2,EUV,Dpl,unitv)
      integer, intent(in) :: nI1,nI2
      double precision, intent(in) :: EUV(nI1,nI2),unitv
      double precision, intent(inout) :: Dpl(nI1,nI2)

      integer :: ix1,ix2

      do ix1=1,nI1
        do ix2=1,nI2
          if (EUV(ix1,ix2)/=zero) then
            Dpl(ix1,ix2)=(Dpl(ix1,ix2)/EUV(ix1,ix2))*unitv
          else
            Dpl(ix1,ix2)=zero
          endif
          if (abs(Dpl(ix1,ix2))<smalldouble) Dpl(ix1,ix2)=zero
        enddo
      enddo
    end subroutine normalize_euv_doppler

    subroutine fill_euv_absorption_fraction(nI1,nI2,EUV,EUVthin,smallflux,Absorption,cap_to_one)
      integer, intent(in) :: nI1,nI2
      double precision, intent(in) :: EUV(nI1,nI2),EUVthin(nI1,nI2),smallflux
      double precision, intent(out) :: Absorption(nI1,nI2)
      logical, intent(in), optional :: cap_to_one

      integer :: ix1,ix2
      logical :: cap_absorption

      Absorption=zero
      cap_absorption=.false.
      if (present(cap_to_one)) cap_absorption=cap_to_one
      do ix1=1,nI1
        do ix2=1,nI2
          if (EUVthin(ix1,ix2)>smallflux) then
            Absorption(ix1,ix2)=max(zero,(EUVthin(ix1,ix2)-EUV(ix1,ix2))/EUVthin(ix1,ix2))
            if (cap_absorption) Absorption(ix1,ix2)=min(one,Absorption(ix1,ix2))
          endif
        enddo
      enddo
    end subroutine fill_euv_absorption_fraction

    subroutine pack_euv_image_outputs(nI1,nI2,EUV,wI,smallflux,has_doppler,has_thick,Dpl,Tau,EUVthin,&
                                      cap_absorption)
      integer, intent(in) :: nI1,nI2
      double precision, intent(in) :: EUV(nI1,nI2),smallflux
      double precision, intent(inout) :: wI(:,:,:)
      logical, intent(in) :: has_doppler,has_thick
      double precision, intent(in), optional :: Dpl(nI1,nI2),Tau(nI1,nI2),EUVthin(nI1,nI2)
      logical, intent(in), optional :: cap_absorption

      integer :: iw
      double precision, allocatable :: Absorption(:,:)

      wI=zero
      wI(:,:,1)=EUV(:,:)
      iw=1
      if (has_doppler) then
        if (.not. present(Dpl)) call mpistop("Doppler output requested without Doppler image")
        iw=iw+1
        wI(:,:,iw)=Dpl(:,:)
      endif
      if (has_thick .and. output_tau) then
        if (.not. present(Tau)) call mpistop("tau output requested without tau image")
        iw=iw+1
        wI(:,:,iw)=Tau(:,:)
      endif
      if (has_thick .and. output_absorption_fraction) then
        if (.not. present(EUVthin)) call mpistop("absorption output requested without thin image")
        allocate(Absorption(nI1,nI2))
        call fill_euv_absorption_fraction(nI1,nI2,EUV,EUVthin,smallflux,Absorption,cap_absorption)
        iw=iw+1
        wI(:,:,iw)=Absorption(:,:)
        deallocate(Absorption)
      endif
    end subroutine pack_euv_image_outputs

    subroutine radsyn_get_segment_batch_limits(pixel_batch_target,segment_batch_target,segment_comm_target)
      integer, intent(out) :: pixel_batch_target,segment_batch_target,segment_comm_target

      pixel_batch_target=max(1,radsyn_pixel_batch)
      if (radsyn_segment_batch_factor>0) then
        segment_batch_target=max(128,radsyn_segment_batch_factor*pixel_batch_target)
      else
        segment_batch_target=max(128,int(min(dble(huge(segment_batch_target)),&
             max(128.d0,radsyn_segment_memory_mb*1048576.d0/256.d0))))
      endif
      segment_comm_target=max(128,radsyn_segment_comm_factor*pixel_batch_target)
    end subroutine radsyn_get_segment_batch_limits

    double precision function transfer_attenuation(tau)
      double precision, intent(in) :: tau

      if (tau<=zero) then
        transfer_attenuation=one
      else
        transfer_attenuation=exp_clamped(-tau)
      endif
    end function transfer_attenuation

    double precision function exp_clamped(argument)
      double precision, intent(in) :: argument

      if (argument<-700.d0) then
        exp_clamped=zero
      else if (argument>700.d0) then
        exp_clamped=huge(one)
      else
        exp_clamped=exp(argument)
      endif
    end function exp_clamped

    double precision function pow10_clamped(exponent)
      double precision, intent(in) :: exponent

      if (exponent>300.d0) then
        pow10_clamped=1.d300
      else if (exponent<-300.d0) then
        pow10_clamped=zero
      else
        pow10_clamped=10.d0**exponent
      endif
    end function pow10_clamped

    double precision function interpolate_response_value(temperature,t_table,f_table,n_table,log_temperature,log_response)
      double precision, intent(in) :: temperature
      integer, intent(in) :: n_table
      double precision, intent(in) :: t_table(n_table),f_table(n_table)
      logical, intent(in) :: log_temperature,log_response

      integer :: ilo,ihi,imid
      double precision :: temp_lookup,response_lookup,flo,fhi

      interpolate_response_value=zero
      if (temperature<=zero) return
      if (log_temperature) then
        temp_lookup=log10(temperature)
      else
        temp_lookup=temperature
      endif
      if (temp_lookup<t_table(1) .or. temp_lookup>t_table(n_table)) return
      if (temp_lookup==t_table(n_table)) then
        if (log_response) then
          response_lookup=log10(max(f_table(n_table),1.d-99))
        else
          response_lookup=f_table(n_table)
        endif
      else
        ilo=1
        ihi=n_table
        do while (ihi-ilo>1)
          imid=(ilo+ihi)/2
          if (temp_lookup>=t_table(imid)) then
            ilo=imid
          else
            ihi=imid
          endif
        enddo
        if (log_response) then
          flo=log10(max(f_table(ilo),1.d-99))
          fhi=log10(max(f_table(ilo+1),1.d-99))
        else
          flo=f_table(ilo)
          fhi=f_table(ilo+1)
        endif
        response_lookup=flo*(temp_lookup-t_table(ilo+1))/(t_table(ilo)-t_table(ilo+1))+&
                        fhi*(temp_lookup-t_table(ilo))/(t_table(ilo+1)-t_table(ilo))
      endif

      if (log_response) then
        if (response_lookup>-99.d0) interpolate_response_value=10.d0**response_lookup
      else
        interpolate_response_value=response_lookup
      endif
      if (interpolate_response_value<zero) interpolate_response_value=zero
    end function interpolate_response_value

    subroutine apply_temperature_response(ixI^L,ixO^L,Te,flux,t_table,f_table,n_table,log_temperature,log_response)
      integer, intent(in) :: ixI^L, ixO^L, n_table
      double precision, intent(in) :: Te(ixI^S),t_table(n_table),f_table(n_table)
      double precision, intent(inout) :: flux(ixI^S)
      logical, intent(in) :: log_temperature,log_response

      integer :: ix^D
      double precision :: GT

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        GT=interpolate_response_value(Te(ix^D),t_table,f_table,n_table,log_temperature,log_response)
        flux(ix^D)=flux(ix^D)*GT
        if (flux(ix^D)<zero) flux(ix^D)=zero
      {enddo\}
    end subroutine apply_temperature_response

    subroutine get_EUV_HHe_opacity(wl,ixI^L,ixO^L,w,x,fl,kappa)
      ! H I + He I + He II photoionization opacity in cm^-1.
      use mod_constants, only: kB_cgs
      use mod_eos, only: eos

      integer, intent(in) :: wl
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: kappa(ixI^S)

      integer :: ix^D
      double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)
      double precision :: wave_ratio,s_H1,s_He1,s_He2,Pe
      double precision :: log_H21,log_He21,log_He32,log_He321,logScaleHe,w_H21
      double precision :: term0,term1,term2,denHe,i0,j1,j2,be
      double precision :: N_H,N_H1,N_He1,N_He2
      double precision, parameter :: Xe_H21=13.6d0, Xe_He21=24.587d0, Xe_He32=54.416d0
      double precision, parameter :: rHe=0.1d0
      double precision, parameter :: sigma_H1=5.16d-20, sigma_He1=9.25d-19, sigma_He2=7.17d-19

      call fl%get_pthermal(w,x,ixI^L,ixO^L,pth)
      call fl%get_rho(w,x,ixI^L,ixO^L,Ne)
      call fl%get_var_Rfactor(w,x,ixI^L,ixO^L,Te)
      Te(ixO^S)=pth(ixO^S)/(Ne(ixO^S)*Te(ixO^S))*unit_temperature
      block
        double precision :: nH_dummy(ixI^S)
        call eos%get_ne_nH(ixI^L, ixO^L, w, Ne, nH_dummy)
      end block
      if (SI_unit) then
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity/1.d6
      else
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity
      endif

      wave_ratio=dble(wl)/171.d0
      s_H1=zero
      s_He1=zero
      s_He2=zero
      if (wl<=912) s_H1=wave_ratio**3*sigma_H1
      if (wl<=504) s_He1=wave_ratio**2*sigma_He1
      if (wl<=228) s_He2=wave_ratio**2.75d0*sigma_He2
      kappa(ixO^S)=zero

      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if (Te(ix^D)>zero .and. Ne(ix^D)>zero) then
          Pe=Ne(ix^D)*kB_cgs*Te(ix^D)
          if (Pe>zero) then
            log_H21=2.5d0*log10(Te(ix^D))-5040.d0*Xe_H21/Te(ix^D)-log10(Pe)-0.48d0
            log_He21=log10(4.d0)+2.5d0*log10(Te(ix^D))-5040.d0*Xe_He21/Te(ix^D)-log10(Pe)-0.48d0
            log_He32=2.5d0*log10(Te(ix^D))-5040.d0*Xe_He32/Te(ix^D)-log10(Pe)-0.48d0
            w_H21=pow10_clamped(log_H21)
            i0=w_H21/(1.d0+w_H21)
            log_He321=log_He21+log_He32
            logScaleHe=max(zero,log_He21,log_He321)
            term0=pow10_clamped(-logScaleHe)
            term1=pow10_clamped(log_He21-logScaleHe)
            term2=pow10_clamped(log_He321-logScaleHe)
            denHe=term0+term1+term2
            if (denHe>zero) then
              j1=term1/denHe
              j2=term2/denHe
            else
              j1=zero
              j2=zero
            endif
            be=i0+rHe*(j1+2.d0*j2)
            if (be>smalldouble) then
              N_H=Ne(ix^D)/be
              N_H1=N_H*(1.d0-i0)
              N_He1=(1.d0-j1-j2)*rHe*N_H
              N_He2=j1*rHe*N_H
              kappa(ix^D)=max(zero,N_H1*s_H1+N_He1*s_He1+N_He2*s_He2)
            endif
          endif
        endif
      {enddo\}
    end subroutine get_EUV_HHe_opacity

    subroutine get_pseudo_current(igrid,ixI^L,ixO^L,w,source)
      integer, intent(in) :: igrid
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(out) :: source(ixI^S)

      integer :: ix^D,idir,idirmin,idirmin0
      double precision :: current(ixI^S,7-2*ndir:3)

      if (.not. allocated(iw_mag)) then
        call mpistop("emission_model=pseudo_current needs magnetic-field variables")
      endif

      idirmin0=7-2*ndir
      current=zero
      call curlvector(w(ixI^S,iw_mag(1:ndir)),ixI^L,ixO^L,current,idirmin,idirmin0,ndir)
      if (B0field) then
        current(ixO^S,idirmin0:3)=current(ixO^S,idirmin0:3)+ps(igrid)%J0(ixO^S,idirmin0:3)
      endif

      source(ixI^S)=zero
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        do idir=idirmin0,3
          source(ix^D)=source(ix^D)+current(ix^D,idir)**2
        enddo
      {enddo\}
    end subroutine get_pseudo_current

    subroutine get_radio_ff_source_opacity(ixI^L,ixO^L,w,x,fl,source,kappa)
      use mod_eos, only: eos

      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: source(ixI^S),kappa(ixI^S)

      integer :: ix^D
      double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)
      double precision :: nH_dummy(ixI^S),gff

      call fl%get_pthermal(w,x,ixI^L,ixO^L,pth)
      call fl%get_rho(w,x,ixI^L,ixO^L,Ne)
      call fl%get_var_Rfactor(w,x,ixI^L,ixO^L,Te)
      Te(ixO^S)=pth(ixO^S)/(Ne(ixO^S)*Te(ixO^S))*unit_temperature
      call eos%get_ne_nH(ixI^L,ixO^L,w,Ne,nH_dummy)
      if (SI_unit) then
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity/1.d6
      else
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity
      endif

      source(ixI^S)=zero
      kappa(ixI^S)=zero
      {do ix^DB=ixOmin^DB,ixOmax^DB\}
        if (Te(ix^D)>zero .and. Ne(ix^D)>zero) then
          if (Te(ix^D)<2.d5) then
            gff=18.2d0+1.5d0*log(Te(ix^D))-log(radio_frequency)
          else
            gff=24.5d0+log(Te(ix^D))-log(radio_frequency)
          endif
          gff=max(one,gff)
          kappa(ix^D)=9.78d-3*Ne(ix^D)**2*gff/(radio_frequency**2*Te(ix^D)**1.5d0)
          source(ix^D)=Te(ix^D)*kappa(ix^D)
        endif
      {enddo\}
    end subroutine get_radio_ff_source_opacity

    subroutine get_line_info(wl,ion,mass,logTe,line_center,spatial_px,spectral_px,sigma_PSF,width_slit)
      ! get information of the spectral line
      ! wl: wavelength
      ! mass: ion mass, unit -- proton mass
      ! logTe: peak temperature of emission line in logarithm
      ! line_center: center wavelength of emission line, unit -- Angstrom (0.1 nm) 
      ! spatial_px: pixel size in space of instrument (for image), unit -- arcsec
      ! spectral_px: pixel size in wagelength of instrument (for spectrum), unit -- Angstrom
      ! sigma_PSF: width of point spread function core (for instrument), unit -- pixel
      ! width_slit: width of slit for spectrograph, unit -- arcsec
      use mod_global_parameters

      integer, intent(in) :: wl
      integer, intent(out) :: mass
      character(len=30), intent(out) :: ion
      double precision, intent(out) :: logTe,line_center,spatial_px,spectral_px
      double precision, intent(out) :: sigma_PSF,width_slit

      select case(wl)
      case(304)
        ion='He II'
        mass=4
        logTe=4.7d0
        line_center=303.8d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.895d0
        width_slit=0.6d0
      case(171)
        ion='Fe IX'
        mass=56
        logTe=5.8d0
        line_center=171.1d0
        spatial_px=0.6d0
        spectral_px=0.02d0 
        sigma_PSF=1.019d0
        width_slit=0.6d0
      case(193)
        ion='Fe XXIV'
        mass=56
        logTe=7.3d0
        line_center=193.5d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.813d0
        width_slit=0.6d0
      case(211)
        ion='Fe XIV'
        mass=56
        logTe=6.3d0
        line_center=211.3d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.913d0
        width_slit=0.6d0
      case(335)
        ion='Fe XVI'
        mass=56
        logTe=6.4d0
        line_center=335.4d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=1.019d0
        width_slit=0.6d0
      case(94)
        ion='Fe XVIII'
        mass=56
        logTe=6.8d0
        line_center=93.9d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=1.025d0
        width_slit=0.6d0
      case(131)
        ion='Fe XXI'
        mass=56
        logTe=7.0d0
        line_center=131.0d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.984d0
        width_slit=0.6d0
      case(1354)
        ion='Fe XXI'
        mass=56
        logTe=7.0d0
        line_center=1354.1d0
        spatial_px=0.1663d0
        spectral_px=12.98d-3
        sigma_PSF=1.d0
        width_slit=0.33d0
      case(263)
        ion='Fe XVI'
        mass=56
        logTe=6.4d0
        line_center=262.976d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case(264)
        ion='Fe XXIII'
        mass=56
        logTe=7.1d0
        line_center=263.765d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case(192)
        ion='Fe XXIV'
        mass=56
        logTe=7.2d0
        line_center=192.028d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case(255)
        ion='Fe XXIV'
        mass=56
        logTe=7.2d0
        line_center=255.113d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case default
        call mpistop("No information about this line")
      end select

      spatial_px=spatial_px/instrument_resolution_factor
    end subroutine get_line_info
    
    subroutine get_EUV(wl,ixI^L,ixO^L,w,x,fl,flux)
      ! calculate the local emission intensity of given EUV line (optically thin)
      ! wavelength is the wave length of the emission line
      ! unit [DN cm^-1 s^-1 pixel^-1]
      ! ingrate flux along line of sight: DN s^-1 pixel^-1
      use mod_global_parameters
      use mod_eos, only: eos

      integer, intent(in) :: wl
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: flux(ixI^S)

      integer :: ix^D
      double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)

      call fl%get_pthermal(w,x,ixI^L,ixO^L,pth)
      call fl%get_rho(w,x,ixI^L,ixO^L,Ne)
      call fl%get_var_Rfactor(w,x,ixI^L,ixO^L,Te)
      Te(ixO^S)=pth(ixO^S)/(Ne(ixO^S)*Te(ixO^S))*unit_temperature
      ! get actual electron density from EoS (replaces rho with ne)
      block
        double precision :: nH_dummy(ixI^S)
        call eos%get_ne_nH(ixI^L, ixO^L, w, Ne, nH_dummy)
      end block
      if (SI_unit) then
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity/1.d6 ! m^-3 -> cm-3
        flux(ixO^S)=Ne(ixO^S)**2
      else
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity
        flux(ixO^S)=Ne(ixO^S)**2
      endif

      select case(wl)
      case(94)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_aia,f_94,n_aia,.true.,.true.)
      case(131)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_aia,f_131,n_aia,.true.,.true.)
      case(171)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_aia,f_171,n_aia,.true.,.true.)
      case(193)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_aia,f_193,n_aia,.true.,.true.)
      case(211)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_aia,f_211,n_aia,.true.,.true.)
      case(304)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_aia,f_304,n_aia,.true.,.true.)
      case(335)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_aia,f_335,n_aia,.true.,.true.)
      case(1354)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_iris,f_1354,n_iris,.true.,.true.)
      case(263)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_eis1,f_263,n_eis,.false.,.false.)
      case(264)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_eis2,f_264,n_eis,.false.,.false.)
      case(192)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_eis2,f_192,n_eis,.false.,.false.)
      case(255)
        call apply_temperature_response(ixI^L,ixO^L,Te,flux,t_eis2,f_255,n_eis,.false.,.false.)
      case default
        call mpistop("Unknown wavelength")
      end select
    end subroutine get_EUV

    subroutine get_SXR(ixI^L,ixO^L,w,x,fl,flux,El,Eu)
      !synthesize thermal SXR from El keV to Eu keV released by cm^-3/m^-3
      ! volume of plasma in 1 s
      !flux (cgs): photons cm^-3 s^-1
      !flux (SI): photons m^-3 s^-1
      use mod_global_parameters
      use mod_eos, only: eos

      integer, intent(in)           :: ixI^L,ixO^L
      integer, intent(in)           :: El,Eu
      double precision, intent(in)  :: x(ixI^S,1:ndim)
      double precision, intent(in)  :: w(ixI^S,nw)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: flux(ixI^S)

      integer :: ix^D,ixO^D
      integer :: iE,numE
      double precision :: I0,kb,keV,dE,Ei
      double precision :: pth(ixI^S),Te(ixI^S),kbT(ixI^S)
      double precision :: Ne(ixI^S),gff(ixI^S),fi(ixI^S)
      double precision :: EM(ixI^S)

      I0=3.01d-15   ! I0*4*pi*AU**2, I0 from Pinto (2015)
      kb=const_kb
      keV=1.0d3*const_ev
      dE=0.1
      numE=floor((Eu-El)/dE)
      call fl%get_pthermal(w,x,ixI^L,ixO^L,pth)
      call fl%get_rho(w,x,ixI^L,ixO^L,Ne)
      call fl%get_var_Rfactor(w,x,ixI^L,ixO^L,Te)
      Te(ixO^S)=pth(ixO^S)/(Ne(ixO^S)*Te(ixO^S))*unit_temperature
      ! get actual electron density from EoS (replaces rho with ne)
      block
        double precision :: nH_dummy(ixI^S)
        call eos%get_ne_nH(ixI^L, ixO^L, w, Ne, nH_dummy)
      end block
      if (SI_unit) then
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity/1.d6 ! m^-3 -> cm-3
        EM(ixO^S)=(Ne(ixO^S))**2*1.d6 ! cm^-3 m^-3
      else
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity
        EM(ixO^S)=(Ne(ixO^S))**2
      endif
      kbT(ixO^S)=kb*Te(ixO^S)/keV
      flux(ixO^S)=0.0d0
      do iE=0,numE-1
        Ei=dE*iE+El*1.d0
        gff(ixO^S)=1.d0
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if (kbT(ix^D)>0.01*Ei) then
            if(kbT(ix^D)<Ei) gff(ix^D)=(kbT(ix^D)/Ei)**0.4
            fi(ix^D)=(EM(ix^D)*gff(ix^D))*exp_clamped(-Ei/(kbT(ix^D)))/(Ei*dsqrt(kbT(ix^D)))
          else
            fi(ix^D)=zero
          endif
        {enddo\}
        flux(ixO^S)=flux(ixO^S)+fi(ixO^S)*dE
      enddo
      flux(ixO^S)=flux(ixO^S)*I0
    end subroutine get_SXR

    subroutine get_GOES_SXR_flux(xbox^L,fl,eflux)
      !get GOES SXR 1-8A flux observing at 1AU from given box [w/m^2]
      use mod_global_parameters

      double precision, intent(in) :: xbox^L
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: eflux

      double precision :: dxb^D,xb^L
      integer :: iigrid,igrid,j
      integer :: ixO^L,ixI^L,ix^D
      double precision :: eflux_grid,eflux_pe

      ^D&ixImin^D=ixglo^D;
      ^D&ixImax^D=ixghi^D;
      ^D&ixOmin^D=ixmlo^D;
      ^D&ixOmax^D=ixmhi^D;
      eflux_pe=zero
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
        ^D&xbmin^D=rnode(rpxmin^D_,igrid);
        ^D&xbmax^D=rnode(rpxmax^D_,igrid);
        call get_GOES_flux_grid(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,ps(igrid)%dvolume(ixI^S),xbox^L,xb^L,fl,eflux_grid)
        eflux_pe=eflux_pe+eflux_grid
      enddo
      call MPI_ALLREDUCE(eflux_pe,eflux,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    end subroutine get_GOES_SXR_flux

    subroutine get_GOES_flux_grid(ixI^L,ixO^L,w,x,dV,xbox^L,xb^L,fl,eflux_grid)
      use mod_global_parameters
      use mod_eos, only: eos

      integer, intent(in)           :: ixI^L,ixO^L
      double precision, intent(in)  :: x(ixI^S,1:ndim),dV(ixI^S)
      double precision, intent(in)  :: w(ixI^S,nw)
      double precision, intent(in)  :: xbox^L,xb^L
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: eflux_grid

      integer :: ix^D,ixO^D,ixb^L
      integer :: iE,numE,j,inbox
      double precision :: I0,kb,keV,dE,Ei,El,Eu,A_cgs
      double precision :: pth(ixI^S),Te(ixI^S),kbT(ixI^S)
      double precision :: Ne(ixI^S),EM(ixI^S)
      double precision :: gff,fi,erg_SI

      ! check whether the grid is inside given box
      inbox=0
      {if (xbmin^D<xboxmax^D .and. xbmax^D>xboxmin^D) inbox=inbox+1\}

      if (inbox==ndim) then
        ! indexes for cells inside given box
        ^D&ixbmin^D=ixOmin^D;
        ^D&ixbmax^D=ixOmax^D;
        {if (xbmax^D>xboxmax^D) ixbmax^D=ixOmax^D-ceiling((xbmax^D-xboxmax^D)/dxlevel(^D))\}
        {if (xbmin^D<xboxmin^D) ixbmin^D=ceiling((xboxmin^D-xbmin^D)/dxlevel(^D))+ixOmin^D\}

        I0=1.07d-38 ! photon flux index for observed at 1AU [photon cm^3 m^-2 s^-1 keV^-1]
        kb=const_kb
        keV=1.0d3*const_ev
        erg_SI=1.d-7
        A_cgs=1.d-8 ! Angstrom
        El=const_h*const_c/(8.d0*A_cgs)/keV ! 8 A
        Eu=const_h*const_c/(1.d0*A_cgs)/keV ! 1 A
        dE=0.1  ! keV
        numE=floor((Eu-El)/dE)
        call fl%get_pthermal(w,x,ixI^L,ixb^L,pth)
        call fl%get_rho(w,x,ixI^L,ixb^L,Ne)
        call fl%get_var_Rfactor(w,x,ixI^L,ixb^L,Te)
        Te(ixb^S)=pth(ixb^S)/(Ne(ixb^S)*Te(ixb^S))*unit_temperature
        ! get actual electron density from EoS (replaces rho with ne)
        block
          double precision :: nH_dummy(ixI^S)
          call eos%get_ne_nH(ixI^L, ixb^L, w, Ne, nH_dummy)
        end block
        if (SI_unit) then
          Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity/1.d6 ! m^-3 -> cm-3
          EM(ixb^S)=(I0*(Ne(ixb^S))**2)*dV(ixb^S)*(unit_length*1.d2)**3 ! cm^-3
        else
          Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity
          EM(ixb^S)=(I0*(Ne(ixb^S))**2)*dV(ixb^S)*unit_length**3
        endif
        kbT(ixb^S)=kb*Te(ixb^S)/keV
        eflux_grid=0.0d0

        do iE=0,numE-1
          Ei=dE*iE+El
          {do ix^DB=ixbmin^DB,ixbmax^DB\}
            if (kbT(ix^D)>1.d-2*Ei) then
              if(kbT(ix^D)<Ei) then
                gff=(kbT(ix^D)/Ei)**0.4
              else
                gff=1.d0
              endif
              fi=(EM(ix^D)*gff)*exp_clamped(-Ei/(kbT(ix^D)))/(Ei*dsqrt(kbT(ix^D)))
              eflux_grid=eflux_grid+fi*dE*Ei
            endif
          {enddo\}
        enddo
        eflux_grid=eflux_grid*keV*erg_SI
      endif

    end subroutine get_GOES_flux_grid

  {^IFTHREED
    subroutine get_EUV_spectrum(qunit,fl)
      use mod_global_parameters

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20) :: datatype

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl,wslit
      double precision :: xslit,arcsec

      datatype='spectrum_euv'
      call check_synthetic_emission_options(datatype)
      arcsec=7.25d7/unit_length
      call get_line_info(spectrum_wl,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)

      if (mype==0) print *, '###################################################'
      select case(spectrum_wl)
      case (1354)
        if (mype==0) print *, 'Systhesizing EUV spectrum (observed by IRIS).'
      case (263,264,192,255)
        if (mype==0) print *, 'Systhesizing EUV spectrum (observed by Hinode/EIS).'
      case default
        call MPISTOP('Wrong wavelength!')
      end select

      if (spectrum_window_max<=spectrum_window_min) then
        call MPISTOP('Wrong spectrum window!')
      endif

      if (mype==0) write(*,'(a,f8.3,a)') ' Wavelength: ',lineCent,' Angstrom'
      if (mype==0) print *, 'Unit of EUV flux: DN s^-1 pixel^-1'

      if (dat_resolution) then
        if (mype==0) then
          write(*,'(a,f5.3,a,f5.1,a)') ' Supposed pixel: ',wlRsl,' Angstrom x ',spaceRsl*725.0, ' km'
          print *, 'Unit of wavelength: Angstrom (0.1 nm) '
          if (SI_unit) then
            write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
          else
            write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
          endif
          write(*,'(a,f8.1,a)') ' Supposed width of slit: ',wslit*725.0,' km'
        endif
        call get_spectrum_datresol(qunit,datatype,fl)
      else
        if (mype==0) then
          print *, 'Unit of wavelength: Angstrom (0.1 nm) '
          if (activate_unit_arcsec) then
            write(*,'(a,f5.3,a,f5.1,a)') ' Pixel: ',wlRsl,' Angstrom x ',spaceRsl*725.0, ' km'
            print *, 'Unit of length: arcsec (~725 km)'
            write(*,'(a,f8.1,a)') ' Location of slit: xI1 = ',location_slit,' arcsec'
            write(*,'(a,f8.1,a)') ' Width of slit: ',wslit,' arcsec'
          else
            if (SI_unit) then
              if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
            else
              if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
            endif
            write(*,'(a,f8.1,a)') ' Location of slit: xI1 = ',location_slit,' Unit_length'
            write(*,'(a,f8.1,a)') ' Width of slit: ',wslit*725.d0,' km'
          endif
        endif
        if (mype==0) print *, 'Direction of the slit: parallel to xI2 vector'
        if (coordinate==cartesian .or. coordinate==spherical) then
          call get_spectrum(qunit,datatype,fl)
        else
          call MPISTOP("EUV spectrum synthesis: support for sperical coordinates is to be added!")
        endif
      endif

      if (mype==0) print *, '###################################################'

    end subroutine get_EUV_spectrum

    subroutine get_spectrum_datresol(qunit,datatype,fl)

      integer, intent(in) :: qunit
      character(20), intent(in) :: datatype
      type(te_fluid), intent(in) :: fl

      integer :: numWL,numXS,iwL,ixS,numWI,numS
      double precision :: dwLg,xSmin,xSmax,wLmin,wLmax
      double precision, allocatable :: wL(:),xS(:),dwL(:),dxS(:)
      double precision, allocatable :: wI(:,:,:),spectra(:,:),spectra_rc(:,:)
      integer :: strtype,nstrb,nbb,nuni,nstr,bnx
      double precision :: qs,dxfirst,dxmid,lenstr

      integer :: iigrid,igrid,j,dir_loc
      double precision :: xbmin(1:ndim),xbmax(1:ndim)

      dwLg=1.d-3
      numWL=4*int((spectrum_window_max-spectrum_window_min)/(4.d0*dwLg))
      wLmin=(spectrum_window_max+spectrum_window_min)/2.d0-dwLg*numWL/2
      wLmax=(spectrum_window_max+spectrum_window_min)/2.d0+dwLg*numWL/2
      allocate(wL(numWL),dwL(numWL))
      dwL(:)=dwLg
      do iwL=1,numWL
        wL(iwL)=wLmin+iwL*dwLg-half*dwLg
      enddo

      select case(direction_slit)
      case (1)
        numXS=domain_nx1*2**(refine_max_level-1)
        xSmin=xprobmin1
        xSmax=xprobmax1
        bnx=block_nx1
        nbb=domain_nx1
        strtype=stretch_type(1)
        nstrb=nstretchedblocks_baselevel(1)
        qs=qstretch_baselevel(1)
        if (mype==0) print *, 'Direction of the slit: x'
      case (2)
        numXS=domain_nx2*2**(refine_max_level-1)
        xSmin=xprobmin2
        xSmax=xprobmax2
        bnx=block_nx2
        nbb=domain_nx2
        strtype=stretch_type(2)
        nstrb=nstretchedblocks_baselevel(2)
        qs=qstretch_baselevel(2)
        if (mype==0) print *, 'Direction of the slit: y'
      case (3)
        numXS=domain_nx3*2**(refine_max_level-1)
        xSmin=xprobmin3
        xSmax=xprobmax3
        bnx=block_nx3
        nbb=domain_nx3
        strtype=stretch_type(3)
        nstrb=nstretchedblocks_baselevel(3)
        qs=qstretch_baselevel(3)
        if (mype==0) print *, 'Direction of the slit: z'
      case default
        call MPISTOP('Wrong direction_slit')
      end select

      allocate(xS(numXS),dxS(numXS),spectra(numWL,numXS),spectra_rc(numWL,numXS))
      numWI=1
      allocate(wI(numWL,numXS,numWI))

      select case(strtype)
      case(0) ! uniform
        dxS(:)=(xSmax-xSmin)/numXS
        do ixS=1,numXS
          xS(ixS)=xSmin+dxS(ixS)*(ixS-half)
        enddo
      case(1) ! uni stretch
        qs=qs**(one/2**(refine_max_level-1))
        dxfirst=(xSmax-xSmin)*(one-qs)/(one-qs**numXS)
        dxS(1)=dxfirst
        do ixS=2,numXS
          dxS(ixS)=dxfirst*qs**(ixS-1)
          xS(ixS)=dxS(1)/(one-qs)*(one-qs**(ixS-1))+half*dxS(ixS)
        enddo
      case(2) ! symm stretch
        ! base level, nbb = nstr + nuni + nstr
        nstr=nstrb*bnx/2
        nuni=nbb-nstrb*bnx
        lenstr=(xSmax-xSmin)/(2.d0+nuni*(one-qs)/(one-qs**nstr))
        dxfirst=(xSmax-xSmin)/(dble(nuni)+2.d0/(one-qs)*(one-qs**nstr))
        dxmid=dxfirst
        ! refine_max level, numXI = nstr + nuni + nstr
        nstr=nstr*2**(refine_max_level-1)
        nuni=nuni*2**(refine_max_level-1)
        qs=qs**(one/2**(refine_max_level-1))
        dxfirst=lenstr*(one-qs)/(one-qs**nstr)
        dxmid=dxmid/2**(refine_max_level-1)
        ! uniform center
        if(nuni .gt. 0) then
          do ixS=nstr+1,nstr+nuni
            dxS(ixS)=dxmid
            xS(ixS)=lenstr+(dble(ixS)-0.5d0-nstr)*dxS(ixS)+xSmin
          enddo
        endif
        ! left half
        do ixS=nstr,1,-1
          dxS(ixS)=dxfirst*qs**(nstr-ixS)
          xS(ixS)=xSmin+lenstr-dxS(ixS)*half-dxfirst*(one-qs**(nstr-ixS))/(one-qs)
        enddo
        ! right half
        do ixS=nstr+nuni+1,numXS
          dxS(ixS)=dxfirst*qs**(ixS-nstr-nuni-1)
          xS(ixS)=xSmax-lenstr+dxS(ixS)*half+dxfirst*(one-qs**(ixS-nstr-nuni-1))/(one-qs)
        enddo
      case default
        call mpistop("unknown stretch type")
      end select

      if (LOS_phi==0 .and. LOS_theta==90 .and. direction_slit==2) then
      ! LOS->x slit->y
        dir_loc=3
      else if (LOS_phi==0 .and. LOS_theta==90 .and. direction_slit==3) then
      ! LOS->x slit->z
        dir_loc=2
      else if (LOS_phi==90 .and. LOS_theta==90 .and. direction_slit==1) then
      ! LOS->y slit->x
        dir_loc=3
      else if (LOS_phi==90 .and. LOS_theta==90 .and. direction_slit==3) then
      ! LOS->y slit->z
        dir_loc=1
      else if (LOS_theta==0 .and. direction_slit==1) then
      ! LOS->z slit->x
        dir_loc=2
      else if (LOS_theta==0 .and. direction_slit==2) then
      ! LOS->z slit->y
        dir_loc=1
      else
        call MPISTOP('Wrong combination of LOS and slit direction!')
      endif

      if (dir_loc==1) then
        if (location_slit>xprobmax1 .or. location_slit<xprobmin1) then
          call MPISTOP('Wrong value for location_slit!')
        endif
        if(mype==0) write(*,'(a,f8.1,a)') ' Location of slit: x = ',location_slit,' Unit_length'
      else if (dir_loc==2) then
        if (location_slit>xprobmax2 .or. location_slit<xprobmin2) then
          call MPISTOP('Wrong value for location_slit!')
        endif
        if(mype==0) write(*,'(a,f8.1,a)') ' Location of slit: y = ',location_slit,' Unit_length'
      else
        if (location_slit>xprobmax3 .or. location_slit<xprobmin3) then
          call MPISTOP('Wrong value for location_slit!')
        endif
        if(mype==0) write(*,'(a,f8.1,a)') ' Location of slit: z = ',location_slit,' Unit_length'
      endif

      ! find slit and do integration
      spectra=zero
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&xbmin(^D)=rnode(rpxmin^D_,igrid);
        ^D&xbmax(^D)=rnode(rpxmax^D_,igrid);
        if (location_slit>=xbmin(dir_loc) .and. location_slit<xbmax(dir_loc)) then
          call integrate_spectra_datresol(igrid,wL,dwL,spectra,numWL,numXS,dir_loc,fl)
        endif
      enddo

      numS=numWL*numXS
      call MPI_ALLREDUCE(spectra,spectra_rc,numS,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,icomm,ierrmpi)
      do iwL=1,numWL
        do ixS=1,numXS
          if (spectra_rc(iwL,ixS)>smalldouble) then
            wI(iwL,ixS,1)=spectra_rc(iwL,ixS)
          else
            wI(iwL,ixS,1)=zero
          endif
        enddo
      enddo

      call output_data(qunit,wL,xS,dwL,dxS,wI,numWL,numXS,numWI,datatype)

      deallocate(wL,xS,dwL,dxS,spectra,spectra_rc,wI)

    end subroutine get_spectrum_datresol

    subroutine integrate_spectra_datresol(igrid,wL,dwL,spectra,numWL,numXS,dir_loc,fl)
      use mod_constants

      integer, intent(in) :: igrid,numWL,numXS,dir_loc
      type(te_fluid), intent(in) :: fl
      double precision, intent(in) :: wL(numWL),dwL(numWL)
      double precision, intent(inout) :: spectra(numWL,numXS)

      integer :: direction_LOS
      integer :: ixO^L,ixI^L,ix^D,ixOnew
      double precision, allocatable :: flux(:^D&),v(:^D&),pth(:^D&),Te(:^D&),rho(:^D&)
      double precision :: wlc,wlwd

      integer :: mass
      double precision :: logTe,lineCent
      character (30) :: ion
      double precision :: spaceRsl,wlRsl,sigma_PSF,wslit

      integer :: levelg,rft,ixSmin,ixSmax,iwL
      double precision :: flux_pix,dL

      call get_line_info(spectrum_wl,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)

      if (LOS_phi==0 .and. LOS_theta==90) then
        direction_LOS=1
      else if (LOS_phi==90 .and. LOS_theta==90) then
        direction_LOS=2
      else
        direction_LOS=3
      endif

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      allocate(flux(ixI^S),v(ixI^S),pth(ixI^S),Te(ixI^S),rho(ixI^S))

      ^D&ix^D=ixOmin^D;
      if (dir_loc==1) then
        do ix1=ixOmin1,ixOmax1
          if (location_slit>=(ps(igrid)%x(ix^D,1)-half*ps(igrid)%dx(ix^D,1)) .and. &
              location_slit<(ps(igrid)%x(ix^D,1)+half*ps(igrid)%dx(ix^D,1))) then
            ixOnew=ix1
          endif
        enddo
        ixOmin1=ixOnew
        ixOmax1=ixOnew
      else if (dir_loc==2) then
        do ix2=ixOmin2,ixOmax2
          if (location_slit>=(ps(igrid)%x(ix^D,2)-half*ps(igrid)%dx(ix^D,2)) .and. &
              location_slit<(ps(igrid)%x(ix^D,2)+half*ps(igrid)%dx(ix^D,2))) then
            ixOnew=ix2
          endif
        enddo
        ixOmin2=ixOnew
        ixOmax2=ixOnew
      else
        do ix3=ixOmin3,ixOmax3
          if (location_slit>=(ps(igrid)%x(ix^D,3)-half*ps(igrid)%dx(ix^D,3)) .and. &
              location_slit<(ps(igrid)%x(ix^D,3)+half*ps(igrid)%dx(ix^D,3))) then
            ixOnew=ix3
          endif
        enddo
        ixOmin3=ixOnew
        ixOmax3=ixOnew
      endif

      call get_EUV(spectrum_wl,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
      flux(ixO^S)=flux(ixO^S)/instrument_resolution_factor**2   ! adjust flux due to artifical change of resolution
      call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
      v(ixO^S)=-ps(igrid)%w(ixO^S,iw_mom(direction_LOS))/rho(ixO^S)
      call fl%get_pthermal(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,pth)
      call fl%get_var_Rfactor(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,Te)
      Te(ixO^S)=pth(ixO^S)/(Te(ixO^S)*rho(ixO^S))

      ! grid parameters
      levelg=ps(igrid)%level
      rft=2**(refine_max_level-levelg)

      {do ix^D=ixOmin^D,ixOmax^D\}
        if (flux(ix^D)>smalldouble) then
          if (SI_unit) then
            wlc=lineCent*(1.d0+v(ix^D)*unit_velocity*1.d2/const_c)
          else
            wlc=lineCent*(1.d0+v(ix^D)*unit_velocity/const_c)
          endif
          wlwd=sqrt(kb_cgs*Te(ix^D)*unit_temperature/(mass*mp_cgs))
          wlwd=wlwd*lineCent/const_c
          ! involved pixel
          select case(direction_slit)
          case(1)
            ixSmin=(block_nx1*(node(pig1_,igrid)-1)+(ix1-ixOmin1))*rft+1
            ixSmax=(block_nx1*(node(pig1_,igrid)-1)+(ix1-ixOmin1+1))*rft
          case(2)
            ixSmin=(block_nx2*(node(pig2_,igrid)-1)+(ix2-ixOmin2))*rft+1
            ixSmax=(block_nx2*(node(pig2_,igrid)-1)+(ix2-ixOmin2+1))*rft
          case(3)
            ixSmin=(block_nx3*(node(pig3_,igrid)-1)+(ix3-ixOmin3))*rft+1
            ixSmax=(block_nx3*(node(pig3_,igrid)-1)+(ix3-ixOmin3+1))*rft
          end select
          ! LOS depth
          select case(direction_LOS)
          case(1)
            dL=ps(igrid)%dx(ix^D,1)*unit_length
          case(2)
            dL=ps(igrid)%dx(ix^D,2)*unit_length
          case default
            dL=ps(igrid)%dx(ix^D,3)*unit_length
          end select
          if (SI_unit) dL=dL*1.d2
          ! integral pixel flux
          do iwL=1,numWL
            flux_pix=flux(ix^D)*wlRsl*dL*exp_clamped(-(wL(iwL)-wlc)**2/(2*wlwd**2))/(sqrt(2*dpi)*wlwd)
            if (flux_pix>smalldouble) then
              flux_pix=flux_pix*wslit/spaceRsl
              spectra(iwL,ixSmin:ixSmax)=spectra(iwL,ixSmin:ixSmax)+flux_pix
            endif
          enddo
        endif
      {enddo\}

      deallocate(flux,v,pth,Te,rho)

    end subroutine integrate_spectra_datresol

    subroutine get_spectrum(qunit,datatype,fl)

      integer, intent(in) :: qunit
      character(20), intent(in) :: datatype
      type(te_fluid), intent(in) :: fl

      integer :: numWL,numXS,iwL,ixS,numWI,ix^D
      double precision :: dwLg,dxSg,xSmin,xSmax,xScent,wLmin,wLmax
      double precision, allocatable :: wL(:),xS(:),dwL(:),dxS(:)
      double precision, allocatable :: wI(:,:,:),spectra(:,:),spectra_rc(:,:)
      double precision :: vec_cor(1:3),xI_cor(1:2)
      double precision :: res,r_loc,r_max

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl,wslit
      double precision :: unitv,arcsec,RHESSI_rsl,pixel
      integer :: iigrid,igrid,i,j,numS
      double precision :: xLmin,xLmax,xslit

      if (coordinate==spherical) then
        call init_vectors_spherical()
      else
        ! cartesian
        call init_vectors_cartesian()
      endif

      ! calculate domain in space
      if (coordinate==spherical) then
        xSmin=-abs(xprobmax1)
        xSmax=abs(xprobmax1)
      else
        do ix1=1,2
          if (ix1==1) vec_cor(1)=xprobmin1
          if (ix1==2) vec_cor(1)=xprobmax1
          do ix2=1,2
            if (ix2==1) vec_cor(2)=xprobmin2
            if (ix2==2) vec_cor(2)=xprobmax2
            do ix3=1,2
              if (ix3==1) vec_cor(3)=xprobmin3
              if (ix3==2) vec_cor(3)=xprobmax3
              if (big_image) then
                r_loc=(vec_cor(1)-x_origin(1))**2
                r_loc=r_loc+(vec_cor(2)-x_origin(2))**2
                r_loc=r_loc+(vec_cor(3)-x_origin(3))**2
                r_loc=sqrt(r_loc)
                if (ix1==1 .and. ix2==1 .and. ix3==1) then
                  r_max=r_loc
                else
                  r_max=max(r_max,r_loc)
                endif
              else
                call get_cor_image(vec_cor,xI_cor)
                if (ix1==1 .and. ix2==1 .and. ix3==1) then
                  xSmin=xI_cor(2)
                  xSmax=xI_cor(2)
                else
                  xSmin=min(xSmin,xI_cor(2))
                  xSmax=max(xSmax,xI_cor(2))
                endif
              endif
            enddo
          enddo
        enddo
        if (big_image) then
          xSmin=-r_max
          xSmax=r_max
        endif
      endif
      xScent=(xSmin+xSmax)/2.d0

      ! tables for storing spectra data
      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif
      call get_line_info(spectrum_wl,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
      dxSg=spaceRsl*arcsec
      numXS=ceiling((xSmax-xScent)/dxSg)
      xSmin=xScent-numXS*dxSg
      xSmax=xScent+numXS*dxSg
      numXS=numXS*2
      dwLg=wlRsl
      numWL=2*int((spectrum_window_max-spectrum_window_min)/(2.d0*dwLg))
      wLmin=(spectrum_window_max+spectrum_window_min)/2.d0-dwLg*numWL/2
      wLmax=(spectrum_window_max+spectrum_window_min)/2.d0+dwLg*numWL/2
      allocate(wL(numWL),dwL(numWL),xS(numXS),dxS(numXS))
      numWI=1
      allocate(wI(numWL,numXS,numWI),spectra(numWL,numXS),spectra_rc(numWL,numXS))
      do iWL=1,numWL
        wL(iwL)=wLmin+iwL*dwLg-half*dwLg
        dwL=dwLg
      enddo
      do ixS=1,numXS
        xS(ixS)=xSmin+dxSg*(ixS-half)
        dxS(ixS)=dxSg
      enddo

      ! find slit and do integration
      spectra=zero
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        do ix1=1,2
          if (ix1==1) vec_cor(1)=rnode(rpxmin1_,igrid)
          if (ix1==2) vec_cor(1)=rnode(rpxmax1_,igrid)
          do ix2=1,2
            if (ix2==1) vec_cor(2)=rnode(rpxmin2_,igrid)
            if (ix2==2) vec_cor(2)=rnode(rpxmax2_,igrid)
            do ix3=1,2
              if (ix3==1) vec_cor(3)=rnode(rpxmin3_,igrid)
              if (ix3==2) vec_cor(3)=rnode(rpxmax3_,igrid)
              call get_cor_image(vec_cor,xI_cor)
              if (ix1==1 .and. ix2==1 .and. ix3==1) then
                xLmin=xI_cor(1)
                xLmax=xI_cor(1)
              else
                xLmin=min(xLmin,xI_cor(1))
                xLmax=max(xLmax,xI_cor(1))
              endif
            enddo
          enddo
        enddo

        if (activate_unit_arcsec) then 
          xslit=location_slit*arcsec
        else
          xslit=location_slit
        endif
        if (xslit>=xLmin-wslit*arcsec .and. xslit<=xLmax+wslit*arcsec) then
          call integrate_spectra_cartesian(igrid,wL,dwLg,xS,dxSg,spectra,numWL,numXS,fl)
        endif
      enddo

      numS=numWL*numXS
      call MPI_ALLREDUCE(spectra,spectra_rc,numS,MPI_DOUBLE_PRECISION, &
                         MPI_SUM,icomm,ierrmpi)
      do iwL=1,numWL
        do ixS=1,numXS
          if (spectra_rc(iwL,ixS)>smalldouble) then
            wI(iwL,ixS,1)=spectra_rc(iwL,ixS)
          else
            wI(iwL,ixS,1)=zero
          endif
        enddo
      enddo

      if (activate_unit_arcsec) then 
        xS=xS/arcsec
        dxS=dxS/arcsec
      endif

      call output_data(qunit,wL,xS,dwL,dxS,wI,numWL,numXS,numWI,datatype)

      deallocate(wL,xS,dwL,dxS,spectra,spectra_rc,wI)

    end subroutine get_spectrum

    subroutine integrate_spectra_cartesian(igrid,wL,dwLg,xS,dxSg,spectra,numWL,numXS,fl)

      integer, intent(in) :: igrid,numWL,numXS
      double precision, intent(in) :: wL(numWL),xS(numXS)
      double precision, intent(in) :: dwLg,dxSg
      double precision, intent(inout) :: spectra(numWL,numXS)
      type(te_fluid), intent(in) :: fl

      integer :: ixO^L,ixI^L,ix^D,ixOnew,j
      double precision, allocatable :: flux(:^D&),v(:^D&),pth(:^D&),Te(:^D&),rho(:^D&)
      double precision :: wlc,wlwd,res,dst_slit,xslit,arcsec
      double precision :: vloc(1:3),xloc(1:3),dxloc(1:3),xIloc(1:2),dxIloc(1:2)
      integer :: nSubC^D,iSubC^D,iwL,ixS,ixSmin,ixSmax,iwLmin,iwLmax,nwL
      double precision :: slit_width,dxSubC^D,xerf^L,fluxSubC
      double precision :: xSubC(1:3),xCent(1:2)

      integer :: mass
      double precision :: logTe,lineCent
      character (30) :: ion
      double precision :: spaceRsl,wlRsl,sigma_PSF,wslit
      double precision :: sigma_wl,sigma_xs,factor

      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif
      if (activate_unit_arcsec) then 
        xslit=location_slit*arcsec
      else
        xslit=location_slit
      endif

      call get_line_info(spectrum_wl,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      allocate(flux(ixI^S),v(ixI^S),pth(ixI^S),Te(ixI^S),rho(ixI^S))
      ! get local EUV flux and velocity
      call get_EUV(spectrum_wl,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
      flux(ixO^S)=flux(ixO^S)/instrument_resolution_factor**2   ! adjust flux due to artifical change of resolution
      call fl%get_pthermal(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,pth)
      call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
      call fl%get_var_Rfactor(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,Te)
      Te(ixO^S)=pth(ixO^S)/(Te(ixO^S)*rho(ixO^S))
      {do ix^D=ixOmin^D,ixOmax^D\}
        do j=1,3
          vloc(j)=ps(igrid)%w(ix^D,iw_mom(j))/rho(ix^D)
        enddo
        call dot_product_loc(vloc,vec_LOS,res)
        v(ix^D)=res
      {enddo\}

      deallocate(rho)

      slit_width=wslit*arcsec
      sigma_wl=sigma_PSF*dwLg
      sigma_xs=sigma_PSF*dxSg
      {do ix^D=ixOmin^D,ixOmax^D\}
        if (flux(ix^D)>smalldouble) then
          xloc(1:3)=ps(igrid)%x(ix^D,1:3)
          dxloc(1:3)=ps(igrid)%dx(ix^D,1:3)
          call get_cor_image(xloc,xIloc)
          call dot_product_loc(dxloc,vec_xI1,res)
          dxIloc(1)=abs(res)
          if (xIloc(1)>=xslit-half*(slit_width+dxIloc(1)) .and. & 
              xIloc(1)<=xslit+half*(slit_width+dxIloc(1))) then
            ^D&nSubC^D=1;
            ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI1(^D))/(slit_width/16.d0)));
            ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI2(^D))/(dxSg/4.d0)));
            ^D&dxSubC^D=ps(igrid)%dx(ix^DD,^D)/nSubC^D;
            ! local line center and line width
            if (SI_unit) then
              fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length*1.d2/dxSg/dxSg  ! DN s^-1
              wlc=lineCent*(1.d0+v(ix^D)*unit_velocity*1.d2/const_c)
            else
              fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length/dxSg/dxSg  ! DN s^-1
              wlc=lineCent*(1.d0+v(ix^D)*unit_velocity/const_c)
            endif
            wlwd=sqrt(kb_cgs*Te(ix^D)*unit_temperature/(mass*mp_cgs))
            wlwd=wlwd*lineCent/const_c
            ! dividing a cell to several parts to get more accurate integrating values
            {do iSubC^D=1,nSubC^D\}
              ^D&xSubC(^D)=xloc(^D)-half*dxloc(^D)+(iSubC^D-half)*dxSubC^D;
              call get_cor_image(xSubC,xCent)
              dst_slit=abs(xCent(1)-xslit)  ! space distance to slit center
              if (dst_slit<=half*slit_width) then
                ixS=floor((xCent(2)-(xS(1)-half*dxSg))/dxSg)+1
                ixSmin=max(1,ixS-3)
                ixSmax=min(ixS+3,numXS)
                iwL=floor((wlc-(wL(1)-half*dwLg))/dwLg)+1
                nwL=3*ceiling(wlwd/dwLg+1)
                iwLmin=max(1,iwL-nwL)
                iwLmax=min(iwL+nwL,numWL)
                ! calculate the contribution to nearby pixels
                do iwL=iwLmin,iwLmax
                  do ixS=ixSmin,ixSmax
                    xerfmin1=(wL(iwL)-half*dwLg-wlc)/sqrt(2.d0*(sigma_wl**2+wlwd**2))
                    xerfmax1=(wL(iwL)+half*dwLg-wlc)/sqrt(2.d0*(sigma_wl**2+wlwd**2))
                    xerfmin2=(xS(ixS)-half*dxSg-xCent(2))/(sqrt(2.d0)*sigma_xs)
                    xerfmax2=(xS(ixS)+half*dxSg-xCent(2))/(sqrt(2.d0)*sigma_xs)
                    factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
                    spectra(iwL,ixS)=spectra(iwL,ixS)+fluxSubC*factor
                  enddo
                enddo
                ! nearby pixels
              endif
            {enddo\}
          endif
        endif
      {enddo\}

      deallocate(flux,v,pth,Te)
    end subroutine integrate_spectra_cartesian
  }

  {^IFTHREED
    subroutine get_EUV_image(qunit,fl)
      use mod_global_parameters

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20) :: datatype

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl,wslit
      double precision :: t0,t1

      t0=MPI_WTIME()
      datatype='image_euv'
      call check_synthetic_emission_options(datatype)
      call get_line_info(wavelength,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)

      if (mype==0) then
        print *, '###################################################'
        print *, 'Systhesizing EUV image'
        write(*,'(a,f8.3,a)') ' Wavelength: ',lineCent,' Angstrom'
        print *, 'Unit of EUV flux: DN s^-1 pixel^-1'
      endif

      if (dat_resolution) then
        if (.not. slab .and. .not. (coordinate==spherical .and. trim(ray_method_active)=='spherical')) &
          call MPISTOP('EUV dat-resolution needs Cartesian or spherical native rays')
        if (mype==0) then
          print *, 'Data-resolution image requested; native output pixel sizes are reported below.'
          if (SI_unit) then
            write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
          else
            write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
          endif
        endif
        if (coordinate==spherical .and. trim(ray_method_active)=='spherical') then
          call get_image_datresol(qunit,datatype,fl)
        else if (trim(ray_method_active)=='cart') then
          call get_image_datresol(qunit,datatype,fl)
        else if (LOS_phi==0 .and. LOS_theta==90) then
          call get_image_datresol(qunit,datatype,fl)
        else if (LOS_phi==90 .and. LOS_theta==90) then
          call get_image_datresol(qunit,datatype,fl)
        else if (LOS_theta==0) then
          call get_image_datresol(qunit,datatype,fl)
        else
          call MPISTOP('ERROR: Wrong LOS for synthesizing emission!')
        endif
      else
        if (mype==0) then
          write(*,'(a,f7.1,a,f7.1,a,f5.1,a,f5.1,a)') ' Pixel: ',spaceRsl*725.0,' km x ',spaceRsl*725.0, ' km  (', &
                                                                  spaceRsl, ' arcsec x ', spaceRsl, ' arcsec)'
          if (activate_unit_arcsec) then
            print *, 'Unit of length: arcsec (~725 km)'
          else
            if (SI_unit) then
              write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
            else
              write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
            endif
          endif
        endif
        if (coordinate==cartesian) then
          if (mype==0) write(*,'(a,f8.3,f8.3,f8.3,a)') ' Mapping: [',x_origin(1),x_origin(2),x_origin(3), &
                                                     '] of the simulation box is located at [X=0,Y=0] of the image'
            call get_image(qunit,datatype,fl)
        else if (coordinate==spherical) then
          if (mype==0) write(*,'(a,f6.3,f8.3,f8.3,a)') ' Mapping: R=0 of the simulation box is located at [X=0,Y=0] of the image'
          call get_image(qunit,datatype,fl)
        else
          call MPISTOP("EUV synthesis: this coordinate is not supported!")
        endif
      endif

      t1=MPI_WTIME()
      if (mype==0) print *, 'time comsuming: ',t1-t0,' s'
      if (mype==0) print *, '###################################################'
      
    end subroutine get_EUV_image

    subroutine get_SXR_image(qunit,fl)
      use mod_global_parameters

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20) :: datatype
      double precision :: RHESSI_rsl
      double precision :: t0,t1

      t0=MPI_WTIME()
      datatype='image_sxr'
      call check_synthetic_emission_options(datatype)
      RHESSI_rsl=2.3/instrument_resolution_factor

      if (mype==0) then
        print *, '###################################################'
        print *, 'Systhesizing SXR image (observed at 1 AU).'
        write(*,'(a,i2,a,i2,a)') ' Passband: ',emin_sxr,' - ',emax_sxr,' keV'
      endif

      if (dat_resolution) then
        if (coordinate/=cartesian) call MPISTOP('SXR synthesis: only cartesian is supported for .dat resolution!')
        if (mype==0) then
          print *, 'Unit of SXR flux: photons cm^-2 s^-1 pixel^-1'
          write(*,'(a,f5.1,a,f5.1,a,f5.1,a,f5.1,a)') ' Supposed Pixel: ',RHESSI_rsl*0.725, ' Mm x ',RHESSI_rsl*0.725, &
                                                     ' Mm  (', RHESSI_rsl, ' arcsec x ', RHESSI_rsl, ' arcsec)'
          if (SI_unit) then
            write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
          else
            write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
          endif
        endif
        if (LOS_phi==0 .and. LOS_theta==90) then
          call get_image_datresol(qunit,datatype,fl)
        else if (LOS_phi==90 .and. LOS_theta==90) then
          call get_image_datresol(qunit,datatype,fl)
        else if (LOS_theta==0) then
          call get_image_datresol(qunit,datatype,fl)
        else
          call MPISTOP('ERROR: Wrong LOS for synthesizing emission!')
        endif
      else
        if (mype==0) then
          print *, 'Unit of SXR flux: photons cm^-2 s^-1 pixel^-1'
          write(*,'(a,f5.1,a,f5.1,a,f5.1,a,f5.1,a)') ' Pixel: ',RHESSI_rsl*0.725, ' Mm x ',RHESSI_rsl*0.725, &
                                                      ' Mm  (', RHESSI_rsl, ' arcsec x ', RHESSI_rsl, ' arcsec)'
          if (activate_unit_arcsec) then
            print *, 'Unit of length: arcsec (~725 km)'
          else
            if (SI_unit) then
              write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
            else
              write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
            endif
          endif
        endif
        if (coordinate==cartesian) then
          if (mype==0) write(*,'(a,f8.3,f8.3,f8.3,a)') ' Mapping: [',x_origin(1),x_origin(2),x_origin(3), &
                                                     '] of the simulation box is located at [X=0,Y=0] of the image'
            call get_image(qunit,datatype,fl)
        else if (coordinate==spherical) then
          if (mype==0) write(*,'(a,f6.3,f8.3,f8.3,a)') ' Mapping: R=0 of the simulation box is located at [X=0,Y=0] of the image'
          call get_image(qunit,datatype,fl)
        else
          call MPISTOP("SXR synthesis: this coordinate is not supported!")
        endif
      endif

      t1=MPI_WTIME()
      if (mype==0) print *, 'time comsuming:',t1-t0
      if (mype==0) print *, '###################################################'

    end subroutine get_SXR_image

    subroutine get_whitelight_image(qunit,fl)
      use mod_global_parameters

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20) :: datatype
      double precision :: LASCO_rsl

      if (mype==0) print *, '###################################################'

      if (whitelight_instrument=='LASCO/C1') then
        LASCO_rsl=5.6d0/instrument_resolution_factor
        if (mype==0) print *, 'Systhesizing white light image (observed by LASCO/C1).'
      else if (whitelight_instrument=='LASCO/C2') then
        LASCO_rsl=11.4d0/instrument_resolution_factor
        if (mype==0) print *, 'Systhesizing white light image (observed by LASCO/C2).'
      else if (whitelight_instrument=='LASCO/C3') then
        LASCO_rsl=56.d0/instrument_resolution_factor
        if (mype==0) print *, 'Systhesizing white light image (observed by LASCO/C3).'
      else
        call MPISTOP('Whitelight synthesis: instrument is not supported!')
      endif

      if (mype==0) write(*,'(a,f5.1,a,f5.1,a,f5.1,a,f5.1,a)') ' Pixel: ',LASCO_rsl*0.725,' Mm x ',LASCO_rsl*0.725, ' Mm  (', &
                                                                LASCO_rsl, ' arcsec x ', LASCO_rsl, ' arcsec) '
      if (mype==0) print *, 'Unit of white light flux: average Sun brightness'

      datatype='image_whitelight'
      call check_synthetic_emission_options(datatype)

      if (mype==0) then
        if (activate_unit_arcsec) then
          print *, 'Unit of length: arcsec (~725 km)'
        else
          if (SI_unit) then
            if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
          else
            if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
          endif
        endif
      endif

      if (coordinate==spherical) then
        if (mype==0) write(*,'(a,f6.3,f8.3,f8.3,a)') ' Mapping: R=0 of the simulation box is located at [X=0,Y=0] of the image'
        call get_image(qunit,datatype,fl)
      else
        call MPISTOP("Whitelight synthesis: this coordinate is not supported!")
      endif

      if (mype==0) print *, '###################################################'

    end subroutine get_whitelight_image

    subroutine postprocess_euv_instrument_image(nSrc1,nSrc2,xSrc1,xSrc2,dxSrc1,dxSrc2,&
                                                EUV,Dpl,nOut1,nOut2,xOut1,xOut2,&
                                                dxOut1,dxOut2,wOut,numWOut,Tau,EUVthin)
      use mod_global_parameters
      use mod_constants

      integer, intent(in) :: nSrc1,nSrc2
      double precision, intent(in) :: xSrc1(nSrc1),xSrc2(nSrc2)
      double precision, intent(in) :: dxSrc1(nSrc1),dxSrc2(nSrc2)
      double precision, intent(in) :: EUV(nSrc1,nSrc2),Dpl(nSrc1,nSrc2)
      integer, intent(out) :: nOut1,nOut2,numWOut
      double precision, allocatable, intent(out) :: xOut1(:),xOut2(:),dxOut1(:),dxOut2(:)
      double precision, allocatable, intent(out) :: wOut(:,:,:)
      double precision, intent(in), optional :: Tau(nSrc1,nSrc2),EUVthin(nSrc1,nSrc2)

      integer :: mass,ixS1,ixS2,ixP1,ixP2,ixC1,ixC2,iw
      integer :: ixPmin1,ixPmax1,ixPmin2,ixPmax2
      character(30) :: ion
      double precision :: logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit
      double precision :: arcsec,dxInst,xMin1,xMax1,xMin2,xMax2,xCent1,xCent2
      double precision :: sigma0,xerfmin1,xerfmax1,xerfmin2,xerfmax2
      double precision :: factor,weightSum,weightNorm,thinVal,tauVal
      double precision, allocatable :: dplNum(:,:),thinOut(:,:),tauOut(:,:),tauWeight(:,:)

      call get_line_info(wavelength,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif
      dxInst=spaceRsl*arcsec
      if (dxInst<=zero) call mpistop("instrument_postprocess has non-positive pixel size")

      xMin1=minval(xSrc1-half*dxSrc1)
      xMax1=maxval(xSrc1+half*dxSrc1)
      xMin2=minval(xSrc2-half*dxSrc2)
      xMax2=maxval(xSrc2+half*dxSrc2)
      xCent1=half*(xMin1+xMax1)
      xCent2=half*(xMin2+xMax2)
      nOut1=16*max(1,ceiling((xMax1-xMin1)/(16.d0*dxInst)))
      nOut2=16*max(1,ceiling((xMax2-xMin2)/(16.d0*dxInst)))
      xMin1=xCent1-half*dble(nOut1)*dxInst
      xMin2=xCent2-half*dble(nOut2)*dxInst

      allocate(xOut1(nOut1),xOut2(nOut2),dxOut1(nOut1),dxOut2(nOut2))
      do ixP1=1,nOut1
        xOut1(ixP1)=xMin1+dxInst*(dble(ixP1)-half)
        dxOut1(ixP1)=dxInst
      enddo
      do ixP2=1,nOut2
        xOut2(ixP2)=xMin2+dxInst*(dble(ixP2)-half)
        dxOut2(ixP2)=dxInst
      enddo

      numWOut=2
      if (present(Tau) .and. output_tau) numWOut=numWOut+1
      if (present(EUVthin) .and. output_absorption_fraction) numWOut=numWOut+1
      allocate(wOut(nOut1,nOut2,numWOut),dplNum(nOut1,nOut2))
      wOut=zero
      dplNum=zero
      if (present(EUVthin)) then
        allocate(thinOut(nOut1,nOut2))
        thinOut=zero
      endif
      if (present(Tau) .and. output_tau) then
        allocate(tauOut(nOut1,nOut2),tauWeight(nOut1,nOut2))
        tauOut=zero
        tauWeight=zero
      endif

      sigma0=sigma_PSF*dxInst
      do ixS1=1,nSrc1
        do ixS2=1,nSrc2
          thinVal=zero
          tauVal=zero
          if (present(EUVthin)) thinVal=EUVthin(ixS1,ixS2)
          if (present(Tau)) tauVal=Tau(ixS1,ixS2)
          if (abs(EUV(ixS1,ixS2))<=smalldouble .and. abs(thinVal)<=smalldouble .and. &
              abs(tauVal)<=smalldouble) cycle

          ixC1=floor((xSrc1(ixS1)-(xOut1(1)-half*dxInst))/dxInst)+1
          ixC2=floor((xSrc2(ixS2)-(xOut2(1)-half*dxInst))/dxInst)+1
          ixPmin1=max(1,ixC1-3)
          ixPmax1=min(nOut1,ixC1+3)
          ixPmin2=max(1,ixC2-3)
          ixPmax2=min(nOut2,ixC2+3)

          weightSum=zero
          do ixP1=ixPmin1,ixPmax1
            do ixP2=ixPmin2,ixPmax2
              xerfmin1=((xOut1(ixP1)-half*dxInst)-xSrc1(ixS1))/(sqrt(2.d0)*sigma0)
              xerfmax1=((xOut1(ixP1)+half*dxInst)-xSrc1(ixS1))/(sqrt(2.d0)*sigma0)
              xerfmin2=((xOut2(ixP2)-half*dxInst)-xSrc2(ixS2))/(sqrt(2.d0)*sigma0)
              xerfmax2=((xOut2(ixP2)+half*dxInst)-xSrc2(ixS2))/(sqrt(2.d0)*sigma0)
              factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
              weightSum=weightSum+factor
            enddo
          enddo
          if (weightSum<=zero) cycle

          do ixP1=ixPmin1,ixPmax1
            do ixP2=ixPmin2,ixPmax2
              xerfmin1=((xOut1(ixP1)-half*dxInst)-xSrc1(ixS1))/(sqrt(2.d0)*sigma0)
              xerfmax1=((xOut1(ixP1)+half*dxInst)-xSrc1(ixS1))/(sqrt(2.d0)*sigma0)
              xerfmin2=((xOut2(ixP2)-half*dxInst)-xSrc2(ixS2))/(sqrt(2.d0)*sigma0)
              xerfmax2=((xOut2(ixP2)+half*dxInst)-xSrc2(ixS2))/(sqrt(2.d0)*sigma0)
              factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
              weightNorm=factor/weightSum
              wOut(ixP1,ixP2,1)=wOut(ixP1,ixP2,1)+EUV(ixS1,ixS2)*weightNorm
              dplNum(ixP1,ixP2)=dplNum(ixP1,ixP2)+EUV(ixS1,ixS2)*Dpl(ixS1,ixS2)*weightNorm
              if (present(EUVthin)) thinOut(ixP1,ixP2)=thinOut(ixP1,ixP2)+thinVal*weightNorm
              if (present(Tau) .and. output_tau) then
                tauOut(ixP1,ixP2)=tauOut(ixP1,ixP2)+tauVal*weightNorm
                tauWeight(ixP1,ixP2)=tauWeight(ixP1,ixP2)+weightNorm
              endif
            enddo
          enddo
        enddo
      enddo

      do ixP1=1,nOut1
        do ixP2=1,nOut2
          if (wOut(ixP1,ixP2,1)>smalldouble) then
            wOut(ixP1,ixP2,2)=dplNum(ixP1,ixP2)/wOut(ixP1,ixP2,1)
          else
            wOut(ixP1,ixP2,2)=zero
          endif
        enddo
      enddo
      iw=2
      if (present(Tau) .and. output_tau) then
        iw=iw+1
        do ixP1=1,nOut1
          do ixP2=1,nOut2
            if (tauWeight(ixP1,ixP2)>zero) then
              wOut(ixP1,ixP2,iw)=tauOut(ixP1,ixP2)/tauWeight(ixP1,ixP2)
            else
              wOut(ixP1,ixP2,iw)=zero
            endif
          enddo
        enddo
      endif
      if (present(EUVthin) .and. output_absorption_fraction) then
        iw=iw+1
        do ixP1=1,nOut1
          do ixP2=1,nOut2
            if (thinOut(ixP1,ixP2)>smalldouble) then
              wOut(ixP1,ixP2,iw)=min(one,max(zero,(thinOut(ixP1,ixP2)-wOut(ixP1,ixP2,1))/thinOut(ixP1,ixP2)))
            else
              wOut(ixP1,ixP2,iw)=zero
            endif
          enddo
        enddo
      endif

      if (mype==0) then
        write(*,'(a,2(i8,1x),a,2(i8,1x),a,1pe12.5)') &
          ' instrument_postprocess EUV grid src/out: ',nSrc1,nSrc2,' -> ',nOut1,nOut2,' dx=',dxInst
      endif

      deallocate(dplNum)
      if (allocated(thinOut)) deallocate(thinOut)
      if (allocated(tauOut)) deallocate(tauOut,tauWeight)
    end subroutine postprocess_euv_instrument_image

    subroutine postprocess_radio_beam_image(nSrc1,nSrc2,xSrc1,xSrc2,dxSrc1,dxSrc2,&
                                            Bright,nOut1,nOut2,xOut1,xOut2,&
                                            dxOut1,dxOut2,wOut,numWOut,Tau,BrightThin)
      use mod_global_parameters

      integer, intent(in) :: nSrc1,nSrc2
      double precision, intent(in) :: xSrc1(nSrc1),xSrc2(nSrc2)
      double precision, intent(in) :: dxSrc1(nSrc1),dxSrc2(nSrc2)
      double precision, intent(in) :: Bright(nSrc1,nSrc2)
      integer, intent(out) :: nOut1,nOut2,numWOut
      double precision, allocatable, intent(out) :: xOut1(:),xOut2(:),dxOut1(:),dxOut2(:)
      double precision, allocatable, intent(out) :: wOut(:,:,:)
      double precision, intent(in), optional :: Tau(nSrc1,nSrc2),BrightThin(nSrc1,nSrc2)

      integer :: ixS1,ixS2,ixP1,ixP2,ixC1,ixC2,iw,nStencil
      integer :: ixPmin1,ixPmax1,ixPmin2,ixPmax2
      double precision :: arcsec,beamPixel,beamSigma,xMin1,xMax1,xMin2,xMax2,xCent1,xCent2
      double precision :: distance1,distance2,weight,cellArea,thinVal,tauVal
      double precision, allocatable :: norm(:,:),thinOut(:,:),tauOut(:,:),tauNorm(:,:)

      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif
      beamSigma=radio_beam_fwhm*arcsec/sqrt(8.d0*log(2.d0))
      if (radio_beam_pixel_size>zero) then
        beamPixel=radio_beam_pixel_size*arcsec
      else
        beamPixel=radio_beam_fwhm*arcsec/3.d0
      endif
      if (beamSigma<=zero .or. beamPixel<=zero) then
        call mpistop("radio beam postprocess has non-positive beam or pixel size")
      endif

      xMin1=minval(xSrc1-half*dxSrc1)
      xMax1=maxval(xSrc1+half*dxSrc1)
      xMin2=minval(xSrc2-half*dxSrc2)
      xMax2=maxval(xSrc2+half*dxSrc2)
      xCent1=half*(xMin1+xMax1)
      xCent2=half*(xMin2+xMax2)
      nOut1=16*max(1,ceiling((xMax1-xMin1)/(16.d0*beamPixel)))
      nOut2=16*max(1,ceiling((xMax2-xMin2)/(16.d0*beamPixel)))
      xMin1=xCent1-half*dble(nOut1)*beamPixel
      xMin2=xCent2-half*dble(nOut2)*beamPixel

      allocate(xOut1(nOut1),xOut2(nOut2),dxOut1(nOut1),dxOut2(nOut2))
      do ixP1=1,nOut1
        xOut1(ixP1)=xMin1+beamPixel*(dble(ixP1)-half)
        dxOut1(ixP1)=beamPixel
      enddo
      do ixP2=1,nOut2
        xOut2(ixP2)=xMin2+beamPixel*(dble(ixP2)-half)
        dxOut2(ixP2)=beamPixel
      enddo

      numWOut=1
      if (present(Tau) .and. output_tau) numWOut=numWOut+1
      if (present(BrightThin) .and. output_absorption_fraction) numWOut=numWOut+1
      allocate(wOut(nOut1,nOut2,numWOut),norm(nOut1,nOut2))
      wOut=zero
      norm=zero
      if (present(BrightThin)) then
        allocate(thinOut(nOut1,nOut2))
        thinOut=zero
      endif
      if (present(Tau) .and. output_tau) then
        allocate(tauOut(nOut1,nOut2),tauNorm(nOut1,nOut2))
        tauOut=zero
        tauNorm=zero
      endif

      nStencil=max(3,ceiling(4.d0*beamSigma/beamPixel)+1)
      do ixS1=1,nSrc1
        do ixS2=1,nSrc2
          thinVal=zero
          tauVal=zero
          if (present(BrightThin)) thinVal=BrightThin(ixS1,ixS2)
          if (present(Tau)) tauVal=Tau(ixS1,ixS2)
          if (abs(Bright(ixS1,ixS2))<=smalldouble .and. abs(thinVal)<=smalldouble .and. &
              abs(tauVal)<=smalldouble) cycle

          ixC1=floor((xSrc1(ixS1)-(xOut1(1)-half*beamPixel))/beamPixel)+1
          ixC2=floor((xSrc2(ixS2)-(xOut2(1)-half*beamPixel))/beamPixel)+1
          ixPmin1=max(1,ixC1-nStencil)
          ixPmax1=min(nOut1,ixC1+nStencil)
          ixPmin2=max(1,ixC2-nStencil)
          ixPmax2=min(nOut2,ixC2+nStencil)
          cellArea=max(smalldouble,dxSrc1(ixS1)*dxSrc2(ixS2))

          do ixP1=ixPmin1,ixPmax1
            distance1=xOut1(ixP1)-xSrc1(ixS1)
            do ixP2=ixPmin2,ixPmax2
              distance2=xOut2(ixP2)-xSrc2(ixS2)
              weight=exp_clamped(-half*(distance1**2+distance2**2)/beamSigma**2)*cellArea
              wOut(ixP1,ixP2,1)=wOut(ixP1,ixP2,1)+Bright(ixS1,ixS2)*weight
              norm(ixP1,ixP2)=norm(ixP1,ixP2)+weight
              if (present(BrightThin)) thinOut(ixP1,ixP2)=thinOut(ixP1,ixP2)+thinVal*weight
              if (present(Tau) .and. output_tau) then
                tauOut(ixP1,ixP2)=tauOut(ixP1,ixP2)+tauVal*weight
                tauNorm(ixP1,ixP2)=tauNorm(ixP1,ixP2)+weight
              endif
            enddo
          enddo
        enddo
      enddo

      do ixP1=1,nOut1
        do ixP2=1,nOut2
          if (norm(ixP1,ixP2)>zero) then
            wOut(ixP1,ixP2,1)=wOut(ixP1,ixP2,1)/norm(ixP1,ixP2)
            if (present(BrightThin)) thinOut(ixP1,ixP2)=thinOut(ixP1,ixP2)/norm(ixP1,ixP2)
          else
            wOut(ixP1,ixP2,1)=zero
            if (present(BrightThin)) thinOut(ixP1,ixP2)=zero
          endif
        enddo
      enddo

      iw=1
      if (present(Tau) .and. output_tau) then
        iw=iw+1
        do ixP1=1,nOut1
          do ixP2=1,nOut2
            if (tauNorm(ixP1,ixP2)>zero) then
              wOut(ixP1,ixP2,iw)=tauOut(ixP1,ixP2)/tauNorm(ixP1,ixP2)
            else
              wOut(ixP1,ixP2,iw)=zero
            endif
          enddo
        enddo
      endif
      if (present(BrightThin) .and. output_absorption_fraction) then
        iw=iw+1
        do ixP1=1,nOut1
          do ixP2=1,nOut2
            if (thinOut(ixP1,ixP2)>smalldouble) then
              wOut(ixP1,ixP2,iw)=min(one,max(zero,(thinOut(ixP1,ixP2)-wOut(ixP1,ixP2,1))/thinOut(ixP1,ixP2)))
            else
              wOut(ixP1,ixP2,iw)=zero
            endif
          enddo
        enddo
      endif

      if (mype==0) then
        write(*,'(a,2(i8,1x),a,2(i8,1x),a,2(1pe12.5,1x))') &
          ' radio_beam_postprocess grid src/out: ',nSrc1,nSrc2,' -> ',nOut1,nOut2,&
          ' fwhm/pixel=',radio_beam_fwhm,beamPixel/arcsec
      endif

      deallocate(norm)
      if (allocated(thinOut)) deallocate(thinOut)
      if (allocated(tauOut)) deallocate(tauOut,tauNorm)
    end subroutine postprocess_radio_beam_image

    subroutine get_image_datresol(qunit,datatype,fl)
      ! integrate emission flux along line of sight (LOS)
      ! in a 3D simulation box and get a 2D EUV image
      use mod_global_parameters
      use mod_constants

      integer, intent(in) :: qunit
      character(20), intent(in) :: datatype
      type(te_fluid), intent(in) :: fl

      double precision :: dx^D
      integer :: numX^D,ix^D
      double precision, allocatable :: EUV(:,:),EUVs(:,:),Dpl(:,:),Dpls(:,:)
      double precision, allocatable :: EUVthin(:,:),Tau(:,:)
      double precision, allocatable :: SXR(:,:),SXRs(:,:),wI(:,:,:)
      double precision, allocatable :: xI1(:),xI2(:),dxI1(:),dxI2(:),dxIi
      integer :: numXI1,numXI2,numSI,numWI,iw
      double precision :: xI^L
      integer :: iigrid,igrid,i,j
      double precision, allocatable :: xIF1(:),xIF2(:),dxIF1(:),dxIF2(:)
      double precision, allocatable :: xIP1(:),xIP2(:),dxIP1(:),dxIP2(:),wIP(:,:,:)
      integer :: nXIF1,nXIF2
      integer :: nXIP1,nXIP2,numWIP
      double precision :: xIF^L
      double precision :: vec_cor(1:3),xI_cor(1:2),dxDDA,xIcent1,xIcent2

      double precision :: unitv,arcsec,RHESSI_rsl,length_to_km
      integer :: strtype^D,nstrb^D,nbb^D,nuni^D,nstr^D,bnx^D
      double precision :: qs^D,dxfirst^D,dxmid^D,lenstr^D
      logical :: has_doppler_output,has_thick_output

      numX1=domain_nx1*2**(refine_max_level-1)
      numX2=domain_nx2*2**(refine_max_level-1)
      numX3=domain_nx3*2**(refine_max_level-1)

      if (trim(ray_method_active)=='cart') call init_vectors_cartesian()
      if (coordinate==spherical .and. trim(ray_method_active)=='spherical') call init_vectors_spherical()

      ! parameters for creating table
      if (coordinate==spherical .and. trim(ray_method_active)=='spherical') then
        call get_sph_intersection_image_bounds(xIFmin1,xIFmax1,xIFmin2,xIFmax2)
        call get_sph_intersection_datresol_spacing(dxDDA)
        xIcent1=half*(xIFmin1+xIFmax1)
        xIcent2=half*(xIFmin2+xIFmax2)
        nXIF1=max(1,ceiling((xIFmax1-xIFmin1)/dxDDA))
        nXIF2=max(1,ceiling((xIFmax2-xIFmin2)/dxDDA))
        xIFmin1=xIcent1-half*dble(nXIF1)*dxDDA
        xIFmax1=xIcent1+half*dble(nXIF1)*dxDDA
        xIFmin2=xIcent2-half*dble(nXIF2)*dxDDA
        xIFmax2=xIcent2+half*dble(nXIF2)*dxDDA
        bnx1=1
        bnx2=1
        nbb1=nXIF1
        nbb2=nXIF2
        strtype1=0
        strtype2=0
        nstrb1=0
        nstrb2=0
        qs1=one
        qs2=one
        if (mype==0) write(*,'(a,1pe12.5,a,2(i8,1x))') &
          ' spherical native dat-resolution image-plane dx: ',dxDDA,' n=',nXIF1,nXIF2
      else if (trim(ray_method_active)=='cart' .and. .not. &
          ((LOS_phi==0 .and. LOS_theta==90) .or. &
           (LOS_phi==90 .and. LOS_theta==90) .or. LOS_theta==0)) then
        do ix1=1,2
          if (ix1==1) vec_cor(1)=xprobmin1
          if (ix1==2) vec_cor(1)=xprobmax1
          do ix2=1,2
            if (ix2==1) vec_cor(2)=xprobmin2
            if (ix2==2) vec_cor(2)=xprobmax2
            do ix3=1,2
              if (ix3==1) vec_cor(3)=xprobmin3
              if (ix3==2) vec_cor(3)=xprobmax3
              call get_cor_image(vec_cor,xI_cor)
              if (ix1==1 .and. ix2==1 .and. ix3==1) then
                xIFmin1=xI_cor(1)
                xIFmax1=xI_cor(1)
                xIFmin2=xI_cor(2)
                xIFmax2=xI_cor(2)
              else
                xIFmin1=min(xIFmin1,xI_cor(1))
                xIFmax1=max(xIFmax1,xI_cor(1))
                xIFmin2=min(xIFmin2,xI_cor(2))
                xIFmax2=max(xIFmax2,xI_cor(2))
              endif
            enddo
          enddo
        enddo
        dxDDA=min((xprobmax1-xprobmin1)/dble(numX1),&
                  (xprobmax2-xprobmin2)/dble(numX2),&
                  (xprobmax3-xprobmin3)/dble(numX3))
        xIcent1=half*(xIFmin1+xIFmax1)
        xIcent2=half*(xIFmin2+xIFmax2)
        nXIF1=max(1,ceiling((xIFmax1-xIFmin1)/dxDDA))
        nXIF2=max(1,ceiling((xIFmax2-xIFmin2)/dxDDA))
        xIFmin1=xIcent1-half*dble(nXIF1)*dxDDA
        xIFmax1=xIcent1+half*dble(nXIF1)*dxDDA
        xIFmin2=xIcent2-half*dble(nXIF2)*dxDDA
        xIFmax2=xIcent2+half*dble(nXIF2)*dxDDA
        bnx1=1
        bnx2=1
        nbb1=nXIF1
        nbb2=nXIF2
        strtype1=0
        strtype2=0
        nstrb1=0
        nstrb2=0
        qs1=one
        qs2=one
      else if (LOS_phi==0 .and. LOS_theta==90) then
        nXIF1=domain_nx2*2**(refine_max_level-1)
        nXIF2=domain_nx3*2**(refine_max_level-1)
        xIFmin1=xprobmin2
        xIFmax1=xprobmax2
        xIFmin2=xprobmin3
        xIFmax2=xprobmax3
        bnx1=block_nx2
        bnx2=block_nx3
        nbb1=domain_nx2
        nbb2=domain_nx3
        strtype1=stretch_type(2)
        strtype2=stretch_type(3)
        nstrb1=nstretchedblocks_baselevel(2)
        nstrb2=nstretchedblocks_baselevel(3)
        qs1=qstretch_baselevel(2)
        qs2=qstretch_baselevel(3)
        if (mype==0) write(*,'(a)') ' LOS vector: [-1.00  0.00  0.00]'
        if (mype==0) write(*,'(a)') ' xI1 vector: [ 0.00  1.00  0.00]'
        if (mype==0) write(*,'(a)') ' xI2 vector: [ 0.00  0.00  1.00]'
      else if (LOS_phi==90 .and. LOS_theta==90) then
        nXIF1=domain_nx3*2**(refine_max_level-1)
        nXIF2=domain_nx1*2**(refine_max_level-1)
        xIFmin1=xprobmin3
        xIFmax1=xprobmax3
        xIFmin2=xprobmin1
        xIFmax2=xprobmax1
        bnx1=block_nx3
        bnx2=block_nx1
        nbb1=domain_nx3
        nbb2=domain_nx1
        strtype1=stretch_type(3)
        strtype2=stretch_type(1)
        nstrb1=nstretchedblocks_baselevel(3)
        nstrb2=nstretchedblocks_baselevel(1)
        qs1=qstretch_baselevel(3)
        qs2=qstretch_baselevel(1)
        if (mype==0) write(*,'(a)') ' LOS vector: [ 0.00 -1.00  0.00]'
        if (mype==0) write(*,'(a)') ' xI1 vector: [-1.00  0.00  0.00]'
        if (mype==0) write(*,'(a)') ' xI2 vector: [ 0.00  0.00  1.00]'
      else
        nXIF1=domain_nx1*2**(refine_max_level-1)
        nXIF2=domain_nx2*2**(refine_max_level-1)
        xIFmin1=xprobmin1
        xIFmax1=xprobmax1
        xIFmin2=xprobmin2
        xIFmax2=xprobmax2
        bnx1=block_nx1
        bnx2=block_nx2
        nbb1=domain_nx1
        nbb2=domain_nx2
        strtype1=stretch_type(1)
        strtype2=stretch_type(2)
        nstrb1=nstretchedblocks_baselevel(1)
        nstrb2=nstretchedblocks_baselevel(2)
        qs1=qstretch_baselevel(1)
        qs2=qstretch_baselevel(2)
        if (mype==0) write(*,'(a)') ' LOS vector: [ 0.00  0.00 -1.00]'
        if (mype==0) write(*,'(a)') ' xI1 vector: [ 1.00  0.00  0.00]'
        if (mype==0) write(*,'(a)') ' xI2 vector: [ 0.00  1.00  0.00]'
      endif
      allocate(xIF1(nXIF1),xIF2(nXIF2),dxIF1(nXIF1),dxIF2(nXIF2))

      ! initialize image coordinate
      select case(strtype1)
      case(0) ! uniform
        dxIF1(:)=(xIFmax1-xIFmin1)/nXIF1
        do ix1=1,nXIF1
          xIF1(ix1)=xIFmin1+dxIF1(ix1)*(ix1-half)
        enddo
      case(1) ! uni stretch
        qs1=qs1**(one/2**(refine_max_level-1))
        dxfirst1=(xIFmax1-xIFmin1)*(one-qs1)/(one-qs1**nXIF1)
        dxIF1(1)=dxfirst1
        do ix1=2,nXIF1
          dxIF1(ix1)=dxfirst1*qs1**(ix1-1)
          xIF1(ix1)=dxIF1(1)/(one-qs1)*(one-qs1**(ix1-1))+half*dxIF1(ix1)
        enddo
      case(2) ! symm stretch
        ! base level, nbb = nstr + nuni + nstr
        nstr1=nstrb1*bnx1/2
        nuni1=nbb1-nstrb1*bnx1
        lenstr1=(xIFmax1-xIFmin1)/(2.d0+nuni1*(one-qs1)/(one-qs1**nstr1))
        dxfirst1=(xIFmax1-xIFmin1)/(dble(nuni1)+2.d0/(one-qs1)*(one-qs1**nstr1))
        dxmid1=dxfirst1
        ! refine_max level, numXI = nstr + nuni + nstr
        nstr1=nstr1*2**(refine_max_level-1)
        nuni1=nuni1*2**(refine_max_level-1)
        qs1=qs1**(one/2**(refine_max_level-1))
        dxfirst1=lenstr1*(one-qs1)/(one-qs1**nstr1)
        dxmid1=dxmid1/2**(refine_max_level-1)
        ! uniform center
        if(nuni1 .gt. 0) then
          do ix1=nstr1+1,nstr1+nuni1
            dxIF1(ix1)=dxmid1
            xIF1(ix1)=lenstr1+(dble(ix1)-0.5d0-nstr1)*dxIF1(ix1)+xIFmin1
          enddo
        endif
        ! left half
        do ix1=nstr1,1,-1
          dxIF1(ix1)=dxfirst1*qs1**(nstr1-ix1)
          xIF1(ix1)=xIFmin1+lenstr1-dxIF1(ix1)*half-dxfirst1*(one-qs1**(nstr1-ix1))/(one-qs1)
        enddo
        ! right half
        do ix1=nstr1+nuni1+1,nXIF1
          dxIF1(ix1)=dxfirst1*qs1**(ix1-nstr1-nuni1-1)
          xIF1(ix1)=xIFmax1-lenstr1+dxIF1(ix1)*half+dxfirst1*(one-qs1**(ix1-nstr1-nuni1-1))/(one-qs1)
        enddo
      case default
        call mpistop("unknown stretch type")
      end select

      select case(strtype2)
      case(0) ! uniform
        dxIF2(:)=(xIFmax2-xIFmin2)/nXIF2
        do ix2=1,nXIF2
          xIF2(ix2)=xIFmin2+dxIF2(ix2)*(ix2-half)
        enddo
      case(1) ! uni stretch
        qs2=qs2**(one/2**(refine_max_level-1))
        dxfirst2=(xIFmax2-xIFmin2)*(one-qs2)/(one-qs2**nXIF2)
        dxIF2(1)=dxfirst2
        do ix2=2,nXIF1
          dxIF2(ix2)=dxfirst2*qs2**(ix2-1)
          xIF2(ix2)=dxIF2(1)/(one-qs2)*(one-qs2**(ix2-1))+half*dxIF2(ix2)
        enddo
      case(2) ! symm stretch
        ! base level, nbb = nstr + nuni + nstr
        nstr2=nstrb2*bnx2/2
        nuni2=nbb2-nstrb2*bnx2
        lenstr2=(xIFmax2-xIFmin2)/(2.d0+nuni2*(one-qs2)/(one-qs2**nstr2))
        dxfirst2=(xIFmax2-xIFmin2)/(dble(nuni2)+2.d0/(one-qs2)*(one-qs2**nstr2))
        dxmid2=dxfirst2
        ! refine_max level, numXI = nstr + nuni + nstr
        nstr2=nstr2*2**(refine_max_level-1)
        nuni2=nuni2*2**(refine_max_level-1)
        qs2=qs2**(one/2**(refine_max_level-1))
        dxfirst2=lenstr2*(one-qs2)/(one-qs2**nstr2)
        dxmid2=dxmid2/2**(refine_max_level-1)
        ! uniform center
        if(nuni2 .gt. 0) then
          do ix2=nstr2+1,nstr2+nuni2
            dxIF2(ix2)=dxmid2
            xIF2(ix2)=lenstr2+(dble(ix2)-0.5d0-nstr2)*dxIF2(ix2)+xIFmin2
          enddo
        endif
        ! left half
        do ix2=nstr2,1,-1
          dxIF2(ix2)=dxfirst2*qs2**(nstr2-ix2)
          xIF2(ix2)=xIFmin2+lenstr2-dxIF2(ix2)*half-dxfirst2*(one-qs2**(nstr2-ix2))/(one-qs2)
        enddo
        ! right half
        do ix2=nstr2+nuni2+1,nXIF2
          dxIF2(ix2)=dxfirst2*qs2**(ix2-nstr2-nuni2-1)
          xIF2(ix2)=xIFmax2-lenstr2+dxIF2(ix2)*half+dxfirst2*(one-qs2**(ix2-nstr2-nuni2-1))/(one-qs2)
        enddo
      case default
        call mpistop("unknown stretch type")
      end select

      if (mype==0 .and. datatype=='image_euv') then
        if (SI_unit) then
          length_to_km=unit_length/1.d3
          arcsec=7.25d5/unit_length
        else
          length_to_km=unit_length/1.d5
          arcsec=7.25d7/unit_length
        endif
        write(*,'(a,i8,a,i8)') ' Native data-resolution image grid: ',nXIF1,' x ',nXIF2
        write(*,'(a,f10.3,a,f10.3,a,f8.3,a,f8.3,a)') &
          ' Native xI1 pixel-size range: ',minval(dxIF1)*length_to_km,'--', &
          maxval(dxIF1)*length_to_km,' km (',minval(dxIF1)/arcsec,'--',maxval(dxIF1)/arcsec,' arcsec)'
        write(*,'(a,f10.3,a,f10.3,a,f8.3,a,f8.3,a)') &
          ' Native xI2 pixel-size range: ',minval(dxIF2)*length_to_km,'--', &
          maxval(dxIF2)*length_to_km,' km (',minval(dxIF2)/arcsec,'--',maxval(dxIF2)/arcsec,' arcsec)'
      endif

      ! integrate EUV flux and get cell average flux for image
      if (datatype=='image_euv') then
        if (SI_unit) then
          unitv=unit_velocity/1.0e3 ! km/s
        else
          unitv=unit_velocity/1.0e5 ! km/s
        endif
        has_thick_output=trim(radiation_transfer)=='thick'
        has_doppler_output=radsyn_euv_has_doppler_output()
        numWI=radsyn_euv_num_outputs(has_doppler_output,has_thick_output)
        allocate(wI(nXIF1,nXIF2,numWI))
        allocate(EUV(nXIF1,nXIF2),Dpl(nXIF1,nXIF2))
        if (trim(radiation_transfer)=='thick') then
          allocate(EUVthin(nXIF1,nXIF2),Tau(nXIF1,nXIF2))
          if (coordinate==spherical .and. trim(ray_method_active)=='spherical') then
            call integrate_EUV_sph_intersection_thick(nXIF1,nXIF2,xIF1,xIF2,dxIF1(1),fl,EUV,Tau,EUVthin)
            Dpl=zero
          else if (trim(ray_method_active)=='cart') then
            call integrate_EUV_cart_dda_thick_datresol(nXIF1,nXIF2,xIF1,xIF2,fl,EUV,Dpl,Tau,EUVthin)
          else
            call integrate_EUV_thick_datresol(nXIF1,nXIF2,fl,EUV,Dpl,Tau,EUVthin)
          endif
          if (has_doppler_output) then
            where(EUV<smalldouble) EUV=zero
            call normalize_euv_doppler(nXIF1,nXIF2,EUV,Dpl,unitv)
          endif
        else
          allocate(EUVs(nXIF1,nXIF2),Dpls(nXIF1,nXIF2))
          EUVs=0.0d0
          EUV=0.0d0
          Dpl=0.d0
          Dpls=0.d0
          if (coordinate==spherical .and. trim(ray_method_active)=='spherical') then
            call integrate_EUV_sph_intersection_thin(nXIF1,nXIF2,xIF1,xIF2,dxIF1(1),fl,EUV)
            EUVs=EUV
            numSI=nXIF1*nXIF2
            call MPI_ALLREDUCE(EUVs,EUV,numSI,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,icomm,ierrmpi)
          else if (trim(ray_method_active)=='cart') then
            call integrate_EUV_cart_dda_datresol(nXIF1,nXIF2,xIF1,xIF2,fl,EUV,Dpl)
          else
            do iigrid=1,igridstail; igrid=igrids(iigrid);
              call integrate_EUV_datresol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,EUVs,Dpls)
            enddo
            numSI=nXIF1*nXIF2
            call MPI_ALLREDUCE(EUVs,EUV,numSI,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,icomm,ierrmpi)
            call MPI_ALLREDUCE(Dpls,Dpl,numSI,MPI_DOUBLE_PRECISION, &
                               MPI_SUM,icomm,ierrmpi)
          endif
          if (has_doppler_output) then
            where(EUV<smalldouble) EUV=zero
            call normalize_euv_doppler(nXIF1,nXIF2,EUV,Dpl,unitv)
          endif
          deallocate(EUVs,Dpls)
        endif
        if (has_thick_output) then
          if (has_doppler_output) then
            call pack_euv_image_outputs(nXIF1,nXIF2,EUV,wI,smalldouble,has_doppler_output,&
                                        has_thick_output,Dpl=Dpl,Tau=Tau,EUVthin=EUVthin)
          else
            call pack_euv_image_outputs(nXIF1,nXIF2,EUV,wI,smalldouble,has_doppler_output,&
                                        has_thick_output,Tau=Tau,EUVthin=EUVthin)
          endif
        else if (has_doppler_output) then
          call pack_euv_image_outputs(nXIF1,nXIF2,EUV,wI,smalldouble,has_doppler_output,&
                                      has_thick_output,Dpl=Dpl)
        else
          call pack_euv_image_outputs(nXIF1,nXIF2,EUV,wI,smalldouble,has_doppler_output,&
                                      has_thick_output)
        endif

        if (instrument_postprocess) then
          if (trim(emission_model)=='radio_ff') then
            if (trim(radiation_transfer)=='thick') then
              call postprocess_radio_beam_image(nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,&
                                                EUV,nXIP1,nXIP2,xIP1,xIP2,&
                                                dxIP1,dxIP2,wIP,numWIP,Tau=Tau,BrightThin=EUVthin)
            else
              call postprocess_radio_beam_image(nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,&
                                                EUV,nXIP1,nXIP2,xIP1,xIP2,&
                                                dxIP1,dxIP2,wIP,numWIP)
            endif
          else if (trim(radiation_transfer)=='thick') then
            call postprocess_euv_instrument_image(nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,&
                                                  EUV,Dpl,nXIP1,nXIP2,xIP1,xIP2,&
                                                  dxIP1,dxIP2,wIP,numWIP,Tau=Tau,EUVthin=EUVthin)
          else
            call postprocess_euv_instrument_image(nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,&
                                                  EUV,Dpl,nXIP1,nXIP2,xIP1,xIP2,&
                                                  dxIP1,dxIP2,wIP,numWIP)
          endif
          call output_data(qunit,xIP1,xIP2,dxIP1,dxIP2,wIP,nXIP1,nXIP2,numWIP,datatype)
          deallocate(xIP1,xIP2,dxIP1,dxIP2,wIP)
        else
          call output_data(qunit,xIF1,xIF2,dxIF1,dxIF2,wI,nXIF1,nXIF2,numWI,datatype)
        endif
        if (trim(radiation_transfer)=='thick') then
          deallocate(WI,EUV,Dpl,EUVthin,Tau)
        else
          deallocate(WI,EUV,Dpl)
        endif
      endif

      ! integrate SXR flux and get cell average flux for image
      if (datatype=='image_sxr') then
        if (SI_unit) then
          arcsec=7.25d5
        else
          arcsec=7.25d7
        endif
        RHESSI_rsl=2.3d0/instrument_resolution_factor
        numWI=1
        allocate(wI(nXIF1,nXIF2,numWI))
        allocate(SXRs(nXIF1,nXIF2),SXR(nXIF1,nXIF2))
        SXRs=0.0d0
        SXR=0.0d0
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call integrate_SXR_datresol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,SXRs)
        enddo
        numSI=nXIF1*nXIF2
        call MPI_ALLREDUCE(SXRs,SXR,numSI,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,icomm,ierrmpi)

        SXR=SXR*(RHESSI_rsl*arcsec)**2 ! photons cm^-2 s^-1 pixel^-1
        do ix1=1,nXIF1
          do ix2=1,nXIF2
            if (SXR(ix1,ix2)<smalldouble) SXR(ix1,ix2)=zero
          enddo
        enddo
        wI(:,:,1)=SXR(:,:)

        call output_data(qunit,xIF1,xIF2,dxIF1,dxIF2,wI,nXIF1,nXIF2,numWI,datatype)
        deallocate(WI,SXR,SXRs)
      endif

      deallocate(xIF1,xIF2,dxIF1,dxIF2)

    end subroutine get_image_datresol

    subroutine integrate_SXR_datresol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,SXR)
      use mod_global_parameters

      integer, intent(in) :: igrid,nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      double precision, intent(in) :: dxIF1(nXIF1),dxIF2(nXIF2)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: SXR(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&),opacity(:^D&)
      double precision, allocatable :: dxb1(:^D&),dxb2(:^D&),dxb3(:^D&)
      double precision, allocatable :: SXRg(:,:),xg1(:),xg2(:),dxg1(:),dxg2(:)
      integer :: levelg,nXg1,nXg2,iXgmin1,iXgmax1,iXgmin2,iXgmax2,rft,iXg^D
      double precision :: SXRt,xc^L,xg^L,r2,area_1AU
      integer :: ixP^L,ixP^D
      integer :: direction_LOS

      if (LOS_phi==0 .and. LOS_theta==90) then
        direction_LOS=1
      else if (LOS_phi==90 .and. LOS_theta==90) then
        direction_LOS=2
      else
        direction_LOS=3
      endif

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

      allocate(flux(ixI^S))
      allocate(dxb1(ixI^S),dxb2(ixI^S),dxb3(ixI^S))
      dxb1(ixO^S)=ps(igrid)%dx(ixO^S,1)
      dxb2(ixO^S)=ps(igrid)%dx(ixO^S,2)
      dxb3(ixO^S)=ps(igrid)%dx(ixO^S,3)
      ! get local SXR flux
      call get_SXR(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux,emin_sxr,emax_sxr)

      ! grid parameters
      levelg=ps(igrid)%level
      rft=2**(refine_max_level-levelg)

      ! fine table for storing EUV flux of current grid
      select case(direction_LOS)
      case(1)
        nXg1=ixImax2*rft
        nXg2=ixImax3*rft
      case(2)
        nXg1=ixImax3*rft
        nXg2=ixImax1*rft
      case(3)
        nXg1=ixImax1*rft
        nXg2=ixImax2*rft
      end select
      allocate(SXRg(nXg1,nXg2),xg1(nXg1),xg2(nXg2),dxg1(nXg1),dxg2(nXg2))
      SXRg=zero
      xg1=zero
      xg2=zero

      ! integrate for different direction
      select case(direction_LOS)
      case(1)
        do ix2=ixOmin2,ixOmax2
          iXgmin1=(ix2-1)*rft+1
          iXgmax1=ix2*rft
          do ix3=ixOmin3,ixOmax3
            iXgmin2=(ix3-1)*rft+1
            iXgmax2=ix3*rft
            SXRt=0.d0
            do ix1=ixOmin1,ixOmax1
              SXRt=SXRt+flux(ix^D)*dxb1(ix^D)*unit_length
            enddo
            SXRg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=SXRt
          enddo
        enddo
      case(2)
        do ix3=ixOmin3,ixOmax3
          iXgmin1=(ix3-1)*rft+1
          iXgmax1=ix3*rft
          do ix1=ixOmin1,ixOmax1
            iXgmin2=(ix1-1)*rft+1
            iXgmax2=ix1*rft
            SXRt=0.d0
            do ix2=ixOmin2,ixOmax2
              SXRt=SXRt+flux(ix^D)*dxb2(ix^D)*unit_length
            enddo
            SXRg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=SXRt
          enddo
        enddo
      case(3)
        do ix1=ixOmin1,ixOmax1
          iXgmin1=(ix1-1)*rft+1
          iXgmax1=ix1*rft
          do ix2=ixOmin2,ixOmax2
            iXgmin2=(ix2-1)*rft+1
            iXgmax2=ix2*rft
            SXRt=0.d0
            do ix3=ixOmin3,ixOmax3
              SXRt=SXRt+flux(ix^D)*dxb3(ix^D)*unit_length
            enddo
            SXRg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=SXRt
          enddo
        enddo
      end select

      area_1AU=2.81d27
      SXRg=SXRg/area_1AU

      ! mapping grid data to global table
      ! index ranges in local table
      select case(direction_LOS)
      case(1)
        iXgmin1=(ixOmin2-1)*rft+1
        iXgmax1=ixOmax2*rft
        iXgmin2=(ixOmin3-1)*rft+1
        iXgmax2=ixOmax3*rft
      case(2)
        iXgmin1=(ixOmin3-1)*rft+1
        iXgmax1=ixOmax3*rft
        iXgmin2=(ixOmin1-1)*rft+1
        iXgmax2=ixOmax1*rft
      case(3)
        iXgmin1=(ixOmin1-1)*rft+1
        iXgmax1=ixOmax1*rft
        iXgmin2=(ixOmin2-1)*rft+1
        iXgmax2=ixOmax2*rft
      end select
      ! index ranges in global table & mapping
      select case(direction_LOS)
      case(1)
        ixPmin1=(node(pig2_,igrid)-1)*rft*block_nx2+1
        ixPmax1=node(pig2_,igrid)*rft*block_nx2
        ixPmin2=(node(pig3_,igrid)-1)*rft*block_nx3+1
        ixPmax2=node(pig3_,igrid)*rft*block_nx3
      case(2)
        ixPmin1=(node(pig3_,igrid)-1)*rft*block_nx3+1
        ixPmax1=node(pig3_,igrid)*rft*block_nx3
        ixPmin2=(node(pig1_,igrid)-1)*rft*block_nx1+1
        ixPmax2=node(pig1_,igrid)*rft*block_nx1
      case(3)
        ixPmin1=(node(pig1_,igrid)-1)*rft*block_nx1+1
        ixPmax1=node(pig1_,igrid)*rft*block_nx1
        ixPmin2=(node(pig2_,igrid)-1)*rft*block_nx2+1
        ixPmax2=node(pig2_,igrid)*rft*block_nx2
      end select
      xg1(iXgmin1:iXgmax1)=xIF1(ixPmin1:ixPmax1)
      xg2(iXgmin2:iXgmax2)=xIF2(ixPmin2:ixPmax2)
      dxg1(iXgmin1:iXgmax1)=dxIF1(ixPmin1:ixPmax1)
      dxg2(iXgmin2:iXgmax2)=dxIF2(ixPmin2:ixPmax2)
      SXR(ixPmin1:ixPmax1,ixPmin2:ixPmax2)=SXR(ixPmin1:ixPmax1,ixPmin2:ixPmax2)+&
                                           SXRg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)

      deallocate(flux,dxb1,dxb2,dxb3,SXRg,xg1,xg2,dxg1,dxg2)

    end subroutine integrate_SXR_datresol

    subroutine integrate_EUV_datresol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,EUV,Dpl)
      use mod_global_parameters

      integer, intent(in) :: igrid,nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      double precision, intent(in) :: dxIF1(nXIF1),dxIF2(nXIF2)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: EUV(nXIF1,nXIF2),Dpl(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&),v(:^D&),rho(:^D&),opacity(:^D&)
      double precision, allocatable :: dxb1(:^D&),dxb2(:^D&),dxb3(:^D&)
      double precision, allocatable :: EUVg(:,:),Fvg(:,:),xg1(:),xg2(:),dxg1(:),dxg2(:)
      integer :: levelg,nXg1,nXg2,iXgmin1,iXgmax1,iXgmin2,iXgmax2,rft,iXg^D
      double precision :: EUVt,Fvt,xc^L,xg^L,r2
      integer :: ixP^L,ixP^D
      integer :: direction_LOS

      if (LOS_phi==0 .and. LOS_theta==90) then
        direction_LOS=1
      else if (LOS_phi==90 .and. LOS_theta==90) then
        direction_LOS=2
      else
        direction_LOS=3
      endif

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

      allocate(flux(ixI^S),v(ixI^S),rho(ixI^S),opacity(ixI^S))
      allocate(dxb1(ixI^S),dxb2(ixI^S),dxb3(ixI^S))
      dxb1(ixO^S)=ps(igrid)%dx(ixO^S,1)
      dxb2(ixO^S)=ps(igrid)%dx(ixO^S,2)
      dxb3(ixO^S)=ps(igrid)%dx(ixO^S,3)
      if (trim(emission_model)=='pseudo_current') then
        call get_pseudo_current(igrid,ixI^L,ixO^L,ps(igrid)%w,flux)
        v(ixO^S)=zero
        deallocate(rho)
      else if (trim(emission_model)=='radio_ff') then
        call get_radio_ff_source_opacity(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux,opacity)
        v(ixO^S)=zero
        deallocate(rho)
      else
        ! get local EUV flux and velocity
        call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
        flux(ixO^S)=flux(ixO^S)/instrument_resolution_factor**2   ! adjust flux due to artifical change of resolution
        call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
        v(ixO^S)=-ps(igrid)%w(ixO^S,iw_mom(direction_LOS))/rho(ixO^S)
        deallocate(rho)
      endif

      ! grid parameters
      levelg=ps(igrid)%level
      rft=2**(refine_max_level-levelg)

      ! fine table for storing EUV flux of current grid
      select case(direction_LOS)
      case(1)
        nXg1=ixImax2*rft
        nXg2=ixImax3*rft
      case(2)
        nXg1=ixImax3*rft
        nXg2=ixImax1*rft
      case(3)
        nXg1=ixImax1*rft
        nXg2=ixImax2*rft
      end select
      allocate(EUVg(nXg1,nXg2),Fvg(nXg1,nXg2),xg1(nXg1),xg2(nXg2),dxg1(nXg1),dxg2(nXg2))
      EUVg=zero
      Fvg=zero
      xg1=zero
      xg2=zero

      ! integrate for different direction
      select case(direction_LOS)
      case(1)
        do ix2=ixOmin2,ixOmax2
          iXgmin1=(ix2-1)*rft+1
          iXgmax1=ix2*rft
          do ix3=ixOmin3,ixOmax3
            iXgmin2=(ix3-1)*rft+1
            iXgmax2=ix3*rft
            EUVt=0.d0
            Fvt=0.d0
            do ix1=ixOmin1,ixOmax1
              EUVt=EUVt+flux(ix^D)*dxb1(ix^D)*unit_length
              Fvt=Fvt+flux(ix^D)*dxb1(ix^D)*unit_length*v(ix^D)
            enddo
            EUVg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=EUVt
            Fvg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=Fvt
          enddo
        enddo
      case(2)
        do ix3=ixOmin3,ixOmax3
          iXgmin1=(ix3-1)*rft+1
          iXgmax1=ix3*rft
          do ix1=ixOmin1,ixOmax1
            iXgmin2=(ix1-1)*rft+1
            iXgmax2=ix1*rft
            EUVt=0.d0
            Fvt=0.d0
            do ix2=ixOmin2,ixOmax2
              EUVt=EUVt+flux(ix^D)*dxb2(ix^D)*unit_length
              Fvt=Fvt+flux(ix^D)*dxb2(ix^D)*unit_length*v(ix^D)
            enddo
            EUVg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=EUVt
            Fvg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=Fvt
          enddo
        enddo
      case(3)
        do ix1=ixOmin1,ixOmax1
          iXgmin1=(ix1-1)*rft+1
          iXgmax1=ix1*rft
          do ix2=ixOmin2,ixOmax2
            iXgmin2=(ix2-1)*rft+1
            iXgmax2=ix2*rft
            EUVt=0.d0
            Fvt=0.d0
            do ix3=ixOmin3,ixOmax3
              EUVt=EUVt+flux(ix^D)*dxb3(ix^D)*unit_length
              Fvt=Fvt+flux(ix^D)*dxb3(ix^D)*unit_length*v(ix^D)
            enddo
            EUVg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=EUVt
            Fvg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)=Fvt
          enddo
        enddo
      end select
      if (SI_unit) then
        EUVg=EUVg*1.d2
        Fvg=Fvg*1.d2
      endif

      ! mapping grid data to global table
      ! index ranges in local table
      select case(direction_LOS)
      case(1)
        iXgmin1=(ixOmin2-1)*rft+1
        iXgmax1=ixOmax2*rft
        iXgmin2=(ixOmin3-1)*rft+1
        iXgmax2=ixOmax3*rft
      case(2)
        iXgmin1=(ixOmin3-1)*rft+1
        iXgmax1=ixOmax3*rft
        iXgmin2=(ixOmin1-1)*rft+1
        iXgmax2=ixOmax1*rft
      case(3)
        iXgmin1=(ixOmin1-1)*rft+1
        iXgmax1=ixOmax1*rft
        iXgmin2=(ixOmin2-1)*rft+1
        iXgmax2=ixOmax2*rft
      end select
      ! index ranges in global table & mapping
      select case(direction_LOS)
      case(1)
        ixPmin1=(node(pig2_,igrid)-1)*rft*block_nx2+1
        ixPmax1=node(pig2_,igrid)*rft*block_nx2
        ixPmin2=(node(pig3_,igrid)-1)*rft*block_nx3+1
        ixPmax2=node(pig3_,igrid)*rft*block_nx3
      case(2)
        ixPmin1=(node(pig3_,igrid)-1)*rft*block_nx3+1
        ixPmax1=node(pig3_,igrid)*rft*block_nx3
        ixPmin2=(node(pig1_,igrid)-1)*rft*block_nx1+1
        ixPmax2=node(pig1_,igrid)*rft*block_nx1
      case(3)
        ixPmin1=(node(pig1_,igrid)-1)*rft*block_nx1+1
        ixPmax1=node(pig1_,igrid)*rft*block_nx1
        ixPmin2=(node(pig2_,igrid)-1)*rft*block_nx2+1
        ixPmax2=node(pig2_,igrid)*rft*block_nx2
      end select
      xg1(iXgmin1:iXgmax1)=xIF1(ixPmin1:ixPmax1)
      xg2(iXgmin2:iXgmax2)=xIF2(ixPmin2:ixPmax2)
      dxg1(iXgmin1:iXgmax1)=dxIF1(ixPmin1:ixPmax1)
      dxg2(iXgmin2:iXgmax2)=dxIF2(ixPmin2:ixPmax2)
      EUV(ixPmin1:ixPmax1,ixPmin2:ixPmax2)=EUV(ixPmin1:ixPmax1,ixPmin2:ixPmax2)+&
                                           EUVg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)
      Dpl(ixPmin1:ixPmax1,ixPmin2:ixPmax2)=Dpl(ixPmin1:ixPmax1,ixPmin2:ixPmax2)+&
                                           FVg(iXgmin1:iXgmax1,iXgmin2:iXgmax2)

      deallocate(flux,v,opacity,dxb1,dxb2,dxb3,EUVg,Fvg,xg1,xg2,dxg1,dxg2)

    end subroutine integrate_EUV_datresol

  }

  {^IFTHREED

    subroutine ray_box_intersection_cart(ray_origin,ray_dir,box_min,box_max,hit,t_enter,t_exit)
      double precision, intent(in) :: ray_origin(1:3),ray_dir(1:3),box_min(1:3),box_max(1:3)
      logical, intent(out) :: hit
      double precision, intent(out) :: t_enter,t_exit

      integer :: idir
      double precision :: t1,t2,td

      hit=.true.
      t_enter=-huge(one)
      t_exit=huge(one)
      do idir=1,3
        if (abs(ray_dir(idir))<=smalldouble) then
          if (ray_origin(idir)<box_min(idir) .or. ray_origin(idir)>box_max(idir)) then
            hit=.false.
            return
          endif
        else
          t1=(box_min(idir)-ray_origin(idir))/ray_dir(idir)
          t2=(box_max(idir)-ray_origin(idir))/ray_dir(idir)
          if (t1>t2) then
            td=t1
            t1=t2
            t2=td
          endif
          t_enter=max(t_enter,t1)
          t_exit=min(t_exit,t2)
          if (t_enter>=t_exit) then
            hit=.false.
            return
          endif
        endif
      enddo
    end subroutine ray_box_intersection_cart

    subroutine build_cart_dda_faces(ixI^L,ixO^L,x,dx,xface1,xface2,xface3)
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim),dx(ixI^S,1:ndim)
      double precision, allocatable, intent(out) :: xface1(:),xface2(:),xface3(:)

      integer :: ix^D

      allocate(xface1(ixOmin1:ixOmax1+1),xface2(ixOmin2:ixOmax2+1),xface3(ixOmin3:ixOmax3+1))

      ix2=ixOmin2
      ix3=ixOmin3
      do ix1=ixOmin1,ixOmax1
        xface1(ix1)=x(ix^D,1)-half*dx(ix^D,1)
      enddo
      xface1(ixOmax1+1)=x(ixOmax1,ixOmin2,ixOmin3,1)+half*dx(ixOmax1,ixOmin2,ixOmin3,1)

      ix1=ixOmin1
      ix3=ixOmin3
      do ix2=ixOmin2,ixOmax2
        xface2(ix2)=x(ix^D,2)-half*dx(ix^D,2)
      enddo
      xface2(ixOmax2+1)=x(ixOmin1,ixOmax2,ixOmin3,2)+half*dx(ixOmin1,ixOmax2,ixOmin3,2)

      ix1=ixOmin1
      ix2=ixOmin2
      do ix3=ixOmin3,ixOmax3
        xface3(ix3)=x(ix^D,3)-half*dx(ix^D,3)
      enddo
      xface3(ixOmax3+1)=x(ixOmin1,ixOmin2,ixOmax3,3)+half*dx(ixOmin1,ixOmin2,ixOmax3,3)
    end subroutine build_cart_dda_faces

    integer function cart_dda_locate_index(pos,faces,imin,imax) result(idx)
      integer, intent(in) :: imin,imax
      double precision, intent(in) :: pos,faces(imin:imax+1)

      integer :: ilo,ihi,imid

      if (pos<=faces(imin)) then
        idx=imin
        return
      endif
      if (pos>=faces(imax+1)) then
        idx=imax
        return
      endif

      ilo=imin
      ihi=imax+1
      do while (ihi-ilo>1)
        imid=(ilo+ihi)/2
        if (pos>=faces(imid)) then
          ilo=imid
        else
          ihi=imid
        endif
      enddo
      idx=min(imax,max(imin,ilo))
    end function cart_dda_locate_index

    subroutine cart_dda_init_axis(ray_origin_axis,ray_dir_axis,faces,imin,imax,idx,step,tMax)
      integer, intent(in) :: imin,imax,idx
      double precision, intent(in) :: ray_origin_axis,ray_dir_axis,faces(imin:imax+1)
      integer, intent(out) :: step
      double precision, intent(out) :: tMax

      if (ray_dir_axis>zero) then
        step=1
        tMax=(faces(idx+1)-ray_origin_axis)/ray_dir_axis
      else if (ray_dir_axis<zero) then
        step=-1
        tMax=(faces(idx)-ray_origin_axis)/ray_dir_axis
      else
        step=0
        tMax=huge(one)
      endif
    end subroutine cart_dda_init_axis

    subroutine cart_dda_advance_axis(ray_origin_axis,ray_dir_axis,faces,imin,imax,idx,step,tMax,done)
      integer, intent(in) :: imin,imax,step
      double precision, intent(in) :: ray_origin_axis,ray_dir_axis,faces(imin:imax+1)
      integer, intent(inout) :: idx
      double precision, intent(inout) :: tMax
      logical, intent(out) :: done

      done=.false.
      idx=idx+step
      if (idx<imin .or. idx>imax) then
        done=.true.
        return
      endif
      if (step>0) then
        tMax=(faces(idx+1)-ray_origin_axis)/ray_dir_axis
      else if (step<0) then
        tMax=(faces(idx)-ray_origin_axis)/ray_dir_axis
      else
        tMax=huge(one)
      endif
    end subroutine cart_dda_advance_axis

    subroutine acc_EUV_cart_dda(ixI^L,ixO^L,source,sourcev,&
                                ray_origin,xface1,xface2,xface3,t_enter,t_exit,EUVp,Dplp)
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: source(ixI^S),sourcev(ixI^S)
      double precision, intent(in) :: ray_origin(1:3)
      double precision, intent(in) :: xface1(ixOmin1:ixOmax1+1),xface2(ixOmin2:ixOmax2+1),&
                                      xface3(ixOmin3:ixOmax3+1)
      double precision, intent(in) :: t_enter,t_exit
      double precision, intent(inout) :: EUVp,Dplp

      integer :: ix^D,step(1:3)
      double precision :: pos(1:3),tMax(1:3),tNow,tNext,ds_cm,epsRay
      logical :: done

      if (t_exit<=t_enter) return
      epsRay=max(1.d-12,1.d-10*abs(t_exit-t_enter))
      pos=ray_origin+(t_enter+epsRay)*vec_LOS
      ix1=cart_dda_locate_index(pos(1),xface1,ixOmin1,ixOmax1)
      ix2=cart_dda_locate_index(pos(2),xface2,ixOmin2,ixOmax2)
      ix3=cart_dda_locate_index(pos(3),xface3,ixOmin3,ixOmax3)

      call cart_dda_init_axis(ray_origin(1),vec_LOS(1),xface1,ixOmin1,ixOmax1,ix1,step(1),tMax(1))
      call cart_dda_init_axis(ray_origin(2),vec_LOS(2),xface2,ixOmin2,ixOmax2,ix2,step(2),tMax(2))
      call cart_dda_init_axis(ray_origin(3),vec_LOS(3),xface3,ixOmin3,ixOmax3,ix3,step(3),tMax(3))

      tNow=t_enter
      do
        tNext=min(t_exit,tMax(1),tMax(2),tMax(3))
        if (tNext>tNow) then
          ds_cm=(tNext-tNow)*unit_length
          if (SI_unit) ds_cm=ds_cm*1.d2
          EUVp=EUVp+source(ix^D)*ds_cm
          Dplp=Dplp+sourcev(ix^D)*ds_cm
        endif
        tNow=tNext
        if (tNow>=t_exit-epsRay) exit

        if (tMax(1)<=tNow+epsRay) then
          call cart_dda_advance_axis(ray_origin(1),vec_LOS(1),xface1,ixOmin1,ixOmax1,ix1,step(1),tMax(1),done)
          if (done) exit
        endif
        if (tMax(2)<=tNow+epsRay) then
          call cart_dda_advance_axis(ray_origin(2),vec_LOS(2),xface2,ixOmin2,ixOmax2,ix2,step(2),tMax(2),done)
          if (done) exit
        endif
        if (tMax(3)<=tNow+epsRay) then
          call cart_dda_advance_axis(ray_origin(3),vec_LOS(3),xface3,ixOmin3,ixOmax3,ix3,step(3),tMax(3),done)
          if (done) exit
        endif
      enddo
    end subroutine acc_EUV_cart_dda

    subroutine append_cart_dda_segment(segments,nseg,capacity,pixel_id,tseg,jds,kds,jvds)
      double precision, allocatable, intent(inout) :: segments(:,:)
      integer, intent(inout) :: nseg,capacity
      integer, intent(in) :: pixel_id
      double precision, intent(in) :: tseg,jds,kds,jvds

      double precision, allocatable :: tmp(:,:)
      integer :: new_capacity

      if (capacity<=0) then
        capacity=1024
        allocate(segments(5,capacity))
      else if (nseg>=capacity) then
        new_capacity=2*capacity
        allocate(tmp(5,new_capacity))
        tmp(:,1:capacity)=segments(:,1:capacity)
        call move_alloc(tmp,segments)
        capacity=new_capacity
      endif

      nseg=nseg+1
      segments(1,nseg)=dble(pixel_id)
      segments(2,nseg)=tseg
      segments(3,nseg)=jds
      segments(4,nseg)=kds
      segments(5,nseg)=jvds
    end subroutine append_cart_dda_segment

    subroutine collect_EUV_cart_dda_segments(ixI^L,ixO^L,source,opacity,sourcev,&
                                             pixel_id,ray_origin,xface1,xface2,xface3,t_enter,t_exit,&
                                             segments,nseg,capacity)
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: source(ixI^S),opacity(ixI^S),sourcev(ixI^S)
      integer, intent(in) :: pixel_id
      double precision, intent(in) :: ray_origin(1:3)
      double precision, intent(in) :: xface1(ixOmin1:ixOmax1+1),xface2(ixOmin2:ixOmax2+1),&
                                      xface3(ixOmin3:ixOmax3+1)
      double precision, intent(in) :: t_enter,t_exit
      double precision, allocatable, intent(inout) :: segments(:,:)
      integer, intent(inout) :: nseg,capacity

      integer :: ix^D,step(1:3)
      double precision :: pos(1:3),tMax(1:3),tNow,tNext,ds_cm,epsRay,tseg
      double precision :: jds,kds,jvds
      logical :: done

      if (t_exit<=t_enter) return
      epsRay=max(1.d-12,1.d-10*abs(t_exit-t_enter))
      pos=ray_origin+(t_enter+epsRay)*vec_LOS
      ix1=cart_dda_locate_index(pos(1),xface1,ixOmin1,ixOmax1)
      ix2=cart_dda_locate_index(pos(2),xface2,ixOmin2,ixOmax2)
      ix3=cart_dda_locate_index(pos(3),xface3,ixOmin3,ixOmax3)

      call cart_dda_init_axis(ray_origin(1),vec_LOS(1),xface1,ixOmin1,ixOmax1,ix1,step(1),tMax(1))
      call cart_dda_init_axis(ray_origin(2),vec_LOS(2),xface2,ixOmin2,ixOmax2,ix2,step(2),tMax(2))
      call cart_dda_init_axis(ray_origin(3),vec_LOS(3),xface3,ixOmin3,ixOmax3,ix3,step(3),tMax(3))

      tNow=t_enter
      do
        tNext=min(t_exit,tMax(1),tMax(2),tMax(3))
        if (tNext>tNow) then
          ds_cm=(tNext-tNow)*unit_length
          if (SI_unit) ds_cm=ds_cm*1.d2
          jds=source(ix^D)*ds_cm
          kds=opacity(ix^D)*ds_cm
          jvds=sourcev(ix^D)*ds_cm
          if (jds/=zero .or. kds/=zero .or. jvds/=zero) then
            tseg=half*(tNow+tNext)
            call append_cart_dda_segment(segments,nseg,capacity,pixel_id,tseg,jds,kds,jvds)
          endif
        endif
        tNow=tNext
        if (tNow>=t_exit-epsRay) exit

        if (tMax(1)<=tNow+epsRay) then
          call cart_dda_advance_axis(ray_origin(1),vec_LOS(1),xface1,ixOmin1,ixOmax1,ix1,step(1),tMax(1),done)
          if (done) exit
        endif
        if (tMax(2)<=tNow+epsRay) then
          call cart_dda_advance_axis(ray_origin(2),vec_LOS(2),xface2,ixOmin2,ixOmax2,ix2,step(2),tMax(2),done)
          if (done) exit
        endif
        if (tMax(3)<=tNow+epsRay) then
          call cart_dda_advance_axis(ray_origin(3),vec_LOS(3),xface3,ixOmin3,ixOmax3,ix3,step(3),tMax(3),done)
          if (done) exit
        endif
      enddo
    end subroutine collect_EUV_cart_dda_segments

    subroutine sort_segment_indices_near_to_far(segments,idx,nidx)
      double precision, intent(in) :: segments(:,:)
      integer, intent(inout) :: idx(:)
      integer, intent(in) :: nidx

      if (nidx<=1) return
      if (nidx<=32) then
        call insertion_sort_segment_indices(segments,idx,1,nidx)
      else
        call quicksort_segment_indices(segments,idx,1,nidx)
      endif
    end subroutine sort_segment_indices_near_to_far

    subroutine insertion_sort_segment_indices(segments,idx,ilo,ihi)
      double precision, intent(in) :: segments(:,:)
      integer, intent(inout) :: idx(:)
      integer, intent(in) :: ilo,ihi

      integer :: i,j,key

      do i=ilo+1,ihi
        key=idx(i)
        j=i-1
        do while (j>=ilo .and. segments(2,idx(j))>segments(2,key))
          idx(j+1)=idx(j)
          j=j-1
        enddo
        idx(j+1)=key
      enddo
    end subroutine insertion_sort_segment_indices

    recursive subroutine quicksort_segment_indices(segments,idx,ilo,ihi)
      double precision, intent(in) :: segments(:,:)
      integer, intent(inout) :: idx(:)
      integer, intent(in) :: ilo,ihi

      integer :: i,j,tmp
      double precision :: pivot

      if (ihi-ilo<=32) then
        call insertion_sort_segment_indices(segments,idx,ilo,ihi)
        return
      endif

      i=ilo
      j=ihi
      pivot=segments(2,idx((ilo+ihi)/2))
      do
        do while (segments(2,idx(i))<pivot)
          i=i+1
        enddo
        do while (segments(2,idx(j))>pivot)
          j=j-1
        enddo
        if (i<=j) then
          tmp=idx(i)
          idx(i)=idx(j)
          idx(j)=tmp
          i=i+1
          j=j-1
        endif
        if (i>j) exit
      enddo

      if (ilo<j) call quicksort_segment_indices(segments,idx,ilo,j)
      if (i<ihi) call quicksort_segment_indices(segments,idx,i,ihi)
    end subroutine quicksort_segment_indices

    integer function segment_pixel_owner(pixel_id) result(owner)
      integer, intent(in) :: pixel_id

      owner=mod(pixel_id-1,npe)
    end function segment_pixel_owner

    logical function segment_is_valid(segments,is,nvars) result(valid)
      double precision, intent(in) :: segments(:,:)
      integer, intent(in) :: is,nvars

      integer :: iv

      valid=.true.
      do iv=1,nvars
        if (segments(iv,is)/=segments(iv,is) .or. abs(segments(iv,is))>=1.d90) then
          valid=.false.
          return
        endif
      enddo
    end function segment_is_valid

    subroutine cart_dda_block_pixel_range(box_min,box_max,nXIF1,nXIF2,xIF1,xIF2,ixPmin1,ixPmax1,ixPmin2,ixPmax2,has_pixels)
      double precision, intent(in) :: box_min(1:3),box_max(1:3)
      integer, intent(in) :: nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      integer, intent(out) :: ixPmin1,ixPmax1,ixPmin2,ixPmax2
      logical, intent(out) :: has_pixels

      integer :: i1,i2,i3
      double precision :: vec_cor(1:3),xI_cor(1:2)
      double precision :: xmin1,xmax1,xmin2,xmax2,dx1,dx2

      do i1=1,2
        if (i1==1) vec_cor(1)=box_min(1)
        if (i1==2) vec_cor(1)=box_max(1)
        do i2=1,2
          if (i2==1) vec_cor(2)=box_min(2)
          if (i2==2) vec_cor(2)=box_max(2)
          do i3=1,2
            if (i3==1) vec_cor(3)=box_min(3)
            if (i3==2) vec_cor(3)=box_max(3)
            call get_cor_image(vec_cor,xI_cor)
            if (i1==1 .and. i2==1 .and. i3==1) then
              xmin1=xI_cor(1)
              xmax1=xI_cor(1)
              xmin2=xI_cor(2)
              xmax2=xI_cor(2)
            else
              xmin1=min(xmin1,xI_cor(1))
              xmax1=max(xmax1,xI_cor(1))
              xmin2=min(xmin2,xI_cor(2))
              xmax2=max(xmax2,xI_cor(2))
            endif
          enddo
        enddo
      enddo

      if (nXIF1>1) then
        dx1=abs(xIF1(2)-xIF1(1))
      else
        dx1=max(one,abs(xmax1-xmin1))
      endif
      if (nXIF2>1) then
        dx2=abs(xIF2(2)-xIF2(1))
      else
        dx2=max(one,abs(xmax2-xmin2))
      endif

      ixPmin1=max(1,floor((xmin1-xIF1(1))/dx1)+1-1)
      ixPmax1=min(nXIF1,ceiling((xmax1-xIF1(1))/dx1)+1+1)
      ixPmin2=max(1,floor((xmin2-xIF2(1))/dx2)+1-1)
      ixPmax2=min(nXIF2,ceiling((xmax2-xIF2(1))/dx2)+1+1)
      has_pixels=ixPmin1<=ixPmax1 .and. ixPmin2<=ixPmax2
    end subroutine cart_dda_block_pixel_range

    subroutine integrate_EUV_cart_dda_datresol(nXIF1,nXIF2,xIF1,xIF2,fl,EUV,Dpl)
      use mod_global_parameters

      integer, intent(in) :: nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: EUV(nXIF1,nXIF2),Dpl(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D
      integer :: iigrid,igrid,ixP1,ixP2,numSI,ixPmin1,ixPmax1,ixPmin2,ixPmax2
      double precision :: box_min(1:3),box_max(1:3),ray_origin(1:3)
      double precision :: t_enter,t_exit,vlos
      double precision :: profile_local(2),profile_global(2)
      logical :: hit,has_pixels
      double precision, allocatable :: source(:^D&),sourcev(:^D&),rho(:^D&),opacity(:^D&)
      double precision, allocatable :: xface1(:),xface2(:),xface3(:)
      double precision, allocatable :: EUVs(:,:),Dpls(:,:)

      allocate(EUVs(nXIF1,nXIF2),Dpls(nXIF1,nXIF2))
      EUVs=zero
      Dpls=zero
      profile_local=zero

      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixOmax^D=ixmhi^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\

        box_min(1)=rnode(rpxmin1_,igrid)
        box_min(2)=rnode(rpxmin2_,igrid)
        box_min(3)=rnode(rpxmin3_,igrid)
        box_max(1)=rnode(rpxmax1_,igrid)
        box_max(2)=rnode(rpxmax2_,igrid)
        box_max(3)=rnode(rpxmax3_,igrid)
        call build_cart_dda_faces(ixI^L,ixO^L,ps(igrid)%x,ps(igrid)%dx,xface1,xface2,xface3)
        call cart_dda_block_pixel_range(box_min,box_max,nXIF1,nXIF2,xIF1,xIF2,&
                                        ixPmin1,ixPmax1,ixPmin2,ixPmax2,has_pixels)
        if (.not. has_pixels) then
          deallocate(xface1,xface2,xface3)
          cycle
        endif

        allocate(source(ixI^S),sourcev(ixI^S),rho(ixI^S),opacity(ixI^S))
        source=zero
        sourcev=zero
        if (trim(emission_model)=='pseudo_current') then
          call get_pseudo_current(igrid,ixI^L,ixO^L,ps(igrid)%w,source)
        else if (trim(emission_model)=='radio_ff') then
          call get_radio_ff_source_opacity(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,source,opacity)
        else
          call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,source)
          source(ixO^S)=source(ixO^S)/instrument_resolution_factor**2
          call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
          do ix1=ixOmin1,ixOmax1
            do ix2=ixOmin2,ixOmax2
              do ix3=ixOmin3,ixOmax3
                if (rho(ix^D)>smalldouble) then
                  vlos=(ps(igrid)%w(ix^D,iw_mom(1))*vec_LOS(1)+&
                        ps(igrid)%w(ix^D,iw_mom(2))*vec_LOS(2)+&
                        ps(igrid)%w(ix^D,iw_mom(3))*vec_LOS(3))/rho(ix^D)
                  sourcev(ix^D)=source(ix^D)*vlos
                endif
              enddo
            enddo
          enddo
        endif
        deallocate(rho,opacity)

        do ixP1=ixPmin1,ixPmax1
          do ixP2=ixPmin2,ixPmax2
            ray_origin=x_origin+xIF1(ixP1)*vec_xI1+xIF2(ixP2)*vec_xI2
            profile_local(1)=profile_local(1)+one
            call ray_box_intersection_cart(ray_origin,vec_LOS,box_min,box_max,hit,t_enter,t_exit)
            if (hit) then
              profile_local(2)=profile_local(2)+one
              call acc_EUV_cart_dda(ixI^L,ixO^L,source,sourcev,&
                                    ray_origin,xface1,xface2,xface3,t_enter,t_exit,EUVs(ixP1,ixP2),Dpls(ixP1,ixP2))
            endif
          enddo
        enddo

        deallocate(source,sourcev,xface1,xface2,xface3)
      enddo

      numSI=nXIF1*nXIF2
      call MPI_ALLREDUCE(EUVs,EUV,numSI,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      call MPI_ALLREDUCE(Dpls,Dpl,numSI,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      call MPI_ALLREDUCE(profile_local,profile_global,2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      if (radsyn_verbose .and. mype==0) then
        write(*,'(a,2(es12.5,1x))') ' cart_dda thin profile ray_tests ray_hits: ',profile_global
      endif
      deallocate(EUVs,Dpls)
    end subroutine integrate_EUV_cart_dda_datresol

    subroutine integrate_EUV_cart_dda_thick_datresol(nXIF1,nXIF2,xIF1,xIF2,fl,EUV,Dpl,Tau,EUVthin)
      use mod_global_parameters

      integer, intent(in) :: nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: EUV(nXIF1,nXIF2),Dpl(nXIF1,nXIF2)
      double precision, intent(out) :: Tau(nXIF1,nXIF2),EUVthin(nXIF1,nXIF2)

      integer, parameter :: nSegVars=5
      integer :: ixO^L,ixO^D,ixI^L,ix^D
      integer :: iigrid,igrid,ixP1,ixP2,ipix,ipixStart,ipixEnd,nPixBatch,pixel_id
      integer :: nseg,capacity,totalCount,totalSeg,ipe,is,iseg,nidx,owner,isegDest,nsegBefore
      integer :: ixGlobal,iyGlobal,ixPmin1,ixPmax1,ixPmin2,ixPmax2,iFirst,iLast,iLocal
      integer :: nPixBatchTarget
      integer :: maxSegBatchTarget,maxSegCommTarget,maxNsegBatch,nPixTotal
      integer :: maxOwnerSegCount,maxOwnerSegCountLocal,segOffset,recvFill,totalRoundCount,totalRoundSeg
      integer, allocatable :: sendCounts(:),recvCounts(:),sendDispls(:),recvDispls(:)
      integer, allocatable :: roundSendCounts(:),roundRecvCounts(:)
      integer, allocatable :: roundSendDispls(:),roundRecvDispls(:)
      integer, allocatable :: ownerSegCounts(:),ownerOffsets(:),idx(:)
      integer, allocatable :: bucketCounts(:),bucketOffsets(:),bucketFill(:)
      double precision :: ray_origin(1:3)
      double precision :: t_enter,t_exit,vlos,atten
      double precision :: profile_local(5),profile_global(5),profile_batch(5)
      logical :: hit,has_pixels,batchAccepted,batchReduced
      double precision, allocatable :: rho(:^D&)
      double precision, allocatable :: segments(:,:),segments_send(:,:),segments_recv(:,:)
      double precision, allocatable :: segments_recv_round(:,:)
      double precision, allocatable :: image_reduce(:,:)
      type(radsyn_euv_cache), allocatable :: cache(:)

      EUV=zero
      Dpl=zero
      Tau=zero
      EUVthin=zero
      profile_local=zero
      allocate(sendCounts(0:npe-1),recvCounts(0:npe-1),sendDispls(0:npe-1),recvDispls(0:npe-1))
      allocate(roundSendCounts(0:npe-1),roundRecvCounts(0:npe-1))
      allocate(roundSendDispls(0:npe-1),roundRecvDispls(0:npe-1))
      allocate(ownerSegCounts(0:npe-1),ownerOffsets(0:npe-1))
      allocate(cache(igridstail))
      call radsyn_get_segment_batch_limits(nPixBatchTarget,maxSegBatchTarget,maxSegCommTarget)
      allocate(bucketCounts(nPixBatchTarget),bucketOffsets(nPixBatchTarget+1),&
               bucketFill(nPixBatchTarget))

      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixOmax^D=ixmhi^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\

        cache(iigrid)%igrid=igrid
        allocate(cache(iigrid)%source(ixI^S),cache(iigrid)%opacity(ixI^S),&
                 cache(iigrid)%sourcev(ixI^S),rho(ixI^S))
        cache(iigrid)%source=zero
        cache(iigrid)%opacity=zero
        cache(iigrid)%sourcev=zero
        if (trim(emission_model)=='radio_ff') then
          call get_radio_ff_source_opacity(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,&
                                           cache(iigrid)%source,cache(iigrid)%opacity)
        else
          call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,cache(iigrid)%source)
          cache(iigrid)%source(ixO^S)=cache(iigrid)%source(ixO^S)/instrument_resolution_factor**2
          call get_EUV_HHe_opacity(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,cache(iigrid)%opacity)
          call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
          do ix1=ixOmin1,ixOmax1
            do ix2=ixOmin2,ixOmax2
              do ix3=ixOmin3,ixOmax3
                if (rho(ix^D)>smalldouble) then
                  vlos=(ps(igrid)%w(ix^D,iw_mom(1))*vec_LOS(1)+&
                        ps(igrid)%w(ix^D,iw_mom(2))*vec_LOS(2)+&
                        ps(igrid)%w(ix^D,iw_mom(3))*vec_LOS(3))/rho(ix^D)
                  cache(iigrid)%sourcev(ix^D)=cache(iigrid)%source(ix^D)*vlos
                endif
              enddo
            enddo
          enddo
        endif
        deallocate(rho)
        cache(iigrid)%box_min(1)=rnode(rpxmin1_,igrid)
        cache(iigrid)%box_min(2)=rnode(rpxmin2_,igrid)
        cache(iigrid)%box_min(3)=rnode(rpxmin3_,igrid)
        cache(iigrid)%box_max(1)=rnode(rpxmax1_,igrid)
        cache(iigrid)%box_max(2)=rnode(rpxmax2_,igrid)
        cache(iigrid)%box_max(3)=rnode(rpxmax3_,igrid)
        call build_cart_dda_faces(ixI^L,ixO^L,ps(igrid)%x,ps(igrid)%dx,&
                                  cache(iigrid)%xface1,cache(iigrid)%xface2,&
                                  cache(iigrid)%xface3)
        call cart_dda_block_pixel_range(cache(iigrid)%box_min,cache(iigrid)%box_max,&
             nXIF1,nXIF2,xIF1,xIF2,cache(iigrid)%ixPmin1,cache(iigrid)%ixPmax1,&
             cache(iigrid)%ixPmin2,cache(iigrid)%ixPmax2,cache(iigrid)%has_pixels)
      enddo

      nPixTotal=nXIF1*nXIF2
      ipixStart=1
      do while (ipixStart<=nPixTotal)
        ipixEnd=min(nXIF1*nXIF2,ipixStart+nPixBatchTarget-1)
        nPixBatch=ipixEnd-ipixStart+1
        batchAccepted=.false.
        batchReduced=.false.

        do while (.not. batchAccepted)
          nseg=0
          capacity=0
          profile_batch=zero

          do iigrid=1,igridstail; igrid=igrids(iigrid);
            ^D&ixOmin^D=ixmlo^D\
            ^D&ixOmax^D=ixmhi^D\
            ^D&ixImin^D=ixglo^D\
            ^D&ixImax^D=ixghi^D\

            ixPmin1=cache(iigrid)%ixPmin1
            ixPmax1=cache(iigrid)%ixPmax1
            ixPmin2=cache(iigrid)%ixPmin2
            ixPmax2=cache(iigrid)%ixPmax2
            has_pixels=cache(iigrid)%has_pixels
            if (.not. has_pixels) cycle

            do ixP2=ixPmin2,ixPmax2
              iFirst=max(ipixStart,(ixP2-1)*nXIF1+ixPmin1)
              iLast=min(ipixEnd,(ixP2-1)*nXIF1+ixPmax1)
              if (iFirst>iLast) cycle
              do ipix=iFirst,iLast
                ixP1=1+mod(ipix-1,nXIF1)
                ray_origin=x_origin+xIF1(ixP1)*vec_xI1+xIF2(ixP2)*vec_xI2
                profile_batch(1)=profile_batch(1)+one
                call ray_box_intersection_cart(ray_origin,vec_LOS,cache(iigrid)%box_min,&
                                               cache(iigrid)%box_max,hit,t_enter,t_exit)
                if (hit) then
                  profile_batch(2)=profile_batch(2)+one
                  nsegBefore=nseg
                  call collect_EUV_cart_dda_segments(ixI^L,ixO^L,cache(iigrid)%source,&
                                                     cache(iigrid)%opacity,cache(iigrid)%sourcev,&
                                                     ipix,ray_origin,cache(iigrid)%xface1,&
                                                     cache(iigrid)%xface2,cache(iigrid)%xface3,&
                                                     t_enter,t_exit,&
                                                     segments,nseg,capacity)
                  profile_batch(3)=profile_batch(3)+dble(nseg-nsegBefore)
                endif
              enddo
            enddo
          enddo

          call MPI_ALLREDUCE(nseg,maxNsegBatch,1,MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
          if (maxNsegBatch>maxSegBatchTarget .and. nPixBatch>1) then
            nPixBatch=max(1,nPixBatch/2)
            ipixEnd=ipixStart+nPixBatch-1
            if (allocated(segments)) deallocate(segments)
            batchReduced=.true.
          else
            batchAccepted=.true.
          endif
        enddo

        profile_local=profile_local+profile_batch
        if (radsyn_verbose .and. mype==0 .and. batchReduced) then
          write(*,'(a,3(i0,1x))') ' cart_dda thick adaptive batch: ',&
               ipixStart,ipixEnd,maxNsegBatch
        endif

        if (.not. allocated(segments)) then
          capacity=1
          allocate(segments(nSegVars,capacity))
        endif
        ownerSegCounts=0
        do is=1,nseg
          owner=segment_pixel_owner(nint(segments(1,is)))
          ownerSegCounts(owner)=ownerSegCounts(owner)+1
        enddo
        sendCounts=nSegVars*ownerSegCounts
        sendDispls(0)=0
        do ipe=1,npe-1
          sendDispls(ipe)=sendDispls(ipe-1)+sendCounts(ipe-1)
        enddo

        allocate(segments_send(nSegVars,max(1,nseg)))
        ownerOffsets=0
        do is=1,nseg
          owner=segment_pixel_owner(nint(segments(1,is)))
          isegDest=sendDispls(owner)/nSegVars+ownerOffsets(owner)+1
          segments_send(:,isegDest)=segments(:,is)
          ownerOffsets(owner)=ownerOffsets(owner)+1
        enddo

        call MPI_ALLTOALL(sendCounts,1,MPI_INTEGER,recvCounts,1,MPI_INTEGER,icomm,ierrmpi)
        recvDispls(0)=0
        do ipe=1,npe-1
          recvDispls(ipe)=recvDispls(ipe-1)+recvCounts(ipe-1)
        enddo
        totalCount=sum(recvCounts)
        totalSeg=totalCount/nSegVars
        profile_local(4)=profile_local(4)+dble(totalCount)
        allocate(segments_recv(nSegVars,max(1,totalSeg)))

        recvFill=0
        maxOwnerSegCountLocal=maxval(ownerSegCounts)
        call MPI_ALLREDUCE(maxOwnerSegCountLocal,maxOwnerSegCount,1,MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
        do segOffset=0,maxOwnerSegCount-1,maxSegCommTarget
          roundSendCounts=0
          roundSendDispls=sendDispls
          do ipe=0,npe-1
            if (ownerSegCounts(ipe)>segOffset) then
              roundSendCounts(ipe)=nSegVars*min(maxSegCommTarget,ownerSegCounts(ipe)-segOffset)
              roundSendDispls(ipe)=sendDispls(ipe)+nSegVars*segOffset
            endif
          enddo

          call MPI_ALLTOALL(roundSendCounts,1,MPI_INTEGER,roundRecvCounts,1,MPI_INTEGER,icomm,ierrmpi)
          roundRecvDispls(0)=0
          do ipe=1,npe-1
            roundRecvDispls(ipe)=roundRecvDispls(ipe-1)+roundRecvCounts(ipe-1)
          enddo
          totalRoundCount=sum(roundRecvCounts)
          totalRoundSeg=totalRoundCount/nSegVars
          allocate(segments_recv_round(nSegVars,max(1,totalRoundSeg)))

          call MPI_ALLTOALLV(segments_send,roundSendCounts,roundSendDispls,MPI_DOUBLE_PRECISION,&
                             segments_recv_round,roundRecvCounts,roundRecvDispls,&
                             MPI_DOUBLE_PRECISION,icomm,ierrmpi)

          if (totalRoundSeg>0) then
            segments_recv(:,recvFill+1:recvFill+totalRoundSeg)=segments_recv_round(:,1:totalRoundSeg)
            recvFill=recvFill+totalRoundSeg
          endif
          deallocate(segments_recv_round)
        enddo

        if (recvFill/=totalSeg) call mpistop("cart_dda thick segmented receive mismatch")

        if (totalSeg>0) then
          allocate(idx(totalSeg))
          bucketCounts(1:nPixBatch)=0
          do is=1,totalSeg
            if (segment_is_valid(segments_recv,is,nSegVars)) then
              ipix=nint(segments_recv(1,is))
              if (ipix>=ipixStart .and. ipix<=ipixEnd .and. segment_pixel_owner(ipix)==mype) then
                iLocal=ipix-ipixStart+1
                bucketCounts(iLocal)=bucketCounts(iLocal)+1
              endif
            endif
          enddo

          bucketOffsets(1)=1
          do iLocal=1,nPixBatch
            bucketOffsets(iLocal+1)=bucketOffsets(iLocal)+bucketCounts(iLocal)
          enddo
          bucketFill(1:nPixBatch)=bucketOffsets(1:nPixBatch)
          do is=1,totalSeg
            if (segment_is_valid(segments_recv,is,nSegVars)) then
              ipix=nint(segments_recv(1,is))
              if (ipix>=ipixStart .and. ipix<=ipixEnd .and. segment_pixel_owner(ipix)==mype) then
                iLocal=ipix-ipixStart+1
                idx(bucketFill(iLocal))=is
                bucketFill(iLocal)=bucketFill(iLocal)+1
              endif
            endif
          enddo

          do ipix=ipixStart,ipixEnd
            if (segment_pixel_owner(ipix)/=mype) cycle
            iLocal=ipix-ipixStart+1
            nidx=bucketCounts(iLocal)
            if (nidx>0) then
              profile_local(5)=profile_local(5)+dble(nidx)*dble(nidx)
              call sort_segment_indices_near_to_far(segments_recv,idx(bucketOffsets(iLocal):bucketOffsets(iLocal+1)-1),nidx)
              ixGlobal=1+mod(ipix-1,nXIF1)
              iyGlobal=1+(ipix-1)/nXIF1
              do iseg=bucketOffsets(iLocal),bucketOffsets(iLocal+1)-1
                is=idx(iseg)
                EUVthin(ixGlobal,iyGlobal)=EUVthin(ixGlobal,iyGlobal)+segments_recv(3,is)
                atten=transfer_attenuation(Tau(ixGlobal,iyGlobal))
                EUV(ixGlobal,iyGlobal)=EUV(ixGlobal,iyGlobal)+atten*segments_recv(3,is)
                Dpl(ixGlobal,iyGlobal)=Dpl(ixGlobal,iyGlobal)+atten*segments_recv(5,is)
                Tau(ixGlobal,iyGlobal)=Tau(ixGlobal,iyGlobal)+max(zero,segments_recv(4,is))
              enddo
            endif
          enddo
          deallocate(idx)
        endif

        deallocate(segments_send,segments_recv)
        if (allocated(segments)) deallocate(segments)
        ipixStart=ipixEnd+1
      enddo

      do iigrid=1,igridstail
        if (allocated(cache(iigrid)%source)) deallocate(cache(iigrid)%source)
        if (allocated(cache(iigrid)%opacity)) deallocate(cache(iigrid)%opacity)
        if (allocated(cache(iigrid)%sourcev)) deallocate(cache(iigrid)%sourcev)
        if (allocated(cache(iigrid)%xface1)) deallocate(cache(iigrid)%xface1)
        if (allocated(cache(iigrid)%xface2)) deallocate(cache(iigrid)%xface2)
        if (allocated(cache(iigrid)%xface3)) deallocate(cache(iigrid)%xface3)
      enddo
      deallocate(cache)
      deallocate(sendCounts,recvCounts,sendDispls,recvDispls,roundSendCounts,roundRecvCounts,&
                 roundSendDispls,roundRecvDispls,ownerSegCounts,ownerOffsets,bucketCounts,&
                 bucketOffsets,bucketFill)
      allocate(image_reduce(nXIF1,nXIF2))
      call MPI_ALLREDUCE(EUV,image_reduce,nXIF1*nXIF2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      EUV=image_reduce
      call MPI_ALLREDUCE(Dpl,image_reduce,nXIF1*nXIF2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      Dpl=image_reduce
      call MPI_ALLREDUCE(Tau,image_reduce,nXIF1*nXIF2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      Tau=image_reduce
      call MPI_ALLREDUCE(EUVthin,image_reduce,nXIF1*nXIF2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      EUVthin=image_reduce
      deallocate(image_reduce)
      call MPI_ALLREDUCE(profile_local,profile_global,5,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      if (radsyn_verbose .and. mype==0) then
        write(*,'(a,5(es12.5,1x))') &
          ' cart_dda thick profile: ',profile_global
      endif
    end subroutine integrate_EUV_cart_dda_thick_datresol

    subroutine integrate_EUV_thick_datresol(nXIF1,nXIF2,fl,EUV,Dpl,Tau,EUVthin)
      use mod_global_parameters

      integer, intent(in) :: nXIF1,nXIF2
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: EUV(nXIF1,nXIF2),Dpl(nXIF1,nXIF2)
      double precision, intent(out) :: Tau(nXIF1,nXIF2),EUVthin(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D
      integer :: iigrid,igrid,levelg,rft,direction_LOS,nLOS,numSeg,nLayerVars,nLayerSeg
      integer :: ixP1,ixP2,ixL,iSub1,iSub2,relL
      integer :: nLosBatch,nBatch,iBatch,ixLstart,ixLend,ixLgridStart,ixLgridEnd
      integer :: ixPmin1,ixPmin2
      double precision :: ds_cm,jds,kds,jvds,atten,layerBytes,targetBytes
      double precision, allocatable :: rho(:^D&)
      double precision, allocatable :: layer_ds(:,:,:,:),layer_all(:,:,:,:)
      type(radsyn_euv_cache), allocatable :: cache(:)

      if (LOS_phi==0 .and. LOS_theta==90) then
        direction_LOS=1
        nLOS=domain_nx1*2**(refine_max_level-1)
      else if (LOS_phi==90 .and. LOS_theta==90) then
        direction_LOS=2
        nLOS=domain_nx2*2**(refine_max_level-1)
      else
        direction_LOS=3
        nLOS=domain_nx3*2**(refine_max_level-1)
      endif

      nLayerVars=3
      if (nXIF1>huge(numSeg)/max(1,nXIF2) .or. nXIF1*nXIF2>huge(numSeg)/nLayerVars) then
        call mpistop("thick EUV layer buffer is too large for one MPI reduction")
      endif
      nLayerSeg=nXIF1*nXIF2*nLayerVars
      targetBytes=256.d0*1024.d0*1024.d0
      layerBytes=dble(nLayerSeg)*8.d0*2.d0
      nLosBatch=max(1,min(16,int(targetBytes/max(one,layerBytes))))
      if (nLayerSeg>huge(numSeg)/nLosBatch) then
        call mpistop("thick EUV batched layer buffer is too large for one MPI reduction")
      endif

      allocate(cache(igridstail))
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixOmax^D=ixmhi^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\

        cache(iigrid)%igrid=igrid
        levelg=ps(igrid)%level
        rft=2**(refine_max_level-levelg)
        cache(iigrid)%level=levelg
        cache(iigrid)%rft=rft

        select case(direction_LOS)
        case(1)
          cache(iigrid)%los_min=(node(pig1_,igrid)-1)*rft*block_nx1+1
          cache(iigrid)%los_max=node(pig1_,igrid)*rft*block_nx1
        case(2)
          cache(iigrid)%los_min=(node(pig2_,igrid)-1)*rft*block_nx2+1
          cache(iigrid)%los_max=node(pig2_,igrid)*rft*block_nx2
        case(3)
          cache(iigrid)%los_min=(node(pig3_,igrid)-1)*rft*block_nx3+1
          cache(iigrid)%los_max=node(pig3_,igrid)*rft*block_nx3
        end select

        allocate(cache(iigrid)%source(ixI^S),cache(iigrid)%opacity(ixI^S),&
                 cache(iigrid)%sourcev(ixI^S),rho(ixI^S))
        cache(iigrid)%source=zero
        cache(iigrid)%opacity=zero
        cache(iigrid)%sourcev=zero
        if (trim(emission_model)=='radio_ff') then
          call get_radio_ff_source_opacity(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,&
                                           cache(iigrid)%source,cache(iigrid)%opacity)
        else
          call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,cache(iigrid)%source)
          cache(iigrid)%source(ixO^S)=cache(iigrid)%source(ixO^S)/instrument_resolution_factor**2
          call get_EUV_HHe_opacity(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,cache(iigrid)%opacity)
          call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
          cache(iigrid)%sourcev(ixO^S)=cache(iigrid)%source(ixO^S)*&
                                       (-ps(igrid)%w(ixO^S,iw_mom(direction_LOS))/rho(ixO^S))
        endif
        deallocate(rho)
      enddo

      ! Stream a small batch of finest LOS layers; this avoids a full Npix*Nlos column buffer.
      allocate(layer_ds(nXIF1,nXIF2,nLayerVars,nLosBatch),layer_all(nXIF1,nXIF2,nLayerVars,nLosBatch))
      EUV=zero
      Dpl=zero
      Tau=zero
      EUVthin=zero

      do ixLstart=1,nLOS,nLosBatch
        ixLend=min(nLOS,ixLstart+nLosBatch-1)
        nBatch=ixLend-ixLstart+1
        layer_ds(:,:,:,1:nBatch)=zero

        do iigrid=1,igridstail
          ixLgridStart=max(ixLstart,cache(iigrid)%los_min)
          ixLgridEnd=min(ixLend,cache(iigrid)%los_max)
          if (ixLgridStart>ixLgridEnd) cycle
          igrid=cache(iigrid)%igrid
          rft=cache(iigrid)%rft
          ^D&ixOmin^D=ixmlo^D\
          ^D&ixOmax^D=ixmhi^D\
          ^D&ixImin^D=ixglo^D\
          ^D&ixImax^D=ixghi^D\

          do ixL=ixLgridStart,ixLgridEnd
            iBatch=ixL-ixLstart+1
            relL=ixL-cache(iigrid)%los_min

            select case(direction_LOS)
            case(1)
              ix1=ixOmin1+relL/rft
              do ix2=ixOmin2,ixOmax2
                ixPmin1=(node(pig2_,igrid)-1)*rft*block_nx2+(ix2-ixOmin2)*rft+1
                do ix3=ixOmin3,ixOmax3
                  ixPmin2=(node(pig3_,igrid)-1)*rft*block_nx3+(ix3-ixOmin3)*rft+1
                  ds_cm=ps(igrid)%dx(ix^D,1)*unit_length/dble(rft)
                  if (SI_unit) ds_cm=ds_cm*1.d2
                  jds=cache(iigrid)%source(ix^D)*ds_cm
                  kds=cache(iigrid)%opacity(ix^D)*ds_cm
                  jvds=cache(iigrid)%sourcev(ix^D)*ds_cm
                  do iSub1=0,rft-1
                    ixP1=ixPmin1+iSub1
                    do iSub2=0,rft-1
                      ixP2=ixPmin2+iSub2
                      layer_ds(ixP1,ixP2,1,iBatch)=layer_ds(ixP1,ixP2,1,iBatch)+jds
                      layer_ds(ixP1,ixP2,2,iBatch)=layer_ds(ixP1,ixP2,2,iBatch)+kds
                      layer_ds(ixP1,ixP2,3,iBatch)=layer_ds(ixP1,ixP2,3,iBatch)+jvds
                    enddo
                  enddo
                enddo
              enddo
            case(2)
              ix2=ixOmin2+relL/rft
              do ix3=ixOmin3,ixOmax3
                ixPmin1=(node(pig3_,igrid)-1)*rft*block_nx3+(ix3-ixOmin3)*rft+1
                do ix1=ixOmin1,ixOmax1
                  ixPmin2=(node(pig1_,igrid)-1)*rft*block_nx1+(ix1-ixOmin1)*rft+1
                  ds_cm=ps(igrid)%dx(ix^D,2)*unit_length/dble(rft)
                  if (SI_unit) ds_cm=ds_cm*1.d2
                  jds=cache(iigrid)%source(ix^D)*ds_cm
                  kds=cache(iigrid)%opacity(ix^D)*ds_cm
                  jvds=cache(iigrid)%sourcev(ix^D)*ds_cm
                  do iSub1=0,rft-1
                    ixP1=ixPmin1+iSub1
                    do iSub2=0,rft-1
                      ixP2=ixPmin2+iSub2
                      layer_ds(ixP1,ixP2,1,iBatch)=layer_ds(ixP1,ixP2,1,iBatch)+jds
                      layer_ds(ixP1,ixP2,2,iBatch)=layer_ds(ixP1,ixP2,2,iBatch)+kds
                      layer_ds(ixP1,ixP2,3,iBatch)=layer_ds(ixP1,ixP2,3,iBatch)+jvds
                    enddo
                  enddo
                enddo
              enddo
            case(3)
              ix3=ixOmin3+relL/rft
              do ix1=ixOmin1,ixOmax1
                ixPmin1=(node(pig1_,igrid)-1)*rft*block_nx1+(ix1-ixOmin1)*rft+1
                do ix2=ixOmin2,ixOmax2
                  ixPmin2=(node(pig2_,igrid)-1)*rft*block_nx2+(ix2-ixOmin2)*rft+1
                  ds_cm=ps(igrid)%dx(ix^D,3)*unit_length/dble(rft)
                  if (SI_unit) ds_cm=ds_cm*1.d2
                  jds=cache(iigrid)%source(ix^D)*ds_cm
                  kds=cache(iigrid)%opacity(ix^D)*ds_cm
                  jvds=cache(iigrid)%sourcev(ix^D)*ds_cm
                  do iSub1=0,rft-1
                    ixP1=ixPmin1+iSub1
                    do iSub2=0,rft-1
                      ixP2=ixPmin2+iSub2
                      layer_ds(ixP1,ixP2,1,iBatch)=layer_ds(ixP1,ixP2,1,iBatch)+jds
                      layer_ds(ixP1,ixP2,2,iBatch)=layer_ds(ixP1,ixP2,2,iBatch)+kds
                      layer_ds(ixP1,ixP2,3,iBatch)=layer_ds(ixP1,ixP2,3,iBatch)+jvds
                    enddo
                  enddo
                enddo
              enddo
            end select
          enddo
        enddo

        numSeg=nLayerSeg*nBatch
        call MPI_ALLREDUCE(layer_ds,layer_all,numSeg,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)

        do iBatch=1,nBatch
          !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(atten) SCHEDULE(static)
          do ixP1=1,nXIF1
            do ixP2=1,nXIF2
              EUVthin(ixP1,ixP2)=EUVthin(ixP1,ixP2)+layer_all(ixP1,ixP2,1,iBatch)
              atten=transfer_attenuation(Tau(ixP1,ixP2))
              EUV(ixP1,ixP2)=EUV(ixP1,ixP2)+atten*layer_all(ixP1,ixP2,1,iBatch)
              Dpl(ixP1,ixP2)=Dpl(ixP1,ixP2)+atten*layer_all(ixP1,ixP2,3,iBatch)
              Tau(ixP1,ixP2)=Tau(ixP1,ixP2)+layer_all(ixP1,ixP2,2,iBatch)
            enddo
          enddo
          !$OMP END PARALLEL DO
        enddo
      enddo

      do iigrid=1,igridstail
        if (allocated(cache(iigrid)%source)) deallocate(cache(iigrid)%source)
        if (allocated(cache(iigrid)%opacity)) deallocate(cache(iigrid)%opacity)
        if (allocated(cache(iigrid)%sourcev)) deallocate(cache(iigrid)%sourcev)
      enddo
      deallocate(cache,layer_ds,layer_all)

    end subroutine integrate_EUV_thick_datresol

  }

  {^IFTHREED

    subroutine build_sph_intersection_faces(ixI^L,ixO^L,x,dx,rface,thetaface,phiface)
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim),dx(ixI^S,1:ndim)
      double precision, allocatable, intent(out) :: rface(:),thetaface(:),phiface(:)
      integer :: ix1,ix2,ix3

      allocate(rface(ixOmin1:ixOmax1+1),thetaface(ixOmin2:ixOmax2+1),phiface(ixOmin3:ixOmax3+1))
      do ix1=ixOmin1,ixOmax1
        rface(ix1)=x(ix1,ixOmin2,ixOmin3,1)-half*dx(ix1,ixOmin2,ixOmin3,1)
      enddo
      rface(ixOmax1+1)=x(ixOmax1,ixOmin2,ixOmin3,1)+half*dx(ixOmax1,ixOmin2,ixOmin3,1)
      do ix2=ixOmin2,ixOmax2
        thetaface(ix2)=x(ixOmin1,ix2,ixOmin3,2)-half*dx(ixOmin1,ix2,ixOmin3,2)
      enddo
      thetaface(ixOmax2+1)=x(ixOmin1,ixOmax2,ixOmin3,2)+half*dx(ixOmin1,ixOmax2,ixOmin3,2)
      do ix3=ixOmin3,ixOmax3
        phiface(ix3)=x(ixOmin1,ixOmin2,ix3,3)-half*dx(ixOmin1,ixOmin2,ix3,3)
      enddo
      phiface(ixOmax3+1)=x(ixOmin1,ixOmin2,ixOmax3,3)+half*dx(ixOmin1,ixOmin2,ixOmax3,3)
    end subroutine build_sph_intersection_faces

    subroutine sph_add_t(tvals,nt,capacity,t)
      double precision, allocatable, intent(inout) :: tvals(:)
      integer, intent(inout) :: nt,capacity
      double precision, intent(in) :: t

      double precision, allocatable :: tmp(:)

      if (t /= t .or. abs(t)>1.d90) return
      if (.not. allocated(tvals)) then
        capacity=64
        allocate(tvals(capacity))
      else if (nt>=capacity) then
        allocate(tmp(capacity))
        tmp=tvals
        deallocate(tvals)
        allocate(tvals(2*capacity))
        tvals(1:capacity)=tmp
        deallocate(tmp)
        capacity=2*capacity
      endif
      nt=nt+1
      tvals(nt)=t
    end subroutine sph_add_t

    subroutine sph_add_t_fixed(tvals,nt,t)
      double precision, intent(inout) :: tvals(:)
      integer, intent(inout) :: nt
      double precision, intent(in) :: t

      if (t /= t .or. abs(t)>1.d90) return
      if (nt>=size(tvals)) return
      nt=nt+1
      tvals(nt)=t
    end subroutine sph_add_t_fixed

    subroutine sph_sort_unique_t(tvals,nt)
      double precision, intent(inout) :: tvals(:)
      integer, intent(inout) :: nt

      integer :: i,j,nout
      double precision :: key,epsT

      if (nt<=1) return
      do i=2,nt
        key=tvals(i)
        j=i-1
        do while (j>=1 .and. tvals(j)>key)
          tvals(j+1)=tvals(j)
          j=j-1
        enddo
        tvals(j+1)=key
      enddo
      epsT=max(1.d-12,1.d-10*max(one,abs(tvals(nt)-tvals(1))))
      nout=1
      do i=2,nt
        if (abs(tvals(i)-tvals(nout))>epsT) then
          nout=nout+1
          tvals(nout)=tvals(i)
        endif
      enddo
      nt=nout
    end subroutine sph_sort_unique_t

    subroutine sph_add_sphere_intersections(ray_origin,ray_dir,rface,tvals,nt,capacity)
      double precision, intent(in) :: ray_origin(1:3),ray_dir(1:3),rface
      double precision, allocatable, intent(inout) :: tvals(:)
      integer, intent(inout) :: nt,capacity

      double precision :: aa,bb,cc,disc,root

      aa=sum(ray_dir**2)
      bb=2.d0*sum(ray_origin*ray_dir)
      cc=sum(ray_origin**2)-rface**2
      disc=bb**2-4.d0*aa*cc
      if (disc<zero) return
      root=sqrt(max(zero,disc))
      call sph_add_t(tvals,nt,capacity,(-bb-root)/(2.d0*aa))
      call sph_add_t(tvals,nt,capacity,(-bb+root)/(2.d0*aa))
    end subroutine sph_add_sphere_intersections

    subroutine sph_add_theta_intersections(ray_origin,ray_dir,thetaface,tvals,nt,capacity)
      double precision, intent(in) :: ray_origin(1:3),ray_dir(1:3),thetaface
      double precision, allocatable, intent(inout) :: tvals(:)
      integer, intent(inout) :: nt,capacity

      double precision :: cth,aa,bb,cc,disc,root

      cth=cos(thetaface)
      ! At the equator the cone equation degenerates to z**2=0.  Solve the
      ! corresponding plane intersection directly to avoid losing its double
      ! root through roundoff in the quadratic discriminant.
      if (abs(cth)<1.d-12) then
        if (abs(ray_dir(3))>1.d-14) &
          call sph_add_t(tvals,nt,capacity,-ray_origin(3)/ray_dir(3))
        return
      endif
      aa=ray_dir(3)**2-cth**2*sum(ray_dir**2)
      bb=2.d0*(ray_origin(3)*ray_dir(3)-cth**2*sum(ray_origin*ray_dir))
      cc=ray_origin(3)**2-cth**2*sum(ray_origin**2)
      if (abs(aa)<1.d-14) then
        if (abs(bb)>1.d-14) call sph_add_t(tvals,nt,capacity,-cc/bb)
        return
      endif
      disc=bb**2-4.d0*aa*cc
      if (disc<zero) return
      root=sqrt(max(zero,disc))
      call sph_add_t(tvals,nt,capacity,(-bb-root)/(2.d0*aa))
      call sph_add_t(tvals,nt,capacity,(-bb+root)/(2.d0*aa))
    end subroutine sph_add_theta_intersections

    subroutine sph_add_phi_intersection(ray_origin,ray_dir,phiface,tvals,nt,capacity)
      double precision, intent(in) :: ray_origin(1:3),ray_dir(1:3),phiface
      double precision, allocatable, intent(inout) :: tvals(:)
      integer, intent(inout) :: nt,capacity

      double precision :: normal(1:3),denom,numer

      normal(1)=-sin(phiface)
      normal(2)=cos(phiface)
      normal(3)=zero
      denom=sum(normal*ray_dir)
      if (abs(denom)<1.d-14) return
      numer=sum(normal*ray_origin)
      call sph_add_t(tvals,nt,capacity,-numer/denom)
    end subroutine sph_add_phi_intersection

    subroutine sph_cart_to_coord(pos,sph)
      double precision, intent(in) :: pos(1:3)
      double precision, intent(out) :: sph(1:3)

      sph(1)=sqrt(sum(pos**2))
      if (sph(1)>zero) then
        sph(2)=acos(max(-one,min(one,pos(3)/sph(1))))
      else
        sph(2)=zero
      endif
      sph(3)=atan2(pos(2),pos(1))
    end subroutine sph_cart_to_coord

    integer function sph_locate_index(value,faces,imin,imax) result(idx)
      integer, intent(in) :: imin,imax
      double precision, intent(in) :: value,faces(imin:imax+1)

      integer :: ilo,ihi,imid

      idx=0
      if (value<faces(imin)-1.d-12 .or. value>faces(imax+1)+1.d-12) return
      if (value<=faces(imin)) then
        idx=imin
        return
      endif
      if (value>=faces(imax+1)) then
        idx=imax
        return
      endif

      ilo=imin
      ihi=imax+1
      do while (ihi-ilo>1)
        imid=(ilo+ihi)/2
        if (value>=faces(imid)) then
          ilo=imid
        else
          ihi=imid
        endif
      enddo
      idx=min(imax,max(imin,ilo))
    end function sph_locate_index

    integer function sph_locate_index_desc(value,faces,imin,imax) result(idx)
      integer, intent(in) :: imin,imax
      double precision, intent(in) :: value,faces(imin:imax+1)

      integer :: ilo,ihi,imid

      idx=0
      if (value>faces(imin)+1.d-12 .or. value<faces(imax+1)-1.d-12) return
      if (value>=faces(imin)) then
        idx=imin
        return
      endif
      if (value<=faces(imax+1)) then
        idx=imax
        return
      endif

      ilo=imin
      ihi=imax+1
      do while (ihi-ilo>1)
        imid=(ilo+ihi)/2
        if (value<=faces(imid)) then
          ilo=imid
        else
          ihi=imid
        endif
      enddo
      idx=min(imax,max(imin,ilo))
    end function sph_locate_index_desc

    subroutine sph_locate_cell(pos,rface,thetaface,phiface,ixO^L,ix1,ix2,ix3,inside)
      double precision, intent(in) :: pos(1:3)
      double precision, intent(in) :: rface(ixOmin1:ixOmax1+1),thetaface(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      integer, intent(in) :: ixO^L
      integer, intent(out) :: ix1,ix2,ix3
      logical, intent(out) :: inside

      double precision :: sph(1:3),phi

      call sph_cart_to_coord(pos,sph)
      phi=sph(3)
      if (phi<phiface(ixOmin3)-1.d-12) phi=phi+2.d0*dpi
      if (phi>phiface(ixOmax3+1)+1.d-12) phi=phi-2.d0*dpi
      ix1=sph_locate_index(sph(1),rface,ixOmin1,ixOmax1)
      ix2=sph_locate_index(sph(2),thetaface,ixOmin2,ixOmax2)
      ix3=sph_locate_index(phi,phiface,ixOmin3,ixOmax3)
      inside=ix1>0 .and. ix2>0 .and. ix3>0
    end subroutine sph_locate_cell

    logical function sph_segment_visible(pos,ximg1,ximg2) result(visible)
      use mod_constants
      double precision, intent(in) :: pos(1:3),ximg1,ximg2

      double precision :: dotp,rc,rthick,rloc

      rthick=R_opt_thick*const_Rsun/unit_length
      rc=sqrt(ximg1**2+ximg2**2)
      rloc=sqrt(sum(pos**2))
      call dot_product_loc(vec_LOS,pos,dotp)
      visible=.true.
      if (dotp>=zero) then
        if (rc<=rthick) visible=.false.
      else
        if (rloc<=rthick) visible=.false.
      endif
    end function sph_segment_visible

    subroutine sph_block_pixel_range(rface,thetaface,phiface,ixO^L,nXI1,nXI2,xI1,xI2,dxI,&
                                     ixPmin1,ixPmax1,ixPmin2,ixPmax2,has_pixels)
      double precision, intent(in) :: rface(ixOmin1:ixOmax1+1),thetaface(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      integer, intent(in) :: ixO^L,nXI1,nXI2
      double precision, intent(in) :: xI1(nXI1),xI2(nXI2),dxI
      integer, intent(out) :: ixPmin1,ixPmax1,ixPmin2,ixPmax2
      logical, intent(out) :: has_pixels

      integer, parameter :: nsample=5
      integer :: ir,it,ip
      double precision :: sph(1:3),xcent(1:2)
      double precision :: xmin1,xmax1,xmin2,xmax2
      double precision :: wr,wt,wp,pad

      has_pixels=.false.
      xmin1=huge(one)
      xmax1=-huge(one)
      xmin2=huge(one)
      xmax2=-huge(one)
      do ir=0,nsample-1
        wr=dble(ir)/dble(nsample-1)
        sph(1)=(one-wr)*rface(ixOmin1)+wr*rface(ixOmax1+1)
        do it=0,nsample-1
          wt=dble(it)/dble(nsample-1)
          sph(2)=(one-wt)*thetaface(ixOmin2)+wt*thetaface(ixOmax2+1)
          do ip=0,nsample-1
            wp=dble(ip)/dble(nsample-1)
            if (ir/=0 .and. ir/=nsample-1 .and. it/=0 .and. it/=nsample-1 .and. &
                ip/=0 .and. ip/=nsample-1) cycle
            sph(3)=(one-wp)*phiface(ixOmin3)+wp*phiface(ixOmax3+1)
            call get_cor_image_spherical(sph,xcent)
            xmin1=min(xmin1,xcent(1))
            xmax1=max(xmax1,xcent(1))
            xmin2=min(xmin2,xcent(2))
            xmax2=max(xmax2,xcent(2))
          enddo
        enddo
      enddo
      pad=2.d0*dxI
      xmin1=xmin1-pad
      xmax1=xmax1+pad
      xmin2=xmin2-pad
      xmax2=xmax2+pad
      ixPmin1=max(1,floor((xmin1-(xI1(1)-half*dxI))/dxI)+1)
      ixPmax1=min(nXI1,ceiling((xmax1-(xI1(1)-half*dxI))/dxI))
      ixPmin2=max(1,floor((xmin2-(xI2(1)-half*dxI))/dxI)+1)
      ixPmax2=min(nXI2,ceiling((xmax2-(xI2(1)-half*dxI))/dxI))
      has_pixels=ixPmin1<=ixPmax1 .and. ixPmin2<=ixPmax2
    end subroutine sph_block_pixel_range

    subroutine acc_EUV_sph_intersection(ixI^L,ixO^L,source,ray_origin,ximg1,ximg2,&
                                        rface,thetaface,phiface,EUVp)
      integer, intent(in) :: ixI^L,ixO^L
      double precision, intent(in) :: source(ixI^S),ray_origin(1:3),ximg1,ximg2
      double precision, intent(in) :: rface(ixOmin1:ixOmax1+1),thetaface(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      double precision, intent(inout) :: EUVp

      integer :: nt,capacity,i,ix^D
      double precision, allocatable :: tvals(:)
      double precision :: posMid(1:3),ds_cm,tMid,t0,t1
      logical :: inside

      nt=0
      capacity=0
      do ix1=ixOmin1,ixOmax1+1
        call sph_add_sphere_intersections(ray_origin,vec_LOS,rface(ix1),tvals,nt,capacity)
      enddo
      do ix2=ixOmin2,ixOmax2+1
        call sph_add_theta_intersections(ray_origin,vec_LOS,thetaface(ix2),tvals,nt,capacity)
      enddo
      do ix3=ixOmin3,ixOmax3+1
        call sph_add_phi_intersection(ray_origin,vec_LOS,phiface(ix3),tvals,nt,capacity)
      enddo
      if (nt<2) then
        if (allocated(tvals)) deallocate(tvals)
        return
      endif
      call sph_sort_unique_t(tvals,nt)

      do i=1,nt-1
        t0=tvals(i)
        t1=tvals(i+1)
        if (t1<=t0) cycle
        tMid=half*(t0+t1)
        posMid=ray_origin+tMid*vec_LOS
        if (.not. sph_segment_visible(posMid,ximg1,ximg2)) cycle
        call sph_locate_cell(posMid,rface,thetaface,phiface,ixO^L,ix1,ix2,ix3,inside)
        if (.not. inside) cycle
        ds_cm=(t1-t0)*unit_length
        if (SI_unit) ds_cm=ds_cm*1.d2
        EUVp=EUVp+source(ix^D)*ds_cm
      enddo
      deallocate(tvals)
    end subroutine acc_EUV_sph_intersection

    subroutine collect_EUV_sph_intersection_segments(ixI^L,ixO^L,source,opacity,&
                                                     pixel_id,ray_origin,ximg1,ximg2,&
                                                     rface,thetaface,phiface,rface2,&
                                                     theta_cos,phi_sin,phi_cos,&
                                                     segments,nseg,capacity)
      use mod_constants, only: const_Rsun

      integer, intent(in) :: ixI^L,ixO^L,pixel_id
      double precision, intent(in) :: source(ixI^S),opacity(ixI^S)
      double precision, intent(in) :: ray_origin(1:3),ximg1,ximg2
      double precision, intent(in) :: rface(ixOmin1:ixOmax1+1),thetaface(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      double precision, intent(in) :: rface2(ixOmin1:ixOmax1+1),theta_cos(ixOmin2:ixOmax2+1),&
                                      phi_sin(ixOmin3:ixOmax3+1),phi_cos(ixOmin3:ixOmax3+1)
      double precision, allocatable, intent(inout) :: segments(:,:)
      integer, intent(inout) :: nseg,capacity

      integer :: nt,i,ix^D
      double precision :: tvals(2*(ixOmax1-ixOmin1+2)+2*(ixOmax2-ixOmin2+2)+&
                                (ixOmax3-ixOmin3+2))
      double precision :: posMid(1:3),ds_cm,tMid,t0,t1,jds,kds
      double precision :: dir2,origin2,odotd,aa,bb,cc,disc,root,cth2
      double precision :: denom,numer,r2,mu,phi,dotp,rthick2,rc2
      logical :: inside

      nt=0
      dir2=sum(vec_LOS**2)
      origin2=sum(ray_origin**2)
      odotd=sum(ray_origin*vec_LOS)
      rthick2=(R_opt_thick*const_Rsun/unit_length)**2
      rc2=ximg1**2+ximg2**2

      do ix1=ixOmin1,ixOmax1+1
        aa=dir2
        bb=2.d0*odotd
        cc=origin2-rface2(ix1)
        disc=bb**2-4.d0*aa*cc
        if (disc>=zero) then
          root=sqrt(max(zero,disc))
          call sph_add_t_fixed(tvals,nt,(-bb-root)/(2.d0*aa))
          call sph_add_t_fixed(tvals,nt,(-bb+root)/(2.d0*aa))
        endif
      enddo
      do ix2=ixOmin2,ixOmax2+1
        if (abs(theta_cos(ix2))<1.d-12) then
          if (abs(vec_LOS(3))>1.d-14) &
            call sph_add_t_fixed(tvals,nt,-ray_origin(3)/vec_LOS(3))
          cycle
        endif
        cth2=theta_cos(ix2)**2
        aa=vec_LOS(3)**2-cth2*dir2
        bb=2.d0*(ray_origin(3)*vec_LOS(3)-cth2*odotd)
        cc=ray_origin(3)**2-cth2*origin2
        if (abs(aa)<1.d-14) then
          if (abs(bb)>1.d-14) call sph_add_t_fixed(tvals,nt,-cc/bb)
        else
          disc=bb**2-4.d0*aa*cc
          if (disc>=zero) then
            root=sqrt(max(zero,disc))
            call sph_add_t_fixed(tvals,nt,(-bb-root)/(2.d0*aa))
            call sph_add_t_fixed(tvals,nt,(-bb+root)/(2.d0*aa))
          endif
        endif
      enddo
      do ix3=ixOmin3,ixOmax3+1
        denom=-phi_sin(ix3)*vec_LOS(1)+phi_cos(ix3)*vec_LOS(2)
        if (abs(denom)>=1.d-14) then
          numer=-phi_sin(ix3)*ray_origin(1)+phi_cos(ix3)*ray_origin(2)
          call sph_add_t_fixed(tvals,nt,-numer/denom)
        endif
      enddo
      if (nt<2) return
      call sph_sort_unique_t(tvals,nt)

      do i=1,nt-1
        t0=tvals(i)
        t1=tvals(i+1)
        if (t1<=t0) cycle
        tMid=half*(t0+t1)
        posMid=ray_origin+tMid*vec_LOS
        r2=sum(posMid**2)
        dotp=sum(vec_LOS*posMid)
        if (dotp>=zero) then
          if (rc2<=rthick2) cycle
        else
          if (r2<=rthick2) cycle
        endif

        ix1=sph_locate_index(r2,rface2,ixOmin1,ixOmax1)
        if (r2>zero) then
          mu=posMid(3)/sqrt(r2)
        else
          mu=one
        endif
        ix2=sph_locate_index_desc(mu,theta_cos,ixOmin2,ixOmax2)
        phi=atan2(posMid(2),posMid(1))
        if (phi<phiface(ixOmin3)-1.d-12) phi=phi+2.d0*dpi
        if (phi>phiface(ixOmax3+1)+1.d-12) phi=phi-2.d0*dpi
        ix3=sph_locate_index(phi,phiface,ixOmin3,ixOmax3)
        inside=ix1>0 .and. ix2>0 .and. ix3>0
        if (.not. inside) cycle
        ds_cm=(t1-t0)*unit_length
        if (SI_unit) ds_cm=ds_cm*1.d2
        jds=max(zero,source(ix^D))*ds_cm
        kds=max(zero,opacity(ix^D))*ds_cm
        call append_cart_dda_segment(segments,nseg,capacity,pixel_id,tMid,jds,kds,zero)
      enddo
    end subroutine collect_EUV_sph_intersection_segments

    subroutine sph_locate_cell_fast(pos,rface2,theta_cos,phiface,ixO^L,ix1,ix2,ix3,inside)
      double precision, intent(in) :: pos(1:3)
      double precision, intent(in) :: rface2(ixOmin1:ixOmax1+1),theta_cos(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      integer, intent(in) :: ixO^L
      integer, intent(out) :: ix1,ix2,ix3
      logical, intent(out) :: inside

      double precision :: r2,mu,phi

      r2=sum(pos**2)
      ix1=sph_locate_index(r2,rface2,ixOmin1,ixOmax1)
      if (r2>zero) then
        mu=pos(3)/sqrt(r2)
      else
        mu=one
      endif
      ix2=sph_locate_index_desc(mu,theta_cos,ixOmin2,ixOmax2)
      phi=atan2(pos(2),pos(1))
      if (phi<phiface(ixOmin3)-1.d-12) phi=phi+2.d0*dpi
      if (phi>phiface(ixOmax3+1)+1.d-12) phi=phi-2.d0*dpi
      ix3=sph_locate_index(phi,phiface,ixOmin3,ixOmax3)
      inside=ix1>0 .and. ix2>0 .and. ix3>0
    end subroutine sph_locate_cell_fast

    subroutine sph_try_exit_candidate(t,tNow,tExit,epsRay,tNext,found)
      double precision, intent(in) :: t,tNow,tExit,epsRay
      double precision, intent(inout) :: tNext
      logical, intent(inout) :: found

      if (t>tNow+epsRay .and. t<=tExit+epsRay .and. t<tNext) then
        tNext=t
        found=.true.
      endif
    end subroutine sph_try_exit_candidate

    subroutine sph_try_theta_exit_candidate(t,theta_face_cos,ray_origin,tNow,tExit,epsRay,tNext,found)
      double precision, intent(in) :: t,theta_face_cos,ray_origin(1:3),tNow,tExit,epsRay
      double precision, intent(inout) :: tNext
      logical, intent(inout) :: found

      double precision :: pos(1:3),r2

      if (t<=tNow+epsRay .or. t>tExit+epsRay .or. t>=tNext) return
      pos=ray_origin+t*vec_LOS
      r2=sum(pos**2)
      if (r2<=zero) return
      if (theta_face_cos>1.d-12 .and. pos(3)<-1.d-10) return
      if (theta_face_cos<-1.d-12 .and. pos(3)>1.d-10) return
      if (abs(pos(3)**2-theta_face_cos**2*r2)>1.d-6*max(one,r2)) return
      tNext=t
      found=.true.
    end subroutine sph_try_theta_exit_candidate

    subroutine sph_try_phi_exit_candidate(t,phi_face_sin,phi_face_cos,ray_origin,tNow,tExit,epsRay,tNext,found)
      double precision, intent(in) :: t,phi_face_sin,phi_face_cos,ray_origin(1:3),tNow,tExit,epsRay
      double precision, intent(inout) :: tNext
      logical, intent(inout) :: found

      double precision :: pos(1:3),radialDot

      if (t<=tNow+epsRay .or. t>tExit+epsRay .or. t>=tNext) return
      pos=ray_origin+t*vec_LOS
      radialDot=phi_face_cos*pos(1)+phi_face_sin*pos(2)
      if (radialDot<-1.d-10) return
      tNext=t
      found=.true.
    end subroutine sph_try_phi_exit_candidate

    subroutine sph_next_cell_exit(ray_origin,rface2,theta_cos,phiface,phi_sin,phi_cos,ixO^L,&
                                  ix1,ix2,ix3,tNow,tExit,epsRay,tNext,found)
      double precision, intent(in) :: ray_origin(1:3)
      double precision, intent(in) :: rface2(ixOmin1:ixOmax1+1),theta_cos(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      double precision, intent(in) :: phi_sin(ixOmin3:ixOmax3+1),phi_cos(ixOmin3:ixOmax3+1)
      integer, intent(in) :: ixO^L,ix1,ix2,ix3
      double precision, intent(in) :: tNow,tExit,epsRay
      double precision, intent(out) :: tNext
      logical, intent(out) :: found

      integer :: iface
      double precision :: dir2,origin2,odotd,aa,bb,cc,disc,root,cth2,denom,numer

      found=.false.
      tNext=huge(one)
      dir2=sum(vec_LOS**2)
      origin2=sum(ray_origin**2)
      odotd=sum(ray_origin*vec_LOS)

      do iface=ix1,ix1+1
        aa=dir2
        bb=2.d0*odotd
        cc=origin2-rface2(iface)
        disc=bb**2-4.d0*aa*cc
        if (disc>=zero) then
          root=sqrt(max(zero,disc))
          call sph_try_exit_candidate((-bb-root)/(2.d0*aa),tNow,tExit,epsRay,tNext,found)
          call sph_try_exit_candidate((-bb+root)/(2.d0*aa),tNow,tExit,epsRay,tNext,found)
        endif
      enddo

      do iface=ix2,ix2+1
        if (abs(theta_cos(iface))<1.d-12) then
          if (abs(vec_LOS(3))>1.d-14) call sph_try_theta_exit_candidate(&
            -ray_origin(3)/vec_LOS(3),theta_cos(iface),ray_origin,tNow,tExit,epsRay,&
            tNext,found)
          cycle
        endif
        cth2=theta_cos(iface)**2
        aa=vec_LOS(3)**2-cth2*dir2
        bb=2.d0*(ray_origin(3)*vec_LOS(3)-cth2*odotd)
        cc=ray_origin(3)**2-cth2*origin2
        if (abs(aa)<1.d-14) then
          if (abs(bb)>1.d-14) call sph_try_theta_exit_candidate(-cc/bb,theta_cos(iface),&
               ray_origin,tNow,tExit,epsRay,tNext,found)
        else
          disc=bb**2-4.d0*aa*cc
          if (disc>=zero) then
            root=sqrt(max(zero,disc))
            call sph_try_theta_exit_candidate((-bb-root)/(2.d0*aa),theta_cos(iface),&
                 ray_origin,tNow,tExit,epsRay,tNext,found)
            call sph_try_theta_exit_candidate((-bb+root)/(2.d0*aa),theta_cos(iface),&
                 ray_origin,tNow,tExit,epsRay,tNext,found)
          endif
        endif
      enddo

      do iface=ix3,ix3+1
        denom=-phi_sin(iface)*vec_LOS(1)+phi_cos(iface)*vec_LOS(2)
        if (abs(denom)>=1.d-14) then
          numer=-phi_sin(iface)*ray_origin(1)+phi_cos(iface)*ray_origin(2)
          call sph_try_phi_exit_candidate(-numer/denom,phi_sin(iface),phi_cos(iface),ray_origin,&
               tNow,tExit,epsRay,tNext,found)
        endif
      enddo
    end subroutine sph_next_cell_exit

    subroutine collect_EUV_sph_dda_interval(ixI^L,ixO^L,source,opacity,pixel_id,&
                                            ray_origin,ximg1,ximg2,rface2,theta_cos,phiface,&
                                            phi_sin,phi_cos,&
                                            t_enter,t_exit,segments,nseg,capacity,ok)
      use mod_constants, only: const_Rsun

      integer, intent(in) :: ixI^L,ixO^L,pixel_id
      double precision, intent(in) :: source(ixI^S),opacity(ixI^S)
      double precision, intent(in) :: ray_origin(1:3),ximg1,ximg2
      double precision, intent(in) :: rface2(ixOmin1:ixOmax1+1),theta_cos(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      double precision, intent(in) :: phi_sin(ixOmin3:ixOmax3+1),phi_cos(ixOmin3:ixOmax3+1)
      double precision, intent(in) :: t_enter,t_exit
      double precision, allocatable, intent(inout) :: segments(:,:)
      integer, intent(inout) :: nseg,capacity
      logical, intent(out) :: ok

      integer :: ix^D,nstep,maxSteps
      double precision :: tNow,tNext,tEnd,tMid,epsRay,ds_cm,jds,kds
      double precision :: pos(1:3),r2,dotp,rthick2,rc2
      logical :: inside,found

      ok=.true.
      if (t_exit<=t_enter) return
      epsRay=max(1.d-12,1.d-10*max(one,abs(t_exit-t_enter)))
      pos=ray_origin+(t_enter+epsRay)*vec_LOS
      call sph_locate_cell_fast(pos,rface2,theta_cos,phiface,ixO^L,ix1,ix2,ix3,inside)
      if (.not. inside) then
        ok=.false.
        return
      endif

      rthick2=(R_opt_thick*const_Rsun/unit_length)**2
      rc2=ximg1**2+ximg2**2
      tNow=t_enter
      nstep=0
      maxSteps=8*((ixOmax1-ixOmin1+1)+(ixOmax2-ixOmin2+1)+(ixOmax3-ixOmin3+1)+3)

      do while (tNow<t_exit-epsRay)
        call sph_next_cell_exit(ray_origin,rface2,theta_cos,phiface,phi_sin,phi_cos,ixO^L,&
                                ix1,ix2,ix3,tNow,t_exit,epsRay,tNext,found)
        if (.not. found) then
          ok=.false.
          return
        endif
        tEnd=min(tNext,t_exit)
        if (tEnd>tNow) then
          tMid=half*(tNow+tEnd)
          pos=ray_origin+tMid*vec_LOS
          r2=sum(pos**2)
          dotp=sum(vec_LOS*pos)
          if (.not. ((dotp>=zero .and. rc2<=rthick2) .or. &
                     (dotp<zero .and. r2<=rthick2))) then
            ds_cm=(tEnd-tNow)*unit_length
            if (SI_unit) ds_cm=ds_cm*1.d2
            jds=max(zero,source(ix^D))*ds_cm
            kds=max(zero,opacity(ix^D))*ds_cm
            call append_cart_dda_segment(segments,nseg,capacity,pixel_id,tMid,jds,kds,zero)
          endif
        endif
        tNow=tEnd
        if (tNow>=t_exit-epsRay) exit
        pos=ray_origin+(tNow+epsRay)*vec_LOS
        call sph_locate_cell_fast(pos,rface2,theta_cos,phiface,ixO^L,ix1,ix2,ix3,inside)
        if (.not. inside) then
          ok=.false.
          return
        endif
        nstep=nstep+1
        if (nstep>maxSteps) then
          ok=.false.
          return
        endif
      enddo
    end subroutine collect_EUV_sph_dda_interval

    subroutine collect_EUV_sph_dda_segments(ixI^L,ixO^L,source,opacity,&
                                            pixel_id,ray_origin,ximg1,ximg2,&
                                            rface,thetaface,phiface,rface2,&
                                            theta_cos,phi_sin,phi_cos,&
                                            segments,nseg,capacity,fallback)
      integer, intent(in) :: ixI^L,ixO^L,pixel_id
      double precision, intent(in) :: source(ixI^S),opacity(ixI^S)
      double precision, intent(in) :: ray_origin(1:3),ximg1,ximg2
      double precision, intent(in) :: rface(ixOmin1:ixOmax1+1),thetaface(ixOmin2:ixOmax2+1),&
                                      phiface(ixOmin3:ixOmax3+1)
      double precision, intent(in) :: rface2(ixOmin1:ixOmax1+1),theta_cos(ixOmin2:ixOmax2+1),&
                                      phi_sin(ixOmin3:ixOmax3+1),phi_cos(ixOmin3:ixOmax3+1)
      double precision, allocatable, intent(inout) :: segments(:,:)
      integer, intent(inout) :: nseg,capacity
      logical, intent(out) :: fallback

      integer :: nt,i,nsegStart,ix^D
      double precision :: tvals(12),posMid(1:3),t0,t1,tMid
      double precision :: dir2,origin2,odotd,aa,bb,cc,disc,root,cth2,denom,numer
      logical :: inside,ok

      fallback=.false.
      nsegStart=nseg
      nt=0
      dir2=sum(vec_LOS**2)
      origin2=sum(ray_origin**2)
      odotd=sum(ray_origin*vec_LOS)

      do ix1=ixOmin1,ixOmax1+1,ixOmax1-ixOmin1+1
        aa=dir2
        bb=2.d0*odotd
        cc=origin2-rface2(ix1)
        disc=bb**2-4.d0*aa*cc
        if (disc>=zero) then
          root=sqrt(max(zero,disc))
          call sph_add_t_fixed(tvals,nt,(-bb-root)/(2.d0*aa))
          call sph_add_t_fixed(tvals,nt,(-bb+root)/(2.d0*aa))
        endif
      enddo
      do ix2=ixOmin2,ixOmax2+1,ixOmax2-ixOmin2+1
        if (abs(theta_cos(ix2))<1.d-12) then
          if (abs(vec_LOS(3))>1.d-14) &
            call sph_add_t_fixed(tvals,nt,-ray_origin(3)/vec_LOS(3))
          cycle
        endif
        cth2=theta_cos(ix2)**2
        aa=vec_LOS(3)**2-cth2*dir2
        bb=2.d0*(ray_origin(3)*vec_LOS(3)-cth2*odotd)
        cc=ray_origin(3)**2-cth2*origin2
        if (abs(aa)<1.d-14) then
          if (abs(bb)>1.d-14) call sph_add_t_fixed(tvals,nt,-cc/bb)
        else
          disc=bb**2-4.d0*aa*cc
          if (disc>=zero) then
            root=sqrt(max(zero,disc))
            call sph_add_t_fixed(tvals,nt,(-bb-root)/(2.d0*aa))
            call sph_add_t_fixed(tvals,nt,(-bb+root)/(2.d0*aa))
          endif
        endif
      enddo
      do ix3=ixOmin3,ixOmax3+1,ixOmax3-ixOmin3+1
        denom=-phi_sin(ix3)*vec_LOS(1)+phi_cos(ix3)*vec_LOS(2)
        if (abs(denom)>=1.d-14) then
          numer=-phi_sin(ix3)*ray_origin(1)+phi_cos(ix3)*ray_origin(2)
          call sph_add_t_fixed(tvals,nt,-numer/denom)
        endif
      enddo

      if (nt<2) return
      call sph_sort_unique_t(tvals,nt)

      do i=1,nt-1
        t0=tvals(i)
        t1=tvals(i+1)
        if (t1<=t0) cycle
        tMid=half*(t0+t1)
        posMid=ray_origin+tMid*vec_LOS
        call sph_locate_cell_fast(posMid,rface2,theta_cos,phiface,ixO^L,ix1,ix2,ix3,inside)
        if (.not. inside) cycle
        call collect_EUV_sph_dda_interval(ixI^L,ixO^L,source,opacity,pixel_id,&
             ray_origin,ximg1,ximg2,rface2,theta_cos,phiface,phi_sin,phi_cos,&
             t0,t1,segments,nseg,capacity,ok)
        if (.not. ok) then
          nseg=nsegStart
          fallback=.true.
          call collect_EUV_sph_intersection_segments(ixI^L,ixO^L,source,opacity,&
               pixel_id,ray_origin,ximg1,ximg2,rface,thetaface,phiface,rface2,&
               theta_cos,phi_sin,phi_cos,segments,nseg,capacity)
          return
        endif
      enddo

    end subroutine collect_EUV_sph_dda_segments

    subroutine integrate_EUV_sph_intersection_thin(numXI1,numXI2,xI1,xI2,dxI,fl,EM)
      use mod_global_parameters

      integer, intent(in) :: numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2),dxI
      type(te_fluid), intent(in) :: fl
      double precision, intent(inout) :: EM(numXI1,numXI2)

      integer :: ixO^L,ixI^L,ix^D
      integer :: iigrid,igrid,ixP1,ixP2,ixPmin1,ixPmax1,ixPmin2,ixPmax2
      integer :: iseg,nseg,capacity,sphDdaFallbackLocal,sphDdaFallbackGlobal
      double precision :: ray_origin(1:3),profile_local(3),profile_global(3)
      double precision, allocatable :: source(:^D&)
      double precision, allocatable :: rface(:),thetaface(:),phiface(:)
      double precision, allocatable :: rface2(:),theta_cos(:),phi_sin(:),phi_cos(:)
      double precision, allocatable :: segments(:,:)
      logical :: has_pixels,ddaFallback

      profile_local=zero
      sphDdaFallbackLocal=0
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixOmax^D=ixmhi^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\

        allocate(source(ixI^S))
        call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,source)
        source(ixO^S)=source(ixO^S)/instrument_resolution_factor**2
        call build_sph_intersection_faces(ixI^L,ixO^L,ps(igrid)%x,ps(igrid)%dx,rface,thetaface,phiface)
        if (sph_use_dda) then
          allocate(rface2(ixOmin1:ixOmax1+1),theta_cos(ixOmin2:ixOmax2+1),&
                   phi_sin(ixOmin3:ixOmax3+1),phi_cos(ixOmin3:ixOmax3+1))
          rface2=rface**2
          theta_cos=cos(thetaface)
          phi_sin=sin(phiface)
          phi_cos=cos(phiface)
        endif
        call sph_block_pixel_range(rface,thetaface,phiface,ixO^L,numXI1,numXI2,xI1,xI2,dxI,&
                                   ixPmin1,ixPmax1,ixPmin2,ixPmax2,has_pixels)
        if (has_pixels) then
          if (sph_use_dda) capacity=0
          do ixP1=ixPmin1,ixPmax1
            do ixP2=ixPmin2,ixPmax2
              ray_origin=xI1(ixP1)*vec_xI1+xI2(ixP2)*vec_xI2
              profile_local(1)=profile_local(1)+one
              if (sph_use_dda) then
                nseg=0
                call collect_EUV_sph_dda_segments(ixI^L,ixO^L,source,source,&
                     1,ray_origin,xI1(ixP1),xI2(ixP2),rface,thetaface,phiface,&
                     rface2,theta_cos,phi_sin,phi_cos,segments,nseg,capacity,ddaFallback)
                if (ddaFallback) sphDdaFallbackLocal=sphDdaFallbackLocal+1
                do iseg=1,nseg
                  EM(ixP1,ixP2)=EM(ixP1,ixP2)+segments(3,iseg)
                enddo
              else
                call acc_EUV_sph_intersection(ixI^L,ixO^L,source,ray_origin,xI1(ixP1),xI2(ixP2),&
                                              rface,thetaface,phiface,EM(ixP1,ixP2))
              endif
            enddo
          enddo
          profile_local(2)=profile_local(2)+dble((ixPmax1-ixPmin1+1)*(ixPmax2-ixPmin2+1))
        endif
        if (allocated(segments)) deallocate(segments)
        profile_local(3)=profile_local(3)+one
        if (allocated(rface2)) deallocate(rface2)
        if (allocated(theta_cos)) deallocate(theta_cos)
        if (allocated(phi_sin)) deallocate(phi_sin)
        if (allocated(phi_cos)) deallocate(phi_cos)
        deallocate(source,rface,thetaface,phiface)
      enddo
      call MPI_ALLREDUCE(profile_local,profile_global,3,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      call MPI_ALLREDUCE(sphDdaFallbackLocal,sphDdaFallbackGlobal,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
      if (radsyn_verbose .and. mype==0) then
        if (sph_use_dda) then
          write(*,'(a,3(es12.5,1x))') ' sph_dda thin profile rays pixels blocks: ',profile_global
          write(*,'(a,i0)') ' sph_dda thin fallback rays: ',sphDdaFallbackGlobal
        else
          write(*,'(a,3(es12.5,1x))') ' sph_intersection thin profile rays pixels blocks: ',profile_global
        endif
      endif
    end subroutine integrate_EUV_sph_intersection_thin

    subroutine integrate_EUV_sph_intersection_thick(numXI1,numXI2,xI1,xI2,dxI,fl,EUV,Tau,EUVthin)
      use mod_global_parameters

      integer, intent(in) :: numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2),dxI
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: EUV(numXI1,numXI2),Tau(numXI1,numXI2),EUVthin(numXI1,numXI2)

      integer, parameter :: nSegVars=5
      integer :: ixO^L,ixI^L,ix^D
      integer :: iigrid,igrid,ixP1,ixP2,ipix,ipixStart,ipixEnd,nPixBatch,pixel_id
      integer :: nseg,capacity,totalCount,totalSeg,ipe,is,iseg,nidx,owner,isegDest,nsegBefore
      integer :: ixGlobal,iyGlobal,ixPmin1,ixPmax1,ixPmin2,ixPmax2,iFirst,iLast,iLocal
      integer :: nPixBatchTarget
      integer :: maxSegBatchTarget,maxSegCommTarget,maxNsegBatch,nPixTotal
      integer :: maxOwnerSegCount,maxOwnerSegCountLocal,segOffset,recvFill,totalRoundCount,totalRoundSeg
      integer :: sphDdaFallbackLocal,sphDdaFallbackGlobal
      integer, allocatable :: sendCounts(:),recvCounts(:),sendDispls(:),recvDispls(:)
      integer, allocatable :: roundSendCounts(:),roundRecvCounts(:)
      integer, allocatable :: roundSendDispls(:),roundRecvDispls(:)
      integer, allocatable :: ownerSegCounts(:),ownerOffsets(:),idx(:)
      integer, allocatable :: bucketCounts(:),bucketOffsets(:),bucketFill(:)
      double precision :: ray_origin(1:3),atten
      double precision :: profile_local(5),profile_global(5),profile_batch(5)
      double precision :: phys_max_local(2),phys_max_global(2)
      double precision :: phys_sum_local(2),phys_sum_global(2),phys_sum_batch(2)
      double precision, allocatable :: segments(:,:),segments_send(:,:),segments_recv(:,:)
      double precision, allocatable :: segments_recv_round(:,:)
      double precision, allocatable :: image_reduce(:,:)
      logical :: has_pixels,batchAccepted,batchReduced,ddaFallback
      type(radsyn_euv_cache), allocatable :: cache(:)

      EUV=zero
      Tau=zero
      EUVthin=zero
      profile_local=zero
      phys_max_local=zero
      phys_sum_local=zero
      sphDdaFallbackLocal=0
      allocate(sendCounts(0:npe-1),recvCounts(0:npe-1),sendDispls(0:npe-1),recvDispls(0:npe-1))
      allocate(roundSendCounts(0:npe-1),roundRecvCounts(0:npe-1))
      allocate(roundSendDispls(0:npe-1),roundRecvDispls(0:npe-1))
      allocate(ownerSegCounts(0:npe-1),ownerOffsets(0:npe-1))
      allocate(cache(igridstail))
      call radsyn_get_segment_batch_limits(nPixBatchTarget,maxSegBatchTarget,maxSegCommTarget)
      allocate(bucketCounts(nPixBatchTarget),bucketOffsets(nPixBatchTarget+1),&
               bucketFill(nPixBatchTarget))

      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixOmax^D=ixmhi^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\

        cache(iigrid)%igrid=igrid
        allocate(cache(iigrid)%source(ixI^S),cache(iigrid)%opacity(ixI^S))
        cache(iigrid)%source=zero
        cache(iigrid)%opacity=zero
        call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,cache(iigrid)%source)
        cache(iigrid)%source(ixO^S)=cache(iigrid)%source(ixO^S)/instrument_resolution_factor**2
        call get_EUV_HHe_opacity(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,cache(iigrid)%opacity)
        phys_max_local(1)=max(phys_max_local(1),maxval(cache(iigrid)%source(ixO^S)))
        phys_max_local(2)=max(phys_max_local(2),maxval(cache(iigrid)%opacity(ixO^S)))
        call build_sph_intersection_faces(ixI^L,ixO^L,ps(igrid)%x,ps(igrid)%dx,&
                                          cache(iigrid)%rface,cache(iigrid)%thetaface,&
                                          cache(iigrid)%phiface)
        allocate(cache(iigrid)%rface2(ixOmin1:ixOmax1+1),&
                 cache(iigrid)%theta_cos(ixOmin2:ixOmax2+1),&
                 cache(iigrid)%phi_sin(ixOmin3:ixOmax3+1),&
                 cache(iigrid)%phi_cos(ixOmin3:ixOmax3+1))
        cache(iigrid)%rface2=cache(iigrid)%rface**2
        cache(iigrid)%theta_cos=cos(cache(iigrid)%thetaface)
        cache(iigrid)%phi_sin=sin(cache(iigrid)%phiface)
        cache(iigrid)%phi_cos=cos(cache(iigrid)%phiface)
        call sph_block_pixel_range(cache(iigrid)%rface,cache(iigrid)%thetaface,&
             cache(iigrid)%phiface,ixO^L,numXI1,numXI2,xI1,xI2,dxI,&
             cache(iigrid)%ixPmin1,cache(iigrid)%ixPmax1,&
             cache(iigrid)%ixPmin2,cache(iigrid)%ixPmax2,cache(iigrid)%has_pixels)
      enddo

      nPixTotal=numXI1*numXI2
      ipixStart=1
      do while (ipixStart<=nPixTotal)
        ipixEnd=min(numXI1*numXI2,ipixStart+nPixBatchTarget-1)
        nPixBatch=ipixEnd-ipixStart+1
        batchAccepted=.false.
        batchReduced=.false.

        do while (.not. batchAccepted)
          nseg=0
          capacity=0
          profile_batch=zero
          phys_sum_batch=zero

          do iigrid=1,igridstail; igrid=igrids(iigrid);
            ^D&ixOmin^D=ixmlo^D\
            ^D&ixOmax^D=ixmhi^D\
            ^D&ixImin^D=ixglo^D\
            ^D&ixImax^D=ixghi^D\

            ixPmin1=cache(iigrid)%ixPmin1
            ixPmax1=cache(iigrid)%ixPmax1
            ixPmin2=cache(iigrid)%ixPmin2
            ixPmax2=cache(iigrid)%ixPmax2
            has_pixels=cache(iigrid)%has_pixels
            if (.not. has_pixels) cycle

            do ixP2=ixPmin2,ixPmax2
              iFirst=max(ipixStart,(ixP2-1)*numXI1+ixPmin1)
              iLast=min(ipixEnd,(ixP2-1)*numXI1+ixPmax1)
              if (iFirst>iLast) cycle
              do ipix=iFirst,iLast
                ixP1=1+mod(ipix-1,numXI1)
                pixel_id=ipix
                ray_origin=xI1(ixP1)*vec_xI1+xI2(ixP2)*vec_xI2
                profile_batch(1)=profile_batch(1)+one
                nsegBefore=nseg
                if (sph_use_dda) then
                  call collect_EUV_sph_dda_segments(ixI^L,ixO^L,cache(iigrid)%source,&
                       cache(iigrid)%opacity,pixel_id,ray_origin,xI1(ixP1),xI2(ixP2),&
                       cache(iigrid)%rface,cache(iigrid)%thetaface,cache(iigrid)%phiface,&
                       cache(iigrid)%rface2,cache(iigrid)%theta_cos,&
                       cache(iigrid)%phi_sin,cache(iigrid)%phi_cos,&
                       segments,nseg,capacity,ddaFallback)
                  if (ddaFallback) sphDdaFallbackLocal=sphDdaFallbackLocal+1
                else
                  call collect_EUV_sph_intersection_segments(ixI^L,ixO^L,cache(iigrid)%source,&
                       cache(iigrid)%opacity,pixel_id,ray_origin,xI1(ixP1),xI2(ixP2),&
                       cache(iigrid)%rface,cache(iigrid)%thetaface,cache(iigrid)%phiface,&
                       cache(iigrid)%rface2,cache(iigrid)%theta_cos,&
                       cache(iigrid)%phi_sin,cache(iigrid)%phi_cos,&
                       segments,nseg,capacity)
                endif
                if (nseg>nsegBefore) profile_batch(2)=profile_batch(2)+one
                profile_batch(3)=profile_batch(3)+dble(nseg-nsegBefore)
                do iseg=nsegBefore+1,nseg
                  phys_sum_batch(1)=phys_sum_batch(1)+segments(3,iseg)
                  phys_sum_batch(2)=phys_sum_batch(2)+segments(4,iseg)
                enddo
              enddo
            enddo
          enddo

          call MPI_ALLREDUCE(nseg,maxNsegBatch,1,MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
          if (maxNsegBatch>maxSegBatchTarget .and. nPixBatch>1) then
            nPixBatch=max(1,nPixBatch/2)
            ipixEnd=ipixStart+nPixBatch-1
            if (allocated(segments)) deallocate(segments)
            batchReduced=.true.
          else
            batchAccepted=.true.
          endif
        enddo

        profile_local=profile_local+profile_batch
        phys_sum_local=phys_sum_local+phys_sum_batch
        if (radsyn_verbose .and. mype==0 .and. batchReduced) then
          if (sph_use_dda) then
            write(*,'(a,3(i0,1x))') ' sph_dda thick adaptive batch: ',&
                 ipixStart,ipixEnd,maxNsegBatch
          else
            write(*,'(a,3(i0,1x))') ' sph_intersection thick adaptive batch: ',&
                 ipixStart,ipixEnd,maxNsegBatch
          endif
        endif

        if (.not. allocated(segments)) then
          capacity=1
          allocate(segments(nSegVars,capacity))
        endif
        ownerSegCounts=0
        do is=1,nseg
          owner=segment_pixel_owner(nint(segments(1,is)))
          ownerSegCounts(owner)=ownerSegCounts(owner)+1
        enddo
        sendCounts=nSegVars*ownerSegCounts
        sendDispls(0)=0
        do ipe=1,npe-1
          sendDispls(ipe)=sendDispls(ipe-1)+sendCounts(ipe-1)
        enddo

        allocate(segments_send(nSegVars,max(1,nseg)))
        ownerOffsets=0
        do is=1,nseg
          owner=segment_pixel_owner(nint(segments(1,is)))
          isegDest=sendDispls(owner)/nSegVars+ownerOffsets(owner)+1
          segments_send(:,isegDest)=segments(:,is)
          ownerOffsets(owner)=ownerOffsets(owner)+1
        enddo

        call MPI_ALLTOALL(sendCounts,1,MPI_INTEGER,recvCounts,1,MPI_INTEGER,icomm,ierrmpi)
        recvDispls(0)=0
        do ipe=1,npe-1
          recvDispls(ipe)=recvDispls(ipe-1)+recvCounts(ipe-1)
        enddo
        totalCount=sum(recvCounts)
        totalSeg=totalCount/nSegVars
        profile_local(4)=profile_local(4)+dble(totalCount)
        allocate(segments_recv(nSegVars,max(1,totalSeg)))

        recvFill=0
        maxOwnerSegCountLocal=maxval(ownerSegCounts)
        call MPI_ALLREDUCE(maxOwnerSegCountLocal,maxOwnerSegCount,1,MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
        do segOffset=0,maxOwnerSegCount-1,maxSegCommTarget
          roundSendCounts=0
          roundSendDispls=sendDispls
          do ipe=0,npe-1
            if (ownerSegCounts(ipe)>segOffset) then
              roundSendCounts(ipe)=nSegVars*min(maxSegCommTarget,ownerSegCounts(ipe)-segOffset)
              roundSendDispls(ipe)=sendDispls(ipe)+nSegVars*segOffset
            endif
          enddo

          call MPI_ALLTOALL(roundSendCounts,1,MPI_INTEGER,roundRecvCounts,1,MPI_INTEGER,icomm,ierrmpi)
          roundRecvDispls(0)=0
          do ipe=1,npe-1
            roundRecvDispls(ipe)=roundRecvDispls(ipe-1)+roundRecvCounts(ipe-1)
          enddo
          totalRoundCount=sum(roundRecvCounts)
          totalRoundSeg=totalRoundCount/nSegVars
          allocate(segments_recv_round(nSegVars,max(1,totalRoundSeg)))

          call MPI_ALLTOALLV(segments_send,roundSendCounts,roundSendDispls,MPI_DOUBLE_PRECISION,&
                             segments_recv_round,roundRecvCounts,roundRecvDispls,&
                             MPI_DOUBLE_PRECISION,icomm,ierrmpi)

          if (totalRoundSeg>0) then
            segments_recv(:,recvFill+1:recvFill+totalRoundSeg)=segments_recv_round(:,1:totalRoundSeg)
            recvFill=recvFill+totalRoundSeg
          endif
          deallocate(segments_recv_round)
        enddo

        if (recvFill/=totalSeg) call mpistop("ray-segment receive mismatch")

        if (totalSeg>0) then
          allocate(idx(totalSeg))
          bucketCounts(1:nPixBatch)=0
          do is=1,totalSeg
            if (segment_is_valid(segments_recv,is,4)) then
              ipix=nint(segments_recv(1,is))
              if (ipix>=ipixStart .and. ipix<=ipixEnd .and. segment_pixel_owner(ipix)==mype) then
                iLocal=ipix-ipixStart+1
                bucketCounts(iLocal)=bucketCounts(iLocal)+1
              endif
            endif
          enddo

          bucketOffsets(1)=1
          do iLocal=1,nPixBatch
            bucketOffsets(iLocal+1)=bucketOffsets(iLocal)+bucketCounts(iLocal)
          enddo
          bucketFill(1:nPixBatch)=bucketOffsets(1:nPixBatch)
          do is=1,totalSeg
            if (segment_is_valid(segments_recv,is,4)) then
              ipix=nint(segments_recv(1,is))
              if (ipix>=ipixStart .and. ipix<=ipixEnd .and. segment_pixel_owner(ipix)==mype) then
                iLocal=ipix-ipixStart+1
                idx(bucketFill(iLocal))=is
                bucketFill(iLocal)=bucketFill(iLocal)+1
              endif
            endif
          enddo

          do ipix=ipixStart,ipixEnd
            if (segment_pixel_owner(ipix)/=mype) cycle
            iLocal=ipix-ipixStart+1
            nidx=bucketCounts(iLocal)
            if (nidx>0) then
              profile_local(5)=profile_local(5)+dble(nidx)*dble(nidx)
              call sort_segment_indices_near_to_far(segments_recv,&
                   idx(bucketOffsets(iLocal):bucketOffsets(iLocal+1)-1),nidx)
              ixGlobal=1+mod(ipix-1,numXI1)
              iyGlobal=1+(ipix-1)/numXI1
              do iseg=bucketOffsets(iLocal),bucketOffsets(iLocal+1)-1
                is=idx(iseg)
                EUVthin(ixGlobal,iyGlobal)=EUVthin(ixGlobal,iyGlobal)+segments_recv(3,is)
                atten=transfer_attenuation(Tau(ixGlobal,iyGlobal))
                EUV(ixGlobal,iyGlobal)=EUV(ixGlobal,iyGlobal)+atten*segments_recv(3,is)
                Tau(ixGlobal,iyGlobal)=Tau(ixGlobal,iyGlobal)+max(zero,segments_recv(4,is))
              enddo
            endif
          enddo
          deallocate(idx)
        endif

        deallocate(segments_send,segments_recv)
        if (allocated(segments)) deallocate(segments)
        ipixStart=ipixEnd+1
      enddo

      do iigrid=1,igridstail
        if (allocated(cache(iigrid)%source)) deallocate(cache(iigrid)%source)
        if (allocated(cache(iigrid)%opacity)) deallocate(cache(iigrid)%opacity)
        if (allocated(cache(iigrid)%rface)) deallocate(cache(iigrid)%rface)
        if (allocated(cache(iigrid)%thetaface)) deallocate(cache(iigrid)%thetaface)
        if (allocated(cache(iigrid)%phiface)) deallocate(cache(iigrid)%phiface)
        if (allocated(cache(iigrid)%rface2)) deallocate(cache(iigrid)%rface2)
        if (allocated(cache(iigrid)%theta_cos)) deallocate(cache(iigrid)%theta_cos)
        if (allocated(cache(iigrid)%phi_sin)) deallocate(cache(iigrid)%phi_sin)
        if (allocated(cache(iigrid)%phi_cos)) deallocate(cache(iigrid)%phi_cos)
      enddo
      deallocate(cache)
      deallocate(sendCounts,recvCounts,sendDispls,recvDispls,roundSendCounts,roundRecvCounts,&
                 roundSendDispls,roundRecvDispls,ownerSegCounts,ownerOffsets,bucketCounts,&
                 bucketOffsets,bucketFill)
      allocate(image_reduce(numXI1,numXI2))
      call MPI_ALLREDUCE(EUV,image_reduce,numXI1*numXI2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      EUV=image_reduce
      call MPI_ALLREDUCE(Tau,image_reduce,numXI1*numXI2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      Tau=image_reduce
      call MPI_ALLREDUCE(EUVthin,image_reduce,numXI1*numXI2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      EUVthin=image_reduce
      deallocate(image_reduce)
      call MPI_ALLREDUCE(profile_local,profile_global,5,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      call MPI_ALLREDUCE(phys_max_local,phys_max_global,2,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
      call MPI_ALLREDUCE(phys_sum_local,phys_sum_global,2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      call MPI_ALLREDUCE(sphDdaFallbackLocal,sphDdaFallbackGlobal,1,MPI_INTEGER,MPI_SUM,icomm,ierrmpi)
      if (radsyn_verbose .and. mype==0) then
        if (sph_use_dda) then
          write(*,'(a,5(es12.5,1x))') ' sph_dda thick profile: ',profile_global
          write(*,'(a,i0)') ' sph_dda thick fallback rays: ',sphDdaFallbackGlobal
          write(*,'(a,4(es12.5,1x))') ' sph_dda thick physics maxj maxk sumjds sumkds: ',&
               phys_max_global(1),phys_max_global(2),phys_sum_global(1),phys_sum_global(2)
        else
          write(*,'(a,5(es12.5,1x))') ' sph_intersection thick profile: ',profile_global
          write(*,'(a,4(es12.5,1x))') ' sph_intersection thick physics maxj maxk sumjds sumkds: ',&
               phys_max_global(1),phys_max_global(2),phys_sum_global(1),phys_sum_global(2)
        endif
      endif
    end subroutine integrate_EUV_sph_intersection_thick

  }

  {^IFTHREED

    subroutine get_sph_intersection_image_bounds(xImin1,xImax1,xImin2,xImax2)
      double precision, intent(out) :: xImin1,xImax1,xImin2,xImax2

      integer, parameter :: nsample=5
      integer :: iigrid,igrid,ir,it,ip
      integer :: ixI^L,ixO^L
      double precision, allocatable :: rface(:),thetaface(:),phiface(:)
      double precision :: local_min1,local_max1,local_min2,local_max2
      double precision :: sph(1:3),xcent(1:2),wr,wt,wp

      local_min1=huge(one)
      local_max1=-huge(one)
      local_min2=huge(one)
      local_max2=-huge(one)

      do iigrid=1,igridstail
        igrid=igrids(iigrid)
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixOmax^D=ixmhi^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\

        call build_sph_intersection_faces(ixI^L,ixO^L,ps(igrid)%x,ps(igrid)%dx,&
                                          rface,thetaface,phiface)
        do ir=0,nsample-1
          wr=dble(ir)/dble(nsample-1)
          sph(1)=(one-wr)*rface(ixOmin1)+wr*rface(ixOmax1+1)
          do it=0,nsample-1
            wt=dble(it)/dble(nsample-1)
            sph(2)=(one-wt)*thetaface(ixOmin2)+wt*thetaface(ixOmax2+1)
            do ip=0,nsample-1
              wp=dble(ip)/dble(nsample-1)
              if (ir/=0 .and. ir/=nsample-1 .and. it/=0 .and. it/=nsample-1 .and. &
                  ip/=0 .and. ip/=nsample-1) cycle
              sph(3)=(one-wp)*phiface(ixOmin3)+wp*phiface(ixOmax3+1)
              call get_cor_image_spherical(sph,xcent)
              local_min1=min(local_min1,xcent(1))
              local_max1=max(local_max1,xcent(1))
              local_min2=min(local_min2,xcent(2))
              local_max2=max(local_max2,xcent(2))
            enddo
          enddo
        enddo
        deallocate(rface,thetaface,phiface)
      enddo

      call MPI_ALLREDUCE(local_min1,xImin1,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      call MPI_ALLREDUCE(local_max1,xImax1,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
      call MPI_ALLREDUCE(local_min2,xImin2,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      call MPI_ALLREDUCE(local_max2,xImax2,1,MPI_DOUBLE_PRECISION,MPI_MAX,icomm,ierrmpi)
      if (xImin1>0.5d0*huge(one) .or. xImax1<-0.5d0*huge(one) .or. &
          xImin2>0.5d0*huge(one) .or. xImax2<-0.5d0*huge(one)) then
        call mpistop("sph_intersection could not determine image bounds")
      endif
    end subroutine get_sph_intersection_image_bounds

    subroutine get_sph_intersection_datresol_spacing(dxI)
      double precision, intent(out) :: dxI

      integer :: iigrid,igrid,ixI^L,ixO^L,ix^D
      double precision :: local_min,global_min,dr,ds_theta,ds_phi,rval,theta

      local_min=huge(one)
      do iigrid=1,igridstail
        igrid=igrids(iigrid)
        ^D&ixOmin^D=ixmlo^D\
        ^D&ixOmax^D=ixmhi^D\
        ^D&ixImin^D=ixglo^D\
        ^D&ixImax^D=ixghi^D\

        do ix1=ixOmin1,ixOmax1
          do ix2=ixOmin2,ixOmax2
            do ix3=ixOmin3,ixOmax3
              rval=max(smalldouble,ps(igrid)%x(ix^D,1))
              theta=ps(igrid)%x(ix^D,2)
              dr=ps(igrid)%dx(ix^D,1)
              ds_theta=rval*ps(igrid)%dx(ix^D,2)
              ds_phi=rval*max(smalldouble,sin(theta))*ps(igrid)%dx(ix^D,3)
              local_min=min(local_min,dr,ds_theta,ds_phi)
            enddo
          enddo
        enddo
      enddo

      call MPI_ALLREDUCE(local_min,global_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      if (global_min<=zero .or. global_min>half*huge(one)) then
        call mpistop("sph_intersection could not determine dat-resolution image spacing")
      endif
      dxI=global_min
    end subroutine get_sph_intersection_datresol_spacing

    subroutine get_image(qunit,datatype,fl)
      ! integrate emission flux along line of sight (LOS) 
      ! in a 3D simulation box and get a 2D EUV image
      use mod_global_parameters
      use mod_constants

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20), intent(in) :: datatype

      integer :: ix^D,numXI1,numXI2,numWI
      double precision :: xImin1,xImax1,xImin2,xImax2,xIcent1,xIcent2,dxI
      double precision, allocatable :: xI1(:),xI2(:),dxI1(:),dxI2(:)
      double precision, allocatable :: wI(:,:,:),wIs(:,:,:),EM(:,:),Dpl(:,:),Tau(:,:),EMthin(:,:),WLB(:,:,:)
      double precision :: vec_temp1(1:3),vec_temp2(1:3)
      double precision :: vec_z(1:3),vec_cor(1:3),xI_cor(1:2)
      double precision :: res,LOS_psi,r_max,r_loc

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl,wslit
      double precision :: arcsec,RHESSI_rsl,LASCO_rsl,pixel,R_occult,smallflux
      integer :: iigrid,igrid,i,j,numSI,iw
      logical :: emit,ray_image_global,has_thick_output

      if (coordinate==spherical) then
        call init_vectors_spherical()
      else
        ! cartesian
        call init_vectors_cartesian()
      endif

      ! calculate domain of the image
      if (coordinate==spherical) then
        if (trim(ray_method_active)=='spherical' .and. datatype=='image_euv') then
          call get_sph_intersection_image_bounds(xImin1,xImax1,xImin2,xImax2)
        else
          xImin1=-abs(xprobmax1)
          xImin2=-abs(xprobmax1)
          xImax1=abs(xprobmax1)
          xImax2=abs(xprobmax1)
        endif
      else
        ! calculate domain of the image
        do ix1=1,2
          if (ix1==1) vec_cor(1)=xprobmin1
          if (ix1==2) vec_cor(1)=xprobmax1
          do ix2=1,2
            if (ix2==1) vec_cor(2)=xprobmin2
            if (ix2==2) vec_cor(2)=xprobmax2
            do ix3=1,2
              if (ix3==1) vec_cor(3)=xprobmin3
              if (ix3==2) vec_cor(3)=xprobmax3
              if (big_image) then
                r_loc=(vec_cor(1)-x_origin(1))**2
                r_loc=r_loc+(vec_cor(2)-x_origin(2))**2
                r_loc=r_loc+(vec_cor(3)-x_origin(3))**2
                r_loc=sqrt(r_loc)
                if (ix1==1 .and. ix2==1 .and. ix3==1) then
                  r_max=r_loc
                else
                  r_max=max(r_max,r_loc)
                endif
              else
                call get_cor_image(vec_cor,xI_cor)
                if (ix1==1 .and. ix2==1 .and. ix3==1) then
                  xImin1=xI_cor(1)
                  xImax1=xI_cor(1)
                  xImin2=xI_cor(2)
                  xImax2=xI_cor(2)
                else
                  xImin1=min(xImin1,xI_cor(1))
                  xImax1=max(xImax1,xI_cor(1))
                  xImin2=min(xImin2,xI_cor(2))
                  xImax2=max(xImax2,xI_cor(2))
                endif
              endif
            enddo
          enddo
        enddo
        if (big_image) then
          xImin1=-r_max
          xImin2=-r_max
          xImax1=r_max
          xImax2=r_max
        endif
      endif
      xIcent1=(xImin1+xImax1)/2.d0
      xIcent2=(xImin2+xImax2)/2.d0

      ! tables for image
      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif
      if (datatype=='image_euv') then
        call get_line_info(wavelength,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
        dxI=spaceRsl*arcsec  ! intrument resolution of image
        smallflux=smalldouble
      else if (datatype=='image_sxr') then
        RHESSI_rsl=2.3d0/instrument_resolution_factor
        dxI=RHESSI_rsl*arcsec
        smallflux=1.d-40
      else if (datatype=='image_whitelight') then
        if (whitelight_instrument=='LASCO/C1') then
          LASCO_rsl=5.6d0/instrument_resolution_factor
          R_occult=1.1d0
        else if (whitelight_instrument=='LASCO/C2') then
          LASCO_rsl=11.4d0/instrument_resolution_factor
          R_occult=2.d0
        else if (whitelight_instrument=='LASCO/C3') then
          LASCO_rsl=56.d0/instrument_resolution_factor
          R_occult=3.7d0
        else
          call MPISTOP('Whitelight synthesis: instrument is not supported!')
        endif
        dxI=LASCO_rsl*arcsec
        if (R_occultor>1.d0) R_occult=R_occultor
        R_occult=R_occult*const_Rsun/unit_length
        smallflux=1.d-20
      endif
      numXI1=8*ceiling((xImax1-xIcent1)/dxI/8.d0)
      xImin1=xIcent1-numXI1*dxI
      xImax1=xIcent1+numXI1*dxI
      numXI1=numXI1*2
      numXI2=8*ceiling((xImax2-xIcent2)/dxI/8.d0)
      xImin2=xIcent2-numXI2*dxI
      xImax2=xIcent2+numXI2*dxI
      numXI2=numXI2*2
      allocate(xI1(numXI1),xI2(numXI2),dxI1(numXI1),dxI2(numXI2))
      do ix1=1,numXI1
        xI1(ix1)=xImin1+dxI*(ix1-half)
        dxI1(ix1)=dxI
      enddo
      do ix2=1,numXI2
        xI2(ix2)=xImin2+dxI*(ix2-half)
        dxI2(ix2)=dxI
      enddo

      ! calculate emission
      if (datatype=='image_euv' .or. datatype=='image_sxr') then
        has_thick_output=datatype=='image_euv' .and. trim(radiation_transfer)=='thick' .and. &
            ((coordinate==spherical .and. trim(ray_method_active)=='spherical') .or. &
             (coordinate==cartesian .and. trim(ray_method_active)=='cart'))
        if (datatype=='image_euv') then
          numWI=radsyn_euv_num_outputs(.false.,has_thick_output)
        else
          numWI=1
        endif
        allocate(wI(numXI1,numXI2,numWI),wIs(numXI1,numXI2,numWI),EM(numXI1,numXI2))
        wI=zero
        wIs=zero
        EM=zero
        ray_image_global=.false.
        if (has_thick_output) then
          allocate(Tau(numXI1,numXI2),EMthin(numXI1,numXI2))
          Tau=zero
          EMthin=zero
        endif
        if (coordinate==cartesian .and. datatype=='image_euv' .and. &
            trim(ray_method_active)=='cart') then
          ray_image_global=.true.
          allocate(Dpl(numXI1,numXI2))
          Dpl=zero
          if (trim(radiation_transfer)=='thick') then
            call integrate_EUV_cart_dda_thick_datresol(numXI1,numXI2,xI1,xI2,fl,EM,Dpl,Tau,EMthin)
          else
            call integrate_EUV_cart_dda_datresol(numXI1,numXI2,xI1,xI2,fl,EM,Dpl)
          endif
          deallocate(Dpl)
        else if (coordinate==cartesian) then
          do iigrid=1,igridstail; igrid=igrids(iigrid);
            call integrate_emission_cartesian(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,datatype,EM)
          enddo
        else if (trim(ray_method_active) == 'spherical' .and. datatype == 'image_euv') then
          if (trim(radiation_transfer) == 'thick') then
            ray_image_global=.true.
            call integrate_EUV_sph_intersection_thick(numXI1,numXI2,xI1,xI2,dxI,fl,EM,Tau,EMthin)
          else
            call integrate_EUV_sph_intersection_thin(numXI1,numXI2,xI1,xI2,dxI,fl,EM)
          endif
        else
          do iigrid=1,igridstail; igrid=igrids(iigrid);
            call integrate_emission_spherical(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,datatype,EM)
          enddo
        endif
        if (ray_image_global) then
          if (has_thick_output) then
            call pack_euv_image_outputs(numXI1,numXI2,EM,wI,smallflux,.false.,&
                                        has_thick_output,Tau=Tau,EUVthin=EMthin,&
                                        cap_absorption=.true.)
          else
            call pack_euv_image_outputs(numXI1,numXI2,EM,wI,smallflux,.false.,has_thick_output)
          endif
        else
          do ix1=1,numXI1
            do ix2=1,numXI2
              if (EM(ix1,ix2)>smallflux) wIs(ix1,ix2,1)=EM(ix1,ix2)
            enddo
          enddo
        endif
        if (.not. ray_image_global) then
          numSI=numXI1*numXI2*numWI
          call MPI_ALLREDUCE(wIs,wI,numSI,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
        endif
        if (activate_unit_arcsec) then
          xI1=xI1/arcsec
          dxI1=dxI1/arcsec
          xI2=xI2/arcsec
          dxI2=dxI2/arcsec
        endif
        call output_data(qunit,xI1,xI2,dxI1,dxI2,wI,numXI1,numXI2,numWI,datatype)
        if (allocated(Tau)) deallocate(Tau)
        if (allocated(EMthin)) deallocate(EMthin)
        deallocate(wI,wIs,EM)
      else if (datatype=='image_whitelight') then
        numWI=2
        allocate(wI(numXI1,numXI2,numWI),wIs(numXI1,numXI2,numWI),WLB(numXI1,numXI2,numWI))
        wI=zero
        wIs=zero
        WLB=zero
        if (coordinate==spherical) then
          do iigrid=1,igridstail; igrid=igrids(iigrid);
            call integrate_whitelight_spherical(igrid,numXI1,numXI2,numWI,xI1,xI2,dxI,fl,datatype,WLB)
          enddo
        endif
        do ix1=1,numXI1
          do ix2=1,numXI2
            if (WLB(ix1,ix2,1)>smallflux) then 
              wIs(ix1,ix2,1)=WLB(ix1,ix2,1)
              wIs(ix1,ix2,2)=WLB(ix1,ix2,2)
            endif
          enddo
        enddo
        numSI=numXI1*numXI2*numWI
        call MPI_ALLREDUCE(wIs,wI,numSI,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
        if (activate_unit_arcsec) then
          xI1=xI1/arcsec
          dxI1=dxI1/arcsec
          xI2=xI2/arcsec
          dxI2=dxI2/arcsec
        endif
        call output_data(qunit,xI1,xI2,dxI1,dxI2,wI,numXI1,numXI2,numWI,datatype)
        deallocate(wI,wIs,WLB)
      endif

      deallocate(xI1,xI2,dxI1,dxI2)

    end subroutine get_image

    subroutine integrate_emission_cartesian(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,datatype,EM)
      integer, intent(in) :: igrid,numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2)
      double precision, intent(in) :: dxI
      type(te_fluid), intent(in) :: fl
      character(20), intent(in) :: datatype
      double precision, intent(inout) :: EM(numXI1,numXI2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&),opacity(:^D&)
      double precision :: res
      integer :: ixP^L,ixP^D,nSubC^D,iSubC^D
      double precision :: xSubP1,xSubP2,dxSubP,xerf^L,fluxsubC
      double precision :: xSubC(1:3),dxSubC^D,xCent(1:2)

      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: lineCent
      double precision :: sigma_PSF,spaceRsl,wlRsl,sigma0,factor,wslit
      double precision :: arcsec,pixel,RHESSI_rsl,area_1AU
      double precision :: aa,bb

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif

      allocate(flux(ixI^S),opacity(ixI^S))
      if (datatype=='image_euv') then
        if (trim(emission_model)=='pseudo_current') then
          call get_pseudo_current(igrid,ixI^L,ixO^L,ps(igrid)%w,flux)
        else if (trim(emission_model)=='radio_ff') then
          call get_radio_ff_source_opacity(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux,opacity)
        else
          ! get local EUV flux and velocity
          call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
          flux(ixO^S)=flux(ixO^S)/instrument_resolution_factor**2   ! adjust flux due to artifical change of resolution
        endif
        call get_line_info(wavelength,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
        pixel=spaceRsl*arcsec
        sigma0=sigma_PSF*pixel
      else if (datatype=='image_sxr') then
        ! get local SXR flux photons cm^-3 s^-1 (cgs) or photons m^-3 s^-1 (SI)
        call get_SXR(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux,emin_sxr,emax_sxr) 
        RHESSI_rsl=2.3d0/instrument_resolution_factor
        sigma_PSF=1.d0
        pixel=RHESSI_rsl*arcsec
        sigma0=sigma_PSF*pixel
        area_1AU=2.81d27
      endif

      ! integrate emission
      {do ix^D=ixOmin^D,ixOmax^D\}
        ^D&nSubC^D=1;
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI1(^D))/(dxI/2.d0)));
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI2(^D))/(dxI/2.d0)));
        ^D&dxSubC^D=ps(igrid)%dx(ix^DD,^D)/nSubC^D;
        if (datatype=='image_euv') then
          if (SI_unit) then
            fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length*1.d2/dxI/dxI  ! DN s^-1
          else
            fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length/dxI/dxI  ! DN s^-1
          endif
        else if (datatype=='image_sxr') then
          ! sub-cell SXR flux at 1 AU [photons s^-1 cm^-2]
          fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length**3/area_1AU
        endif
        if (fluxSubC>smalldouble) then
          ! dividing a cell to several parts to get more accurate integrating values
          {do iSubC^D=1,nSubC^D\}
            ^D&xSubC(^D)=ps(igrid)%x(ix^DD,^D)-half*ps(igrid)%dx(ix^DD,^D)+(iSubC^D-half)*dxSubC^D;
            ! mapping the 3D coordinate to location at the image
            call get_cor_image(xSubC,xCent)
            ! distribution at nearby pixels
            ixP1=floor((xCent(1)-(xI1(1)-half*dxI))/dxI)+1
            ixP2=floor((xCent(2)-(xI2(1)-half*dxI))/dxI)+1
            ixPmin1=max(1,ixP1-3)
            ixPmax1=min(ixP1+3,numXI1)
            ixPmin2=max(1,ixP2-3)
            ixPmax2=min(ixP2+3,numXI2)
            do ixP1=ixPmin1,ixPmax1
              do ixP2=ixPmin2,ixPmax2
                xerfmin1=((xI1(ixP1)-half*dxI)-xCent(1))/(sqrt(2.d0)*sigma0) 
                xerfmax1=((xI1(ixP1)+half*dxI)-xCent(1))/(sqrt(2.d0)*sigma0)
                xerfmin2=((xI2(ixP2)-half*dxI)-xCent(2))/(sqrt(2.d0)*sigma0)
                xerfmax2=((xI2(ixP2)+half*dxI)-xCent(2))/(sqrt(2.d0)*sigma0)
                factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
                EM(ixP1,ixP2)=EM(ixP1,ixP2)+fluxSubC*factor
              enddo !ixP2
            enddo !ixP1
          {enddo\} !iSubC
        endif
      {enddo\} !ix

      deallocate(flux,opacity)
    end subroutine integrate_emission_cartesian

    subroutine integrate_emission_spherical(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,datatype,EM)
      integer, intent(in) :: igrid,numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2)
      double precision, intent(in) :: dxI
      type(te_fluid), intent(in) :: fl
      character(20), intent(in) :: datatype
      double precision, intent(inout) :: EM(numXI1,numXI2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision, allocatable :: flux(:^D&),Ne(:^D&),opacity(:^D&)
      integer :: ixP^L,ixP^D,nSubC^D,iSubC^D
      double precision :: xSubP1,xSubP2,dxSubP,xerf^L,fluxsubC,RsubC
      double precision :: TBsubC,PBsubC
      double precision :: xSubC(1:3),dxSubC^D,xCent(1:2),xSubC_car(1:3)
      double precision :: R_thick,dotp,dvolume,R_occult,Rc
      double precision :: dxl(1:3),x_sph(1:3),dx_sph(1:3)
      double precision :: unitv_r(1:3),unitv_theta(1:3),unitv_phi(1:3)
      logical :: sun_back_side,emit

      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: lineCent
      double precision :: sigma_PSF,spaceRsl,wlRsl,sigma0,factor,wslit
      double precision :: RHESSI_rsl,area_1AU,arcsec,pixel

      ^D&ixOmin^D=ixmlo^D;
      ^D&ixOmax^D=ixmhi^D;
      ^D&ixImin^D=ixglo^D;
      ^D&ixImax^D=ixghi^D;

      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif

      allocate(flux(ixI^S),opacity(ixI^S))
      if (datatype=='image_euv') then
        if (trim(emission_model)=='pseudo_current') then
          call get_pseudo_current(igrid,ixI^L,ixO^L,ps(igrid)%w,flux)
        else if (trim(emission_model)=='radio_ff') then
          call get_radio_ff_source_opacity(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux,opacity)
        else
          ! get local EUV flux and velocity
          call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
          flux(ixO^S)=flux(ixO^S)/instrument_resolution_factor**2   ! adjust flux due to artifical change of resolution
        endif
        call get_line_info(wavelength,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
        pixel=spaceRsl*arcsec
        sigma0=sigma_PSF*pixel
      else if (datatype=='image_sxr') then
        ! get local SXR flux photons cm^-3 s^-1 (cgs) or photons m^-3 s^-1 (SI)
        call get_SXR(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux,emin_sxr,emax_sxr) 
        RHESSI_rsl=2.3d0/instrument_resolution_factor
        sigma_PSF=1.d0
        pixel=RHESSI_rsl*arcsec
        sigma0=sigma_PSF*pixel
        area_1AU=2.81d27
      endif

      ! integrate emission
      R_thick=R_opt_thick*const_Rsun/unit_length
      {do ix^D=ixOmin^D,ixOmax^D\}
        x_sph(1:3)=ps(igrid)%x(ix^D,1:3)
        dx_sph(1:3)=ps(igrid)%dx(ix^D,1:3)
        dxl(1)=dx_sph(1) ! cell size in length
        dxl(2)=x_sph(1)*dx_sph(2) ! cell size in length
        dxl(3)=x_sph(1)*dsin(x_sph(2))*dx_sph(3) ! cell size in length
        ! dividing a cell to several sub-cells to get more accurate integrating values
        ^D&nSubC^D=1;
        call get_unit_vector_spherical(x_sph,unitv_r,unitv_theta,unitv_phi)
        call dot_product_loc(unitv_r,vec_xI1,dotp)
        nSubC1=max(nSubC1,ceiling(dxl(1)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_r,vec_xI2,dotp)
        nSubC1=max(nSubC1,ceiling(dxl(1)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_theta,vec_xI1,dotp)
        nSubC2=max(nSubC2,ceiling(dxl(2)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_theta,vec_xI2,dotp)
        nSubC2=max(nSubC2,ceiling(dxl(2)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_phi,vec_xI1,dotp)
        nSubC3=max(nSubC3,ceiling(dxl(3)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_phi,vec_xI2,dotp)
        nSubC3=max(nSubC3,ceiling(dxl(3)*abs(dotp)/(dxI/2.d0)))

        ! integrate sub-cells
        do iSubc1=1,nSubC1
          ! sub-cell center coordinate in spherical
          xSubC(1)=x_sph(1)-half*dx_sph(1)+(iSubC1-half)*dx_sph(1)/nSubC1
          RsubC=xSubC(1)
          dxSubC1=dx_sph(1)/nSubC1 ! sub-cell size in length
          do iSubc2=1,nSubC2
            ! sub-cell center coordinate in spherical
            xSubC(2)=x_sph(2)-half*dx_sph(2)+(iSubC2-half)*dx_sph(2)/nSubC2
            dxSubC2=xSubC(1)*dx_sph(2)/nSubC2 ! sub-cell size in length
            dxSubC3=xSubC(1)*dsin(xSubC(2))*dx_sph(3)/nSubC3  ! sub-cell size in length
            dvolume=dxSubC1*dxSubC2*dxSubC3
            if (datatype=='image_euv') then
              if (SI_unit) then
                fluxSubC=flux(ix^D)*dvolume*unit_length*1.d2/dxI/dxI  ! DN s^-1
              else
                fluxSubC=flux(ix^D)*dvolume*unit_length/dxI/dxI  ! DN s^-1
              endif
            else if (datatype=='image_sxr') then
              ! sub-cell SXR flux at 1 AU [photons s^-1 cm^-2]
              fluxSubC=flux(ix^D)*dvolume*unit_length**3/area_1AU
            endif
            ! enter integration if flux large enough
            if (fluxSubC>smalldouble) then
              do iSubc3=1,nSubC3
                ! sub-cell center coordinate in spherical
                xSubC(3)=x_sph(3)-half*dx_sph(3)+(iSubC3-half)*dx_sph(3)/nSubC3
                call get_cor_image_spherical(xSubC,xCent)
                Rc=dsqrt(xCent(1)**2+xCent(2)**2) ! distance to sun center (on the image plane)
                !
                ! whether the local emitted photons can arrive the telescope
                call spherical_to_cartesian(xSubC,xSubC_car)
                call dot_product_loc(vec_LOS,xSubC_car,dotp)
                sun_back_side=.true.
                if (dotp<0.d0) sun_back_side=.false.
                ! whether the local emission can reach the telescope
                if (sun_back_side) then
                  emit=.false.
                  if (Rc>R_thick) emit=.true.
                else
                  emit=.true.
                  if (xSubC(1)<=R_thick) emit=.false.
                endif
                !
                if (emit) then
                  ! mapping the 3D coordinate to location at the image
                  ! distribution at nearby pixels
                  ixP1=floor((xCent(1)-(xI1(1)-half*dxI))/dxI)+1
                  ixP2=floor((xCent(2)-(xI2(1)-half*dxI))/dxI)+1
                  ixPmin1=max(1,ixP1-3)
                  ixPmax1=min(ixP1+3,numXI1)
                  ixPmin2=max(1,ixP2-3)
                  ixPmax2=min(ixP2+3,numXI2)
                  do ixP1=ixPmin1,ixPmax1
                    do ixP2=ixPmin2,ixPmax2
                      xerfmin1=((xI1(ixP1)-half*dxI)-xCent(1))/(sqrt(2.d0)*sigma0)
                      xerfmax1=((xI1(ixP1)+half*dxI)-xCent(1))/(sqrt(2.d0)*sigma0)
                      xerfmin2=((xI2(ixP2)-half*dxI)-xCent(2))/(sqrt(2.d0)*sigma0)
                      xerfmax2=((xI2(ixP2)+half*dxI)-xCent(2))/(sqrt(2.d0)*sigma0)
                      factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
                      EM(ixP1,ixP2)=EM(ixP1,ixP2)+fluxSubC*factor
                    enddo !ixP2
                  enddo !ixP1
                endif !emit
              enddo !iSubC3
            endif !smallflux
          enddo !iSubC2
        enddo !iSubC1
      {enddo\} !ix

      deallocate(flux,opacity)

    end subroutine integrate_emission_spherical

    subroutine integrate_whitelight_spherical(igrid,numXI1,numXI2,numWI,xI1,xI2,dxI,fl,datatype,WLB)
      use mod_eos, only: eos

      integer, intent(in) :: igrid,numXI1,numXI2,numWI
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2)
      double precision, intent(in) :: dxI
      type(te_fluid), intent(in) :: fl
      character(20), intent(in) :: datatype
      double precision, intent(inout) :: WLB(numXI1,numXI2,numWI)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision, allocatable :: flux(:^D&),Ne(:^D&)
      integer :: ixP^L,ixP^D,nSubC^D,iSubC^D
      double precision :: xSubP1,xSubP2,dxSubP,xerf^L,fluxsubC,RsubC
      double precision :: sigma_PSF,sigma0,arcsec,pixel,LASCO_rsl
      double precision :: A,B,C,D,Rc,Ne0,TBsubC,PBsubC,factor
      double precision :: R_thick,dotp,dvolume,R_occult
      double precision :: xSubC(1:3),dxSubC^D,xCent(1:2),xSubC_car(1:3)
      double precision :: dxl(1:3),x_sph(1:3),dx_sph(1:3)
      double precision :: unitv_r(1:3),unitv_theta(1:3),unitv_phi(1:3)
      logical :: emit

      ^D&ixOmin^D=ixmlo^D;
      ^D&ixOmax^D=ixmhi^D;
      ^D&ixImin^D=ixglo^D;
      ^D&ixImax^D=ixghi^D;

      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif

      allocate(Ne(ixI^S))
      if (whitelight_instrument=='LASCO/C1') then
        LASCO_rsl=5.6d0/instrument_resolution_factor
        R_occult=1.1d0
      else if (whitelight_instrument=='LASCO/C2') then
        LASCO_rsl=11.4d0/instrument_resolution_factor
        R_occult=2.d0
      else if (whitelight_instrument=='LASCO/C3') then
        LASCO_rsl=56.d0/instrument_resolution_factor
        R_occult=3.7d0
      endif
      if (R_occultor>1.d0) R_occult=R_occultor
      R_occult=R_occult*const_Rsun/unit_length
      call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,Ne)
      ! get actual electron density from EoS (replaces rho with ne)
      block
        double precision :: nH_dummy(ixI^S)
        call eos%get_ne_nH(ixI^L, ixO^L, ps(igrid)%w, Ne, nH_dummy)
      end block
      sigma_PSF=1.d0
      pixel=LASCO_rsl*arcsec
      sigma0=sigma_PSF*pixel

      ! integrate emission
      R_thick=R_opt_thick*const_Rsun/unit_length
      {do ix^D=ixOmin^D,ixOmax^D\}
        x_sph(1:3)=ps(igrid)%x(ix^D,1:3)
        dx_sph(1:3)=ps(igrid)%dx(ix^D,1:3)
        dxl(1)=dx_sph(1) ! cell size in length
        dxl(2)=x_sph(1)*dx_sph(2) ! cell size in length
        dxl(3)=x_sph(1)*dsin(x_sph(2))*dx_sph(3) ! cell size in length
        Ne0=Ne(ix^D)*unit_numberdensity
        ! dividing a cell to several sub-cells to get more accurate integrating values
        ^D&nSubC^D=1;
        call get_unit_vector_spherical(x_sph,unitv_r,unitv_theta,unitv_phi)
        call dot_product_loc(unitv_r,vec_xI1,dotp)
        nSubC1=max(nSubC1,ceiling(dxl(1)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_r,vec_xI2,dotp)
        nSubC1=max(nSubC1,ceiling(dxl(1)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_theta,vec_xI1,dotp)
        nSubC2=max(nSubC2,ceiling(dxl(2)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_theta,vec_xI2,dotp)
        nSubC2=max(nSubC2,ceiling(dxl(2)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_phi,vec_xI1,dotp)
        nSubC3=max(nSubC3,ceiling(dxl(3)*abs(dotp)/(dxI/2.d0)))
        call dot_product_loc(unitv_phi,vec_xI2,dotp)
        nSubC3=max(nSubC3,ceiling(dxl(3)*abs(dotp)/(dxI/2.d0)))

        ! integrate sub-cells
        do iSubc1=1,nSubC1
          ! sub-cell center coordinate in spherical
          xSubC(1)=x_sph(1)-half*dx_sph(1)+(iSubC1-half)*dx_sph(1)/nSubC1
          RsubC=xSubC(1)
          dxSubC1=dx_sph(1)/nSubC1 ! sub-cell size in length
          call get_Thomson_parameters(RsubC,A,B,C,D)
          do iSubc2=1,nSubC2
            ! sub-cell center coordinate in spherical
            xSubC(2)=x_sph(2)-half*dx_sph(2)+(iSubC2-half)*dx_sph(2)/nSubC2
            dxSubC2=xSubC(1)*dx_sph(2)/nSubC2 ! sub-cell size in length
            dxSubC3=xSubC(1)*dsin(xSubC(2))*dx_sph(3)/nSubC3  ! sub-cell size in length
            dvolume=dxSubC1*dxSubC2*dxSubC3
            do iSubc3=1,nSubC3
              ! sub-cell center coordinate in spherical
              xSubC(3)=x_sph(3)-half*dx_sph(3)+(iSubC3-half)*dx_sph(3)/nSubC3
              call get_cor_image_spherical(xSubC,xCent)
              Rc=dsqrt(xCent(1)**2+xCent(2)**2) ! distance to sun center (on the image plane)
              ! whether the local emitted photons can arrive the telescope
              emit=.false.
              if (Rc>R_occult) then 
                emit=.true.
                ! scaterring flux from cm^-3 of plasma
                call get_whitelight_Thomson(RsubC,Rc,Ne0,A,B,C,D,TBsubC,PBsubC)
                TBsubC=TBsubC*dvolume*unit_length/dxI/dxI
                PBsubC=PBsubC*dvolume*unit_length/dxI/dxI
                if (TBsubC<1.d-20) emit=.false.
              endif
              if (emit) then
                ! mapping the 3D coordinate to location at the image
                ! distribution at nearby pixels
                ixP1=floor((xCent(1)-(xI1(1)-half*dxI))/dxI)+1
                ixP2=floor((xCent(2)-(xI2(1)-half*dxI))/dxI)+1
                ixPmin1=max(1,ixP1-3)
                ixPmax1=min(ixP1+3,numXI1)
                ixPmin2=max(1,ixP2-3)
                ixPmax2=min(ixP2+3,numXI2)
                do ixP1=ixPmin1,ixPmax1
                  do ixP2=ixPmin2,ixPmax2
                    xerfmin1=((xI1(ixP1)-half*dxI)-xCent(1))/(sqrt(2.d0)*sigma0)
                    xerfmax1=((xI1(ixP1)+half*dxI)-xCent(1))/(sqrt(2.d0)*sigma0)
                    xerfmin2=((xI2(ixP2)-half*dxI)-xCent(2))/(sqrt(2.d0)*sigma0)
                    xerfmax2=((xI2(ixP2)+half*dxI)-xCent(2))/(sqrt(2.d0)*sigma0)
                    factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
                    WLB(ixP1,ixP2,1)=WLB(ixP1,ixP2,1)+TBsubC*factor
                    WLB(ixP1,ixP2,2)=WLB(ixP1,ixP2,2)+PBsubC*factor
                  enddo !ixP2
                enddo !ixP1
              endif
            enddo !iSubC3
          enddo !iSubC2
        enddo !iSubC1
      {enddo\} !ix

      deallocate(Ne)

    end subroutine integrate_whitelight_spherical

    subroutine get_Thomson_parameters(Rl,A,B,C,D)
      ! parameters given in Billings 1968
      use mod_constants
      double precision, intent(in) :: Rl
      double precision, intent(inout) :: A,B,C,D

      double precision :: sinO,cosO,sinO2,cosO2,tmp

      sinO=const_Rsun/(Rl*unit_length)
      sinO2=sinO**2
      cosO2=1.d0-sinO2
      cosO=abs(dsqrt(cosO2))
      tmp=log((1.d0+sinO)/cosO)
      A=cosO*sinO2
      B=-(1.d0-3.d0*sinO2-(cosO2/sinO)*(1.d0+3.d0*sinO2)*tmp)/8.d0
      C=4.d0/3.d0-cosO-cosO*cosO2/3.d0
      D=(5.d0+sinO2-(cosO2/sinO)*(5.d0-sinO2)*tmp)/8.d0

    end subroutine get_Thomson_parameters

    subroutine get_whitelight_Thomson(Rl,Rin,Ne,A,B,C,D,fluxTB,fluxPB)
      ! use the method in SSW/eltheory
      double precision, intent(in) :: Rl,Rin,Ne,A,B,C,D
      double precision, intent(inout) :: fluxTB,fluxPB

      double precision :: const,u,Bt,Br,PB,TB,sinchi2

      u=0.63d0
      const=1.24878d-25/(1.d0-u/3.d0)
      sinchi2=(Rin/Rl)**2
      Bt=const*(C+u*(D-C))
      PB=const*sinchi2*((A+u*(B-A)))
      Br=Bt-pB
      TB=Bt+Br
      fluxTB=TB*Ne
      fluxPB=PB*Ne

    end subroutine get_whitelight_Thomson

    subroutine get_unit_vector_spherical(x_sph,unitv_r,unitv_theta,unitv_phi)
      double precision, intent(in) :: x_sph(1:3)
      double precision, intent(inout) :: unitv_r(1:3),unitv_theta(1:3),unitv_phi(1:3)

      unitv_r(1)=dsin(x_sph(2))*dcos(x_sph(3))
      unitv_r(2)=dsin(x_sph(2))*dsin(x_sph(3))
      unitv_r(3)=dcos(x_sph(2))
      unitv_theta(1)=dcos(x_sph(2))*dcos(x_sph(3))
      unitv_theta(2)=dcos(x_sph(2))*dsin(x_sph(3))
      unitv_theta(3)=-dsin(x_sph(2))
      unitv_phi(1)=-dsin(x_sph(3))
      unitv_phi(2)=dcos(x_sph(3))
      unitv_phi(3)=zero

    end subroutine get_unit_vector_spherical

    subroutine output_data(qunit,xO1,xO2,dxO1,dxO2,wO,nXO1,nXO2,nWO,datatype)
      ! change the format of data and write data
      use mod_global_parameters

      integer, intent(in) :: qunit,nXO1,nXO2,nWO
      double precision, intent(in) :: dxO1(nxO1),dxO2(nxO2)
      double precision, intent(in) :: xO1(nXO1),xO2(nxO2)
      double precision, intent(inout) :: wO(nXO1,nXO2,nWO)
      character(20), intent(in) :: datatype

      integer :: nPiece,nP1,nP2,nC1,nC2,nWC
      integer :: piece_nmax1,piece_nmax2,ix1,ix2,j,ipc,ixc1,ixc2
      double precision :: uniform_tol
      double precision, allocatable :: xC(:,:,:,:),wC(:,:,:,:),dxC(:,:,:,:)

      ! clean small values
      uniform_tol=1.d-10
      do ix1=1,nxO1
        do ix2=1,nxO2
          do j=1,nWO
            if (abs(wO(ix1,ix2,j))<smalldouble) wO(ix1,ix2,j)=zero
          enddo
        enddo
      enddo

      ! how many cells in each grid
      if (dat_resolution) then
        if (datatype=='image_euv' .or. datatype=='image_sxr') then
          if (LOS_phi==0 .and. LOS_theta==90) then
            piece_nmax1=block_nx2
            piece_nmax2=block_nx3
          else if (LOS_phi==90 .and. LOS_theta==90) then
            piece_nmax1=block_nx3
            piece_nmax2=block_nx1
          else
            piece_nmax1=block_nx1
            piece_nmax2=block_nx2
          endif
        else if (datatype=='spectrum_euv') then
          piece_nmax1=16
          if (direction_slit==1) then
            piece_nmax2=block_nx1
          else if (direction_slit==2) then
            piece_nmax2=block_nx2
          else
            piece_nmax2=block_nx3
          endif
        endif
      else
        piece_nmax1=20
        piece_nmax2=20
      endif
      LOOPN1: do j=piece_nmax1,1,-1
        if(mod(nXO1,j)==0) then
          nC1=j
          exit LOOPN1
        endif
      enddo LOOPN1
      LOOPN2: do j=piece_nmax2,1,-1
        if(mod(nXO2,j)==0) then
          nC2=j
          exit LOOPN2
        endif
      enddo LOOPN2
      ! how many grids
      nP1=nXO1/nC1
      nP2=nXO2/nC2
      nPiece=nP1*nP2
      nWC=nWO

      ! output images
      select case(convert_type)
        case('EIvtuCCmpi','ESvtuCCmpi','SIvtuCCmpi','WIvtuCCmpi')
          ! put data into grids
          allocate(xC(nPiece,nC1,nC2,2))
          allocate(dxC(nPiece,nC1,nC2,2))
          allocate(wC(nPiece,nC1,nC2,nWO))
          do ipc=1,nPiece
            do ixc1=1,nC1
              do ixc2=1,nC2
                ix1=mod(ipc-1,nP1)*nC1+ixc1
                ix2=floor(1.0*(ipc-1)/nP1)*nC2+ixc2
                xC(ipc,ixc1,ixc2,1)=xO1(ix1)
                xC(ipc,ixc1,ixc2,2)=xO2(ix2)
                dxC(ipc,ixc1,ixc2,1)=dxO1(ix1)
                dxC(ipc,ixc1,ixc2,2)=dxO2(ix2)
                do j=1,nWC
                  wC(ipc,ixc1,ixc2,j)=wO(ix1,ix2,j)
                enddo
              enddo
            enddo
          enddo
          ! write data into vtu file
          call write_image_vtuCC(qunit,xC,wC,dxC,nPiece,nC1,nC2,nWC,datatype)
          deallocate(xC,dxC,wC)
        case('EIvtiCCmpi','ESvtiCCmpi','SIvtiCCmpi','WIvtiCCmpi')
          if (dat_resolution .and. &
              (maxval(abs(dxO1(:)-dxO1(1)))>uniform_tol*max(one,abs(dxO1(1))) .or. &
               maxval(abs(dxO2(:)-dxO2(1)))>uniform_tol*max(one,abs(dxO2(1))))) then
            call mpistop("vti needs uniform dat-resolution image grids")
          else
            call write_image_vtiCC(qunit,xO1,xO2,dxO1,dxO2,wO,nXO1,nXO2,nWO,nC1,nC2)
          endif
        case default
          call mpistop("Error in synthesize emission: Unknown convert_type")
      end select

    end subroutine output_data
  }

    subroutine write_image_vtiCC(qunit,xO1,xO2,dxO1,dxO2,wO,nXO1,nXO2,nWO,nC1,nC2)
      ! write image data to vti
      use mod_global_parameters

      integer, intent(in) :: qunit,nXO1,nXO2,nWO,nC1,nC2
      double precision, intent(in) :: xO1(nXO1),xO2(nxO2)
      double precision, intent(in) :: dxO1(nxO1),dxO2(nxO2)
      double precision, intent(in) :: wO(nXO1,nXO2,nWO)

      double precision :: origin(1:3), spacing(1:3)
      integer :: wholeExtent(1:6)
      integer :: iw
      integer :: ixC1,ixC2

      integer :: filenr
      logical :: fileopen
      character (70) :: subname,wname,vname,nameL,nameS
      character (len=std_len) :: filename
      logical :: sph_datres_no_doppler


      origin(1)=xO1(1)-0.5d0*dxO1(1)
      origin(2)=xO2(1)-0.5d0*dxO2(1)
      origin(3)=zero
      spacing(1)=dxO1(1)
      spacing(2)=dxO2(1)
      spacing(3)=one
      wholeExtent=0
      wholeExtent(2)=nXO1
      wholeExtent(4)=nXO2
      sph_datres_no_doppler=dat_resolution .and. coordinate==spherical .and. trim(ray_method_active)=='spherical'

      if (mype==0) then
        inquire(qunit,opened=fileopen)
        if(.not.fileopen)then
          ! generate filename 
          filenr=snapshotini
          if (autoconvert) filenr=snapshotnext
          if (convert_type=='EIvtiCCmpi') then
            write(filename,'(a,i4.4,a)') trim(filename_euv),filenr,".vti"
          else if (convert_type=='SIvtiCCmpi') then
            write(filename,'(a,i4.4,a)') trim(filename_sxr),filenr,".vti"
          else if (convert_type=='WIvtiCCmpi') then
            write(filename,'(a,i4.4,a)') trim(filename_whitelight),filenr,".vti"
          else if (convert_type=='ESvtiCCmpi') then
            write(filename,'(a,i4.4,a)') trim(filename_spectrum),filenr,".vti"
          endif
          open(qunit,file=filename,status='unknown',form='formatted')
        endif

        ! generate xml header
        write(qunit,'(a)')'<?xml version="1.0"?>'
        write(qunit,'(a)',advance='no') '<VTKFile type="ImageData"'
        write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
        write(qunit,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')'  <ImageData Origin="',&
              origin,'" WholeExtent="',wholeExtent,'" Spacing="',spacing,'">'
        ! file info        
        write(qunit,'(a)')'<FieldData>'
        write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                           'NumberOfTuples="1" format="ascii">'
        write(qunit,*) real(global_time*time_convert_factor)
        write(qunit,'(a)')'</DataArray>'
        write(qunit,'(a)')'</FieldData>'
        ! pixel/cell data
        write(qunit,'(a,6(i10),a)') '<Piece Extent="',wholeExtent,'">'
        write(qunit,'(a)')'<CellData>'
        do iw=1,nWO
          ! variable name
          if (convert_type=='EIvtiCCmpi') then
            if (wavelength<100) then
              write(vname,'(a,i2)') "AIA",wavelength
            else if (wavelength<1000) then
              write(vname,'(a,i3)') "AIA",wavelength
            else
              write(vname,'(a,i4)') "IRIS",wavelength
            endif
            if (trim(emission_model)=='pseudo_current' .and. iw==1) vname='pseudo_current'
            if (trim(emission_model)=='radio_ff' .and. iw==1) vname='radio_brightness_temperature'
            if (trim(radiation_transfer)=='thick' .and. iw==1) vname=trim(vname)//'_thick'
            if (iw==2 .and. dat_resolution .and. (.not. sph_datres_no_doppler) .and. &
                trim(emission_model)/='radio_ff' .and. &
                trim(emission_model)/='pseudo_current') vname='Doppler_velocity'
            if (output_tau .and. trim(radiation_transfer)=='thick' .and. &
                ((trim(emission_model)=='radio_ff' .and. iw==2) .or. &
                 (trim(emission_model)/='radio_ff' .and. trim(emission_model)/='pseudo_current' .and. &
                  ((dat_resolution .and. ((sph_datres_no_doppler .and. iw==2) .or. &
                                          ((.not. sph_datres_no_doppler) .and. iw==3))) .or. &
                   ((.not. dat_resolution) .and. iw==2))))) then
              vname='tau'
            endif
            if (output_absorption_fraction .and. trim(radiation_transfer)=='thick' .and. &
                ((trim(emission_model)=='radio_ff' .and. ((output_tau .and. iw==3) .or. &
                                                         ((.not. output_tau) .and. iw==2))) .or. &
                 (trim(emission_model)/='radio_ff' .and. trim(emission_model)/='pseudo_current' .and. &
                  ((dat_resolution .and. sph_datres_no_doppler .and. output_tau .and. iw==3) .or. &
                   (dat_resolution .and. sph_datres_no_doppler .and. (.not. output_tau) .and. iw==2) .or. &
                   (dat_resolution .and. (.not. sph_datres_no_doppler) .and. output_tau .and. iw==4) .or. &
                   (dat_resolution .and. (.not. sph_datres_no_doppler) .and. (.not. output_tau) .and. iw==3) .or. &
                   ((.not. dat_resolution) .and. output_tau .and. iw==3) .or. &
                   ((.not. dat_resolution) .and. (.not. output_tau) .and. iw==2))))) then
              vname='absorption_fraction'
            endif
          else if (convert_type=='SIvtiCCmpi') then
            if (emin_sxr<10 .and. emax_sxr<10) then
              write(vname,'(a,i1,a,i1,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
            else if (emin_sxr<10 .and. emax_sxr>=10) then
              write(vname,'(a,i1,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
            else
              write(vname,'(a,i2,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
            endif
          else if (convert_type=='WIvtiCCmpi') then
            if (iw==1) write(vname,'(a)')'B'
            if (iw==2) write(vname,'(a)')'pB'
          else if (convert_type=='ESvtiCCmpi') then
            if (spectrum_wl==1354) then
              write(vname,'(a,i4)') "SG",spectrum_wl
            else
              write(vname,'(a,i3)') "EIS",spectrum_wl
            endif
          endif
          write(qunit,'(a,a,a)')&
            '<DataArray type="Float64" Name="',TRIM(vname),'" format="ascii">'
          write(qunit,'(200(1pe14.6))') ((wO(ixC1,ixC2,iw),ixC1=1,nXO1),ixC2=1,nXO2)
          write(qunit,'(a)')'</DataArray>'
        enddo
        write(qunit,'(a)')'</CellData>'
        write(qunit,'(a)')'</Piece>'
        ! end
        write(qunit,'(a)')'</ImageData>'
        write(qunit,'(a)')'</VTKFile>'
        close(qunit)
      endif

    end subroutine write_image_vtiCC

    subroutine write_image_vtuCC(qunit,xC,wC,dxC,nPiece,nC1,nC2,nWC,datatype)
      ! write image data to vtu
      use mod_global_parameters

      integer, intent(in) :: qunit
      integer, intent(in) :: nPiece,nC1,nC2,nWC
      double precision, intent(in) :: xC(nPiece,nC1,nC2,2),dxC(nPiece,nc1,nc2,2)
      double precision, intent(in) :: wC(nPiece,nC1,nC2,nWC)
      character(20), intent(in) :: datatype

      integer :: nP1,nP2
      double precision :: xP(nPiece,nC1+1,nC2+1,2)
      integer :: filenr
      logical :: fileopen
      character (70) :: subname,wname,vname,nameL,nameS
      character (len=std_len) :: filename
      integer :: ixC1,ixC2,ixP,ix1,ix2,j
      integer :: nc,np,icel,VTK_type
      logical :: sph_datres_no_doppler

      nP1=nC1+1
      nP2=nC2+1
      np=nP1*nP2
      nc=nC1*nC2
      sph_datres_no_doppler=dat_resolution .and. coordinate==spherical .and. trim(ray_method_active)=='spherical'
      ! cell corner location     
      do ixP=1,nPiece
        do ix1=1,nP1
          do ix2=1,nP2
            if (ix1<nP1) xP(ixP,ix1,ix2,1)=xC(ixP,ix1,1,1)-0.5d0*dxC(ixP,ix1,1,1)
            if (ix1==nP1) xP(ixP,ix1,ix2,1)=xC(ixP,ix1-1,1,1)+0.5d0*dxC(ixP,ix1-1,1,1)
            if (ix2<nP2) xP(ixP,ix1,ix2,2)=xC(ixP,1,ix2,2)-0.5d0*dxC(ixP,1,ix2,2)
            if (ix2==nP2) xP(ixP,ix1,ix2,2)=xC(ixP,1,ix2-1,2)+0.5d0*dxC(ixP,1,ix2-1,2)
          enddo
        enddo
      enddo
      if (mype==0) then
        inquire(qunit,opened=fileopen)
        if(.not.fileopen)then
          ! generate filename 
          filenr=snapshotini
          if (autoconvert) filenr=snapshotnext
          if (datatype=='image_euv') then
            write(filename,'(a,i4.4,a)') trim(filename_euv),filenr,".vtu"
          else if (datatype=='image_sxr') then
            write(filename,'(a,i4.4,a)') trim(filename_sxr),filenr,".vtu"
          else if (datatype=='image_whitelight') then
            write(filename,'(a,i4.4,a)') trim(filename_whitelight),filenr,".vtu"
          else if (datatype=='spectrum_euv') then
            write(filename,'(a,i4.4,a)') trim(filename_spectrum),filenr,".vtu"
          endif
          open(qunit,file=filename,status='unknown',form='formatted')
        endif
        ! generate xml header
        write(qunit,'(a)')'<?xml version="1.0"?>'
        write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
        write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
        write(qunit,'(a)')'<UnstructuredGrid>'
        write(qunit,'(a)')'<FieldData>'
        write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                           'NumberOfTuples="1" format="ascii">'
        write(qunit,*) real(global_time*time_convert_factor)
        write(qunit,'(a)')'</DataArray>'
        write(qunit,'(a)')'</FieldData>'
        do ixP=1,nPiece
          write(qunit,'(a,i7,a,i7,a)') &
                '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do j=1,nWC
            if (datatype=='image_euv') then
              if (j==1) then
                if (wavelength<100) then
                  write(vname,'(a,i2)') "AIA",wavelength
                else if (wavelength<1000) then
                  write(vname,'(a,i3)') "AIA",wavelength
                else
                  write(vname,'(a,i4)') "IRIS",wavelength
                endif
                if (trim(emission_model)=='pseudo_current') vname='pseudo_current'
                if (trim(emission_model)=='radio_ff') vname='radio_brightness_temperature'
                if (trim(radiation_transfer)=='thick') vname=trim(vname)//'_thick'
              endif
              if (j==2 .and. dat_resolution .and. (.not. sph_datres_no_doppler) .and. &
                  trim(emission_model)/='radio_ff' .and. &
                  trim(emission_model)/='pseudo_current') vname='Doppler_velocity'
              if (output_tau .and. trim(radiation_transfer)=='thick' .and. &
                  ((trim(emission_model)=='radio_ff' .and. j==2) .or. &
                   (trim(emission_model)/='radio_ff' .and. trim(emission_model)/='pseudo_current' .and. &
                    ((dat_resolution .and. ((sph_datres_no_doppler .and. j==2) .or. &
                                            ((.not. sph_datres_no_doppler) .and. j==3))) .or. &
                     ((.not. dat_resolution) .and. j==2))))) then
                vname='tau'
              endif
              if (output_absorption_fraction .and. trim(radiation_transfer)=='thick' .and. &
                  ((trim(emission_model)=='radio_ff' .and. ((output_tau .and. j==3) .or. &
                                                           ((.not. output_tau) .and. j==2))) .or. &
                   (trim(emission_model)/='radio_ff' .and. trim(emission_model)/='pseudo_current' .and. &
                    ((dat_resolution .and. sph_datres_no_doppler .and. output_tau .and. j==3) .or. &
                     (dat_resolution .and. sph_datres_no_doppler .and. (.not. output_tau) .and. j==2) .or. &
                     (dat_resolution .and. (.not. sph_datres_no_doppler) .and. output_tau .and. j==4) .or. &
                     (dat_resolution .and. (.not. sph_datres_no_doppler) .and. (.not. output_tau) .and. j==3) .or. &
                     ((.not. dat_resolution) .and. output_tau .and. j==3) .or. &
                     ((.not. dat_resolution) .and. (.not. output_tau) .and. j==2))))) then
                vname='absorption_fraction'
              endif
            else if (datatype=='image_sxr') then
              if (emin_sxr<10 .and. emax_sxr<10) then
                write(vname,'(a,i1,a,i1,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
              else if (emin_sxr<10 .and. emax_sxr>=10) then
                write(vname,'(a,i1,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
              else
                write(vname,'(a,i2,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
              endif
            else if (datatype=='image_whitelight') then
              write(vname,'(a)')'whitelight'
            else if (datatype=='spectrum_euv') then
              if (spectrum_wl==1354) then
                write(vname,'(a,i4)') "SG",spectrum_wl
              else
                write(vname,'(a,i3)') "EIS",spectrum_wl
              endif
            endif
            write(qunit,'(a,a,a)')&
              '<DataArray type="Float64" Name="',TRIM(vname),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') ((wC(ixP,ixC1,ixC2,j),ixC1=1,nC1),ixC2=1,nC2)
            write(qunit,'(a)')'</DataArray>'
          enddo
          write(qunit,'(a)')'</CellData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)')'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          do ix2=1,nP2
            do ix1=1,nP1
              if (datatype=='image_euv' .and. dat_resolution) then
                if (LOS_phi==0 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') 0.d0,xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2)
                else if (LOS_phi==90 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,2),0.d0,xP(ixP,ix1,ix2,1)
                else
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2),0.d0
                endif
              else if (datatype=='image_sxr' .and. dat_resolution) then
                if (LOS_phi==0 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') 0.d0,xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2)
                else if (LOS_phi==90 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,2),0.d0,xP(ixP,ix1,ix2,1)
                else
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2),0.d0
                endif
              else
                write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2),0.d0
              endif
            enddo
          enddo
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
          ! connetivity part
          write(qunit,'(a)')'<Cells>'
          write(qunit,'(a)')'<DataArray type="Int32" Name="connectivity" format="ascii">'
          do ix2=1,nC2
            do ix1=1,nC1
              write(qunit,'(4(i7))') ix1-1+(ix2-1)*nP1,ix1+(ix2-1)*nP1,&
                                     ix1-1+ix2*nP1,ix1+ix2*nP1
            enddo
          enddo
          write(qunit,'(a)')'</DataArray>'
          ! offsets data array
          write(qunit,'(a)')'<DataArray type="Int32" Name="offsets" format="ascii">'
          do icel=1,nc
            write(qunit,'(i7)') icel*(2**2)
          enddo
          write(qunit,'(a)')'</DataArray>'
          ! VTK cell type data array
          write(qunit,'(a)')'<DataArray type="Int32" Name="types" format="ascii">'
          ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
          VTK_type=8
          do icel=1,nc
            write(qunit,'(i2)') VTK_type
          enddo
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Cells>'
          write(qunit,'(a)')'</Piece>'
        enddo
        write(qunit,'(a)')'</UnstructuredGrid>'
        write(qunit,'(a)')'</VTKFile>'
        close(qunit)
      endif
    end subroutine write_image_vtuCC

    subroutine dot_product_loc(vec1,vec2,res)
      double precision, intent(in) :: vec1(1:3),vec2(1:3)
      double precision, intent(out) :: res

      res=vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)

    end subroutine dot_product_loc

    subroutine cross_product_loc(vec_in1,vec_in2,vec_out)
      double precision, intent(in) :: vec_in1(1:3),vec_in2(1:3)
      double precision, intent(out) :: vec_out(1:3)

      vec_out(1)=vec_in1(2)*vec_in2(3)-vec_in1(3)*vec_in2(2)
      vec_out(2)=vec_in1(3)*vec_in2(1)-vec_in1(1)*vec_in2(3)
      vec_out(3)=vec_in1(1)*vec_in2(2)-vec_in1(2)*vec_in2(1)

    end subroutine cross_product_loc

    subroutine init_vectors_spherical()
      integer :: j
      double precision :: LOS_psi
      double precision :: vec_car(1:3),vec_z(1:3),vec_temp1(1:3),vec_temp2(1:3)
      double precision :: vec_LOS_sph(1:3),vec_xI1_sph(1:3),vec_xI2_sph(1:3)

      ! antiparallel to LOS in spherical
      vec_LOS(1)=1.d0
      vec_LOS(2)=dpi*LOS_theta/180.d0
      vec_LOS(3)=dpi*LOS_phi/180.d0
      ! LOS in cartesian
      call spherical_to_cartesian(vec_LOS,vec_car)
      vec_LOS=-vec_car

      ! theta=0 in cartesian
      vec_z(:)=zero
      vec_z(3)=1.d0

      ! x direction for image
      if (LOS_theta==zero) then
        vec_temp1(1)=1.d0
        vec_temp1(2)=dpi/2.d0
        vec_temp1(3)=dpi*LOS_phi/180.d0
        call spherical_to_cartesian(vec_temp1,vec_car)
        vec_temp1=-vec_car
        call cross_product_loc(vec_temp1,vec_z,vec_xI1)
      else
        call cross_product_loc(vec_LOS,vec_z,vec_xI1)
      endif

      ! y direction for image
      call cross_product_loc(vec_xI1,vec_LOS,vec_xI2)

      ! rotate the image
      vec_temp1=vec_xI1/sqrt(vec_xI1(1)**2+vec_xI1(2)**2+vec_xI1(3)**2)
      vec_temp2=vec_xI2/sqrt(vec_xI2(1)**2+vec_xI2(2)**2+vec_xI2(3)**2)
      LOS_psi=dpi*image_rotate/180.d0
      vec_xI1=vec_temp1*cos(LOS_psi)-vec_temp2*sin(LOS_psi)
      vec_xI2=vec_temp2*cos(LOS_psi)+vec_temp1*sin(LOS_psi)

      do j=1,3
        if (abs(vec_LOS(j))<smalldouble) vec_LOS(j)=zero
        if (abs(vec_xI1(j))<smalldouble) vec_xI1(j)=zero
        if (abs(vec_xI2(j))<smalldouble) vec_xI2(j)=zero
      enddo

      call cartesian_to_spherical(vec_LOS,vec_LOS_sph)
      call cartesian_to_spherical(vec_xI1,vec_xI1_sph)
      call cartesian_to_spherical(vec_xI2,vec_xI2_sph)
      vec_LOS_sph(2:3)=vec_LOS_sph(2:3)*180.d0/dpi
      vec_xI1_sph(2:3)=vec_xI1_sph(2:3)*180.d0/dpi
      vec_xI2_sph(2:3)=vec_xI2_sph(2:3)*180.d0/dpi

      if (mype==0) write(*,'(a,f3.1,f6.1,f6.1,a)') ' ray direction (spherical): [',vec_LOS_sph(1),vec_LOS_sph(2),vec_LOS_sph(3),']'
      if (mype==0) write(*,'(a,f3.1,f6.1,f6.1,a)') ' xI1 direction (spherical): [',vec_xI1_sph(1),vec_xI1_sph(2),vec_xI1_sph(3),']'
      if (mype==0) write(*,'(a,f3.1,f6.1,f6.1,a)') ' xI2 direction (spherical): [',vec_xI2_sph(1),vec_xI2_sph(2),vec_xI2_sph(3),']'

    end subroutine init_vectors_spherical

    subroutine spherical_to_cartesian(vec_sph,vec_car)
      ! angles in rad
      double precision, intent(in) :: vec_sph(1:3)
      double precision, intent(inout) :: vec_car(1:3)

      vec_car(1)=vec_sph(1)*dsin(vec_sph(2))*dcos(vec_sph(3))
      vec_car(2)=vec_sph(1)*dsin(vec_sph(2))*dsin(vec_sph(3))
      vec_car(3)=vec_sph(1)*dcos(vec_sph(2))

    end subroutine spherical_to_cartesian

    subroutine cartesian_to_spherical(vec_car,vec_sph)
      ! angles in rad
      double precision, intent(in) :: vec_car(1:3)
      double precision, intent(inout) :: vec_sph(1:3)

      vec_sph(1)=dsqrt(vec_car(1)**2+vec_car(2)**2+vec_car(3)**2)
      vec_sph(2)=dacos(vec_car(3)/vec_sph(1))
      vec_sph(3)=atan2(vec_car(2),vec_car(1))

    end subroutine cartesian_to_spherical

    subroutine init_vectors_cartesian()
      integer :: j
      double precision :: LOS_psi
      double precision :: vec_z(1:3),vec_temp1(1:3),vec_temp2(1:3)

      ! vectors for image coordinate
      vec_LOS(1)=-cos(dpi*LOS_phi/180.d0)*sin(dpi*LOS_theta/180.d0)
      vec_LOS(2)=-sin(dpi*LOS_phi/180.d0)*sin(dpi*LOS_theta/180.d0)
      vec_LOS(3)=-cos(dpi*LOS_theta/180.d0)
      do j=1,3
        if (abs(vec_LOS(j))<=smalldouble) vec_LOS(j)=zero
      enddo
      vec_z(:)=zero
      vec_z(3)=1.d0
      if (LOS_theta==zero) then
        vec_xI1(1)=cos(dpi*LOS_phi/180.d0)
        vec_xI1(2)=sin(dpi*LOS_phi/180.d0)
        vec_xI1(3)=zero
      else
        call cross_product_loc(vec_LOS,vec_z,vec_xI1)
      endif
      call cross_product_loc(vec_xI1,vec_LOS,vec_xI2)
      vec_temp1=vec_xI1/sqrt(vec_xI1(1)**2+vec_xI1(2)**2+vec_xI1(3)**2)
      vec_temp2=vec_xI2/sqrt(vec_xI2(1)**2+vec_xI2(2)**2+vec_xI2(3)**2)
      LOS_psi=dpi*image_rotate/180.d0
      vec_xI1=vec_temp1*cos(LOS_psi)-vec_temp2*sin(LOS_psi)
      vec_xI2=vec_temp2*cos(LOS_psi)+vec_temp1*sin(LOS_psi)

      do j=1,3
        if (abs(vec_xI1(j))<smalldouble) vec_xI1(j)=zero
        if (abs(vec_xI2(j))<smalldouble) vec_xI2(j)=zero
      enddo

      if (mype==0) write(*,'(a,f5.2,f6.2,f6.2,a)') ' LOS vector: [',vec_LOS(1),vec_LOS(2),vec_LOS(3),']'
      if (mype==0) write(*,'(a,f5.2,f6.2,f6.2,a)') ' xI1 vector: [',vec_xI1(1),vec_xI1(2),vec_xI1(3),']'
      if (mype==0) write(*,'(a,f5.2,f6.2,f6.2,a)') ' xI2 vector: [',vec_xI2(1),vec_xI2(2),vec_xI2(3),']'

    end subroutine init_vectors_cartesian

    subroutine get_cor_image_spherical(x_3D_sph,x_image)
      double precision, intent(in) :: x_3D_sph(1:3)
      double precision, intent(inout) :: x_image(1:2)
      double precision :: res,res_origin
      double precision :: x_3D(1:3)

      call spherical_to_cartesian(x_3D_sph,x_3D)
      call dot_product_loc(x_3D,vec_xI1,res)
      x_image(1)=res
      call dot_product_loc(x_3D,vec_xI2,res)
      x_image(2)=res

    end subroutine get_cor_image_spherical

    subroutine get_cor_image(x_3D,x_image)
      double precision, intent(in) :: x_3D(1:3)
      double precision, intent(inout) :: x_image(1:2)
      double precision :: res,res_origin

      call dot_product_loc(x_3D,vec_xI1,res)
      call dot_product_loc(x_origin,vec_xI1,res_origin)
      x_image(1)=res-res_origin
      call dot_product_loc(x_3D,vec_xI2,res)
      call dot_product_loc(x_origin,vec_xI2,res_origin)
      x_image(2)=res-res_origin

    end subroutine get_cor_image

end module mod_thermal_emission

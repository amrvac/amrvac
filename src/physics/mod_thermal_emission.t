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
  use mod_physics

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

  type te_fluid

    procedure (get_subr1), pointer, nopass :: get_rho => null()
    procedure (get_subr1), pointer, nopass :: get_pthermal => null()

    ! factor in eq of state p = Rfactor * rho * T
    ! used for getting temperature
    double precision :: Rfactor = 1d0

  end type te_fluid


  contains

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
    end subroutine get_line_info
    
    subroutine get_EUV(wl,ixI^L,ixO^L,w,x,fl,flux)
      ! calculate the local emission intensity of given EUV line (optically thin)
      ! wavelength is the wave length of the emission line
      ! unit [DN cm^-1 s^-1 pixel^-1]
      ! ingrate flux along line of sight: DN s^-1 pixel^-1
      use mod_global_parameters

      integer, intent(in) :: wl
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: flux(ixI^S)

      integer :: n_table
      double precision, allocatable :: t_table(:),f_table(:)
      integer :: ix^D,iTt,i
      double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)
      double precision :: logT,logGT,GT

      ! selecting emission table 
      select case(wl)
      case(94)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_94(1:n_aia)
      case(131)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_131(1:n_aia)
      case(171)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_171(1:n_aia)
      case(193)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_193(1:n_aia)
      case(211)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_211(1:n_aia)
      case(304)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_304(1:n_aia)
      case(335)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_335(1:n_aia)
      case(1354)
        n_table=n_iris
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_iris(1:n_iris)
        f_table(1:n_table)=f_1354(1:n_iris)
      case(263)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis1(1:n_eis)
        f_table(1:n_table)=f_263(1:n_eis)
      case(264)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis2(1:n_eis)
        f_table(1:n_table)=f_264(1:n_eis)
      case(192)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis2(1:n_eis)
        f_table(1:n_table)=f_192(1:n_eis)
      case(255)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis2(1:n_eis)
        f_table(1:n_table)=f_255(1:n_eis)
      case default
        allocate(t_table(1))
        allocate(f_table(1))
        call mpistop("Unknown wavelength")
      end select
      call fl%get_pthermal(w,x,ixI^L,ixO^L,pth)
      call fl%get_rho(w,x,ixI^L,ixO^L,Ne)
      Te(ixO^S)=pth(ixO^S)/(Ne(ixO^S)*fl%Rfactor)*unit_temperature
      if (SI_unit) then
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity/1.d6 ! m^-3 -> cm-3
        flux(ixO^S)=Ne(ixO^S)**2
      else
        Ne(ixO^S)=Ne(ixO^S)*unit_numberdensity
        flux(ixO^S)=Ne(ixO^S)**2
      endif

      select case(wl)
      case(94,131,171,193,211,304,335,1354)
      ! temperature table in log
        do i=1,n_table
          if(f_table(i) .gt. 1.d-99) then
            f_table(i)=log10(f_table(i))
          else
            f_table(i)=-99.d0
          endif
        enddo 
        logGT=zero
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
          logT=log10(Te(ix^D))
          if (logT>=t_table(1) .and. logT<=t_table(n_table)) then
            do iTt=1,n_table-1
              if (logT>=t_table(iTt) .and. logT<t_table(iTt+1)) then
                logGT=f_table(iTt)*(logT-t_table(iTt+1))/(t_table(iTt)-t_table(iTt+1))+&
                      f_table(iTt+1)*(logT-t_table(iTt))/(t_table(iTt+1)-t_table(iTt))
              endif
            enddo
            flux(ix^D)=flux(ix^D)*(10**(logGT))
            if(flux(ix^D)<smalldouble) flux(ix^D)=0.d0
          else
            flux(ix^D)=zero
          endif
        {enddo\}
      case default
      ! temperature table linear
        {do ix^DB=ixOmin^DB,ixOmax^DB\}
          if (Te(ix^D)>=t_table(1) .and. Te(ix^D)<=t_table(n_table)) then
            do iTt=1,n_table-1
              if (Te(ix^D)>=t_table(iTt) .and. Te(ix^D)<t_table(iTt+1)) then
                GT=f_table(iTt)*(Te(ix^D)-t_table(iTt+1))/(t_table(iTt)-t_table(iTt+1))+&
                   f_table(iTt+1)*(Te(ix^D)-t_table(iTt))/(t_table(iTt+1)-t_table(iTt))
              endif
            enddo
            flux(ix^D)=flux(ix^D)*GT
            if(flux(ix^D)<smalldouble) flux(ix^D)=0.d0
          else
            flux(ix^D)=zero
          endif
        {enddo\}
      end select

      deallocate(t_table,f_table)
    end subroutine get_EUV

    subroutine get_SXR(ixI^L,ixO^L,w,x,fl,flux,El,Eu)
      !synthesize thermal SXR from El keV to Eu keV
      !flux (cgs): photons cm^-5 s^-1
      !flux (SI): photons m^-3 cm^-2 s^-1
      !integration of the flux is the SXR flux observed at 1AU [photons cm^-2 s^-1]
      use mod_global_parameters

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

      I0=1.07d-42    ! photon flux index for observed at 1AU [photon cm s^-1 keV^-1]
      kb=const_kb
      keV=1.0d3*const_ev
      dE=0.1
      numE=floor((Eu-El)/dE)
      call fl%get_pthermal(w,x,ixI^L,ixO^L,pth)
      call fl%get_rho(w,x,ixI^L,ixO^L,Ne)

      Te(ixO^S)=pth(ixO^S)/(Ne(ixO^S)*fl%Rfactor)*unit_temperature
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
          if(kbT(ix^D)<Ei) then
            gff(ix^D)=(kbT(ix^D)/Ei)**0.4
          endif
        {enddo\}
        fi(ixO^S)=(EM(ixO^S)*gff(ixO^S))* &
                  exp(-Ei/(kbT(ixO^S)))/(Ei*sqrt(kbT(ixO^S)))
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
        Te(ixb^S)=pth(ixb^S)/(Ne(ixb^S)*fl%Rfactor)*unit_temperature
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
            if(kbT(ix^D)<Ei) then
              gff=(kbT(ix^D)/Ei)**0.4
            else
              gff=1.d0
            endif
            fi=(EM(ix^D)*gff)*exp(-Ei/(kbT(ix^D)))/(Ei*sqrt(kbT(ix^D)))
            eflux_grid=eflux_grid+fi*dE*Ei
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

      arcsec=7.25d7/unit_length

      if (mype==0) print *, '###################################################'
      select case(spectrum_wl)
      case (1354)
        if (mype==0) print *, 'Systhesizing EUV spectrum (observed by IRIS).'
      case (263,264,192,255)
        if (mype==0) print *, 'Systhesizing EUV spectrum (observed by Hinode/EIS).'
      case default
        call MPISTOP('Wrong wavelength!')
      end select

      call get_line_info(spectrum_wl,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
      if (mype==0) write(*,'(a,f8.3,a)') ' Wavelength: ',lineCent,' Angstrom'
      if (mype==0) print *, 'Unit of EUV flux: DN s^-1 pixel^-1'
      if (mype==0) write(*,'(a,f5.3,a,f5.1,a)') ' Pixel: ',wlRsl,' Angstrom x ',spaceRsl*725.0, ' km'

      if (spectrum_window_max<=spectrum_window_min) then
        call MPISTOP('Wrong spectrum window!')
      endif

      datatype='spectrum_euv'

      if (resolution_spectrum=='data') then
        if (mype==0) print *, 'Unit of wavelength: Angstrom (0.1 nm) '
        if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
        if (SI_unit) then
          if (mype==0) write(*,'(a,f8.1,a)') ' Location of slit: ',location_slit*unit_length/1.d6,' Mm'
        else
          if (mype==0) write(*,'(a,f8.1,a)') ' Location of slit: ',location_slit*unit_length/1.d8,' Mm'
        endif
        if (mype==0) write(*,'(a,f8.1,a)') ' Width of slit: ',wslit*725.0,' km'
        call get_spectrum_data_resol(qunit,datatype,fl)
      else if (resolution_spectrum=='instrument') then
        if (mype==0) print *, 'Unit of wavelength: Angstrom (0.1 nm) '
        if (mype==0) print *, 'Unit of length: arcsec (~725 km)'
        if (mype==0) print *, 'Direction of the slit: parallel to xI2 vector'
        if (mype==0) write(*,'(a,f8.1,a)') ' Location of slit: xI1 = ',location_slit,' arcsec'
        if (mype==0) write(*,'(a,f8.1,a)') ' Width of slit: ',wslit,' arcsec'
        call get_spectrum_inst_resol(qunit,datatype,fl)
      else
        call MPISTOP('Wrong resolution for resolution_spectrum!')
      endif

      if (mype==0) print *, '###################################################'

    end subroutine get_EUV_spectrum

    subroutine get_spectrum_inst_resol(qunit,datatype,fl)

      integer, intent(in) :: qunit
      character(20), intent(in) :: datatype
      type(te_fluid), intent(in) :: fl

      integer :: numWL,numXS,iwL,ixS,numWI,ix^D
      double precision :: dwLg,dxSg,xSmin,xSmax,xScent,wLmin,wLmax
      double precision, allocatable :: wL(:),xS(:),dwL(:),dxS(:)
      double precision, allocatable :: wI(:,:,:),spectra(:,:),spectra_rc(:,:)
      double precision :: vec_cor(1:3),xI_cor(1:2)
      double precision :: res

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl,wslit
      double precision :: unitv,arcsec,RHESSI_rsl,pixel
      integer :: iigrid,igrid,i,j,numS
      double precision :: xLmin,xLmax,xslit

      call init_vectors()

      ! calculate domain in space
      do ix1=1,2
        if (ix1==1) vec_cor(1)=xprobmin1
        if (ix1==2) vec_cor(1)=xprobmax1
        do ix2=1,2
          if (ix2==1) vec_cor(2)=xprobmin2
          if (ix2==2) vec_cor(2)=xprobmax2
          do ix3=1,2
            if (ix3==1) vec_cor(3)=xprobmin3
            if (ix3==2) vec_cor(3)=xprobmax3
            call dot_product_loc(vec_cor,vec_xI2,res)
            xI_cor(2)=res
            if (ix1==1 .and. ix2==1 .and. ix3==1) then
              xSmin=xI_cor(2)
              xSmax=xI_cor(2)
            else
              xSmin=min(xSmin,xI_cor(2))
              xSmax=max(xSmax,xI_cor(2))
            endif
          enddo
        enddo
      enddo
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
              call dot_product_loc(vec_cor,vec_xI1,res)
              xI_cor(1)=res
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


        xslit=location_slit*arcsec
        !if (location_slit>=xLmin-dxSg .and. location_slit<=xLmax+dxSg) then
        if (xslit>=xLmin-wslit*arcsec .and. xslit<=xLmax+wslit*arcsec) then
          call integrate_spectra_inst_resol(igrid,wL,dwLg,xS,dxSg,spectra,numWL,numXS,fl)
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

      xS=xS/arcsec
      dxS=dxS/arcsec

      call output_data(qunit,wL,xS,dwL,dxS,wI,numWL,numXS,numWI,datatype)

      deallocate(wL,xS,dwL,dxS,spectra,spectra_rc,wI)

    end subroutine get_spectrum_inst_resol

    subroutine integrate_spectra_inst_resol(igrid,wL,dwLg,xS,dxSg,spectra,numWL,numXS,fl)

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
      double precision :: slit_width,dxSubC^D,xCent1,xCent2,xerf^L,fluxSubC
      double precision :: xSubC(1:3)

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
      xslit=location_slit*arcsec

      call get_line_info(spectrum_wl,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      allocate(flux(ixI^S),v(ixI^S),pth(ixI^S),Te(ixI^S),rho(ixI^S))
      ! get local EUV flux and velocity
      call get_EUV(spectrum_wl,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
      call fl%get_pthermal(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,pth)
      call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
      Te(ixO^S)=pth(ixO^S)/(fl%Rfactor*rho(ixO^S))
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
        xloc(1:3)=ps(igrid)%x(ix^D,1:3)
        dxloc(1:3)=ps(igrid)%dx(ix^D,1:3)
        call dot_product_loc(xloc,vec_xI1,res)
        xIloc(1)=res
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
            call dot_product_loc(xSubC,vec_xI1,xCent1)
            dst_slit=abs(xCent1-xslit)  ! space distance to slit center
            if (dst_slit<=half*slit_width) then
              call dot_product_loc(xSubC,vec_xI2,xCent2)  ! get sub cell center
              ixS=floor((xCent2-(xS(1)-half*dxSg))/dxSg)+1
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
                  xerfmin2=(xS(ixS)-half*dxSg-xCent2)/(sqrt(2.d0)*sigma_xs)
                  xerfmax2=(xS(ixS)+half*dxSg-xCent2)/(sqrt(2.d0)*sigma_xs)
                  factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
                  spectra(iwL,ixS)=spectra(iwL,ixS)+fluxSubC*factor
                enddo
              enddo
              ! nearby pixels
            endif
          {enddo\}
        endif
      {enddo\}

      deallocate(flux,v,pth,Te)
    end subroutine integrate_spectra_inst_resol

    subroutine get_spectrum_data_resol(qunit,datatype,fl)

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
      case (2)
        numXS=domain_nx2*2**(refine_max_level-1)
        xSmin=xprobmin2
        xSmax=xprobmax2
        bnx=block_nx2
        nbb=domain_nx2
        strtype=stretch_type(2)
        nstrb=nstretchedblocks_baselevel(2)
        qs=qstretch_baselevel(2)
      case (3)
        numXS=domain_nx3*2**(refine_max_level-1)
        xSmin=xprobmin3
        xSmax=xprobmax3
        bnx=block_nx3
        nbb=domain_nx3
        strtype=stretch_type(3)
        nstrb=nstretchedblocks_baselevel(3)
        qs=qstretch_baselevel(3)
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
      else if (dir_loc==2) then
        if (location_slit>xprobmax2 .or. location_slit<xprobmin2) then
          call MPISTOP('Wrong value for location_slit!')
        endif
      else
        if (location_slit>xprobmax3 .or. location_slit<xprobmin3) then
          call MPISTOP('Wrong value for location_slit!')
        endif
      endif

      ! find slit and do integration
      spectra=zero
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        ^D&xbmin(^D)=rnode(rpxmin^D_,igrid);
        ^D&xbmax(^D)=rnode(rpxmax^D_,igrid);
        if (location_slit>=xbmin(dir_loc) .and. location_slit<xbmax(dir_loc)) then
          call integrate_spectra_data_resol(igrid,wL,dwL,spectra,numWL,numXS,dir_loc,fl)
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

    end subroutine get_spectrum_data_resol

    subroutine integrate_spectra_data_resol(igrid,wL,dwL,spectra,numWL,numXS,dir_loc,fl)
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
      call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
      v(ixO^S)=-ps(igrid)%w(ixO^S,iw_mom(direction_LOS))/rho(ixO^S)
      call fl%get_pthermal(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,pth)
      Te(ixO^S)=pth(ixO^S)/(fl%Rfactor*rho(ixO^S))

      ! grid parameters
      levelg=ps(igrid)%level
      rft=2**(refine_max_level-levelg)

      {do ix^D=ixOmin^D,ixOmax^D\}
        if (SI_unit) then
          wlc=lineCent*(1.d0+v(ix^D)*unit_velocity*1.d2/const_c)
        else
          wlc=lineCent*(1.d0+v(ix^D)*unit_velocity/const_c)
        endif
        wlwd=sqrt(kb_cgs*Te(ix^D)*unit_temperature/(mass*mp_cgs))
        wlwd=wlwd*lineCent/const_c

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

        select case(direction_LOS)
        case(1)
          dL=ps(igrid)%dx(ix^D,1)*unit_length
        case(2)
          dL=ps(igrid)%dx(ix^D,2)*unit_length
        case default
          dL=ps(igrid)%dx(ix^D,3)*unit_length
        end select
        if (SI_unit) dL=dL*1.d2

        do iwL=1,numWL
          flux_pix=flux(ix^D)*wlRsl*dL*exp(-(wL(iwL)-wlc)**2/(2*wlwd**2))/(sqrt(2*dpi)*wlwd)
          flux_pix=flux_pix*wslit/spaceRsl
          spectra(iwL,ixSmin:ixSmax)=spectra(iwL,ixSmin:ixSmax)+flux_pix
        enddo

      {enddo\}

      deallocate(flux,v,pth,Te,rho)

    end subroutine integrate_spectra_data_resol

    subroutine get_EUV_image(qunit,fl)
      use mod_global_parameters

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20) :: datatype

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl,wslit

      call get_line_info(wavelength,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
      if (mype==0) print *, '###################################################'
      if (mype==0) print *, 'Systhesizing EUV image'
      if (mype==0) write(*,'(a,f8.3,a)') ' Wavelength: ',lineCent,' Angstrom'
      if (mype==0) print *, 'Unit of EUV flux: DN s^-1 pixel^-1'
      if (mype==0) write(*,'(a,f5.1,a,f5.1,a)') ' Pixel: ',spaceRsl*725.0,' km x ',spaceRsl*725.0, ' km'

      datatype='image_euv'

      if (resolution_euv=='data') then
        if (SI_unit) then
          if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d6,' Mm'
        else
          if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d8,' Mm'
        endif
        if (LOS_phi==0 .and. LOS_theta==90) then
          call get_image_data_resol(qunit,datatype,fl)
        else if (LOS_phi==90 .and. LOS_theta==90) then
          call get_image_data_resol(qunit,datatype,fl)
        else if (LOS_theta==0) then
          call get_image_data_resol(qunit,datatype,fl)
        else
          call MPISTOP('ERROR: Wrong LOS for synthesizing emission!')
        endif
      else if (resolution_euv=='instrument') then
        if (mype==0) print *, 'Unit of length: arcsec (~725 km)'
        call get_image_inst_resol(qunit,datatype,fl)
      else
        call MPISTOP('ERROR: Wrong resolution type')
      endif

      if (mype==0) print *, '###################################################'
      
    end subroutine get_EUV_image

    subroutine get_SXR_image(qunit,fl)
      use mod_global_parameters

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20) :: datatype

      if (mype==0) print *, '###################################################'
      if (mype==0) print *, 'Systhesizing SXR image (observed at 1 AU).'
      if (mype==0) write(*,'(a,i2,a,i2,a)') ' Passband: ',emin_sxr,' - ',emax_sxr,' keV'
      if (mype==0) print *, 'Unit of SXR flux: photons cm^-2 s^-1 pixel^-1'
      if (mype==0) write(*,'(a,f6.1,a,f6.1,a)') ' Pixel: ',2.3*725.0,' km x ',2.3*725.0, ' km'

      datatype='image_sxr'

      if (resolution_sxr=='data') then
        if (SI_unit) then
          if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d3,' km'
        else
          if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d5,' km'
        endif
        if (LOS_phi==0 .and. LOS_theta==90) then
          call get_image_data_resol(qunit,datatype,fl)
        else if (LOS_phi==90 .and. LOS_theta==90) then
          call get_image_data_resol(qunit,datatype,fl)
        else if (LOS_theta==0) then
          call get_image_data_resol(qunit,datatype,fl)
        else
          call MPISTOP('ERROR: Wrong LOS for synthesizing emission!')
        endif
      else if (resolution_sxr=='instrument') then
        if (mype==0) print *, 'Unit of length: arcsec (~725 km)'
        call get_image_inst_resol(qunit,datatype,fl)
      else
        call MPISTOP('ERROR: Wrong resolution type')
      endif

      if (mype==0) print *, '###################################################'

    end subroutine get_SXR_image

    subroutine get_image_inst_resol(qunit,datatype,fl)
      ! integrate emission flux along line of sight (LOS) 
      ! in a 3D simulation box and get a 2D EUV image
      use mod_global_parameters
      use mod_constants

      integer, intent(in) :: qunit
      type(te_fluid), intent(in) :: fl
      character(20), intent(in) :: datatype

      integer :: ix^D,numXI1,numXI2,numWI
      double precision :: xImin1,xImax1,xImin2,xImax2,xIcent1,xIcent2,dxI
      double precision, allocatable :: xI1(:),xI2(:),dxI1(:),dxI2(:),wI(:,:,:)
      double precision, allocatable :: EUVs(:,:),EUV(:,:),Dpls(:,:),Dpl(:,:)
      double precision, allocatable :: SXRs(:,:),SXR(:,:)
      double precision :: vec_temp1(1:3),vec_temp2(1:3)
      double precision :: vec_z(1:3),vec_cor(1:3),xI_cor(1:2)
      double precision :: res,LOS_psi

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl,wslit
      double precision :: unitv,arcsec,RHESSI_rsl,pixel
      integer :: iigrid,igrid,i,j,numSI

      call init_vectors()

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
            call dot_product_loc(vec_cor,vec_xI1,res)
            xI_cor(1)=res
            call dot_product_loc(vec_cor,vec_xI2,res)
            xI_cor(2)=res
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
          enddo
        enddo
      enddo
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
      else
        RHESSI_rsl=2.3d0
        dxI=RHESSI_rsl*arcsec
      endif
      numXI1=2*ceiling((xImax1-xIcent1)/dxI/2.d0)
      xImin1=xIcent1-numXI1*dxI
      xImax1=xIcent1+numXI1*dxI
      numXI1=numXI1*2
      numXI2=2*ceiling((xImax2-xIcent2)/dxI/2.d0)
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
      if (datatype=='image_euv') then
        if (SI_unit) then
          unitv=unit_velocity/1.0e3 ! km/s
        else
          unitv=unit_velocity/1.0e5 ! km/s
        endif
        numWI=2
        allocate(wI(numXI1,numXI2,numWI))
        allocate(EUV(numXI1,numXI2),EUVs(numXI1,numXI2))
        allocate(Dpl(numXI1,numXI2),Dpls(numXI1,numXI2))
        wI=zero
        EUV=zero
        EUVs=zero
        Dpl=zero
        Dpls=zero
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call integrate_EUV_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,EUVs,Dpls)
        enddo
        numSI=numXI1*numXI2
        call MPI_ALLREDUCE(EUVs,EUV,numSI,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,icomm,ierrmpi)
        call MPI_ALLREDUCE(Dpls,Dpl,numSI,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,icomm,ierrmpi)
        do ix1=1,numXI1
          do ix2=1,numXI2
            if (EUV(ix1,ix2)<smalldouble) EUV(ix1,ix2)=zero
            if(EUV(ix1,ix2)/=0) then
              Dpl(ix1,ix2)=-(Dpl(ix1,ix2)/EUV(ix1,ix2))*unitv
            else
              Dpl(ix1,ix2)=0.d0
            endif
            if (abs(Dpl(ix1,ix2))<smalldouble) Dpl(ix1,ix2)=zero
            wI(ix1,ix2,1)=EUV(ix1,ix2)
            wI(ix1,ix2,2)=Dpl(ix1,ix2)
          enddo
        enddo

        xI1=xI1/arcsec
        dxI1=dxI1/arcsec
        xI2=xI2/arcsec
        dxI2=dxI2/arcsec
        call output_data(qunit,xI1,xI2,dxI1,dxI2,wI,numXI1,numXI2,numWI,datatype)
        deallocate(wI,EUV,EUVs,Dpl,Dpls)
      endif

      if (datatype=='image_sxr') then
        numWI=1
        allocate(wI(numXI1,numXI2,numWI))
        allocate(SXR(numXI1,numXI2),SXRs(numXI1,numXI2))
        wI=zero
        SXR=zero
        SXRs=zero
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call integrate_SXR_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,SXRs)
        enddo
        numSI=numXI1*numXI2
        call MPI_ALLREDUCE(SXRs,SXR,numSI,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,icomm,ierrmpi)

        do ix1=1,numXI1
          do ix2=1,numXI2
            if (SXR(ix1,ix2)<smalldouble) SXR(ix1,ix2)=zero
          enddo
        enddo
        wI(:,:,1)=SXR(:,:)

        xI1=xI1/arcsec
        dxI1=dxI1/arcsec
        xI2=xI2/arcsec
        dxI2=dxI2/arcsec
        call output_data(qunit,xI1,xI2,dxI1,dxI2,wI,numXI1,numXI2,numWI,datatype)
        deallocate(wI,SXR,SXRs)
      endif

      deallocate(xI1,xI2,dxI1,dxI2)

    end subroutine get_image_inst_resol

    subroutine integrate_SXR_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,SXR)
      integer, intent(in) :: igrid,numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2)
      double precision, intent(in) :: dxI
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: SXR(numXI1,numXI2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&)
      double precision :: vloc(1:3),res
      integer :: ixP^L,ixP^D,nSubC^D,iSubC^D
      double precision :: xSubP1,xSubP2,dxSubP,xerf^L,fluxsubC
      double precision :: xSubC(1:3),dxSubC^D,xCent1,xCent2

      double precision :: sigma_PSF,RHESSI_rsl,sigma0,factor
      double precision :: arcsec,pixel

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

      allocate(flux(ixI^S))
      ! get local SXR flux
      call get_SXR(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux,emin_sxr,emax_sxr)

      ! integrate emission
      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif
      RHESSI_rsl=2.3d0      
      sigma_PSF=1.d0
      pixel=RHESSI_rsl*arcsec
      sigma0=sigma_PSF*pixel
      {do ix^D=ixOmin^D,ixOmax^D\}
        ^D&nSubC^D=1;
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI1(^D))/(dxI/4.d0)));
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI2(^D))/(dxI/4.d0)));
        ^D&dxSubC^D=ps(igrid)%dx(ix^DD,^D)/nSubC^D;
        ! dividing a cell to several parts to get more accurate integrating values
        {do iSubC^D=1,nSubC^D\}
          ^D&xSubC(^D)=ps(igrid)%x(ix^DD,^D)-half*ps(igrid)%dx(ix^DD,^D)+(iSubC^D-half)*dxSubC^D;
          fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length**3
          ! mapping the 3D coordinate to location at the image
          call dot_product_loc(xSubC,vec_xI1,xCent1)
          call dot_product_loc(xSubC,vec_xI2,xCent2)
          ixP1=floor((xCent1-(xI1(1)-half*dxI))/dxI)+1
          ixP2=floor((xCent2-(xI2(1)-half*dxI))/dxI)+1
          ixPmin1=max(1,ixP1-3)
          ixPmax1=min(ixP1+3,numXI1)
          ixPmin2=max(1,ixP2-3)
          ixPmax2=min(ixP2+3,numXI2)
          do ixP1=ixPmin1,ixPmax1
            do ixP2=ixPmin2,ixPmax2
              xerfmin1=((xI1(ixP1)-half*dxI)-xCent1)/(sqrt(2.d0)*sigma0)
              xerfmax1=((xI1(ixP1)+half*dxI)-xCent1)/(sqrt(2.d0)*sigma0)
              xerfmin2=((xI2(ixP2)-half*dxI)-xCent2)/(sqrt(2.d0)*sigma0)
              xerfmax2=((xI2(ixP2)+half*dxI)-xCent2)/(sqrt(2.d0)*sigma0)
              factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
              SXR(ixP1,ixP2)=SXR(ixP1,ixP2)+fluxSubC*factor
            enddo !ixP2
          enddo !ixP1
        {enddo\} !iSubC
      {enddo\} !ix

      deallocate(flux)
    end subroutine integrate_SXR_inst_resol

    subroutine integrate_EUV_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,fl,EUV,Dpl)
      integer, intent(in) :: igrid,numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2)
      double precision, intent(in) :: dxI
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: EUV(numXI1,numXI2),Dpl(numXI1,numXI2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&),v(:^D&),rho(:^D&)
      double precision :: vloc(1:3),res
      integer :: ixP^L,ixP^D,nSubC^D,iSubC^D
      double precision :: xSubP1,xSubP2,dxSubP,xerf^L,fluxsubC
      double precision :: xSubC(1:3),dxSubC^D,xCent1,xCent2

      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: lineCent
      double precision :: sigma_PSF,spaceRsl,wlRsl,sigma0,factor,wslit
      double precision :: unitv,arcsec,pixel

      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

      allocate(flux(ixI^S),v(ixI^S),rho(ixI^S))
      ! get local EUV flux and velocity
      call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
      call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
      {do ix^D=ixOmin^D,ixOmax^D\}
        do j=1,3
          vloc(j)=ps(igrid)%w(ix^D,iw_mom(j))/rho(ix^D)
        enddo
        call dot_product_loc(vloc,vec_LOS,res)
        v(ix^D)=res
      {enddo\}
      deallocate(rho)

      ! integrate emission
      call get_line_info(wavelength,ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF,wslit)
      if (SI_unit) then
        arcsec=7.25d5/unit_length
      else
        arcsec=7.25d7/unit_length
      endif
      pixel=spaceRsl*arcsec
      sigma0=sigma_PSF*pixel
      {do ix^D=ixOmin^D,ixOmax^D\}
        ^D&nSubC^D=1;
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI1(^D))/(dxI/2.d0)));
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI2(^D))/(dxI/2.d0)));
        ^D&dxSubC^D=ps(igrid)%dx(ix^DD,^D)/nSubC^D;
        if (SI_unit) then
          fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length*1.d2/dxI/dxI  ! DN s^-1
        else
          fluxSubC=flux(ix^D)*dxSubC1*dxSubC2*dxSubC3*unit_length/dxI/dxI  ! DN s^-1
        endif
        ! dividing a cell to several parts to get more accurate integrating values
        {do iSubC^D=1,nSubC^D\}
          ^D&xSubC(^D)=ps(igrid)%x(ix^DD,^D)-half*ps(igrid)%dx(ix^DD,^D)+(iSubC^D-half)*dxSubC^D;
          ! mapping the 3D coordinate to location at the image
          call dot_product_loc(xSubC,vec_xI1,xCent1)
          call dot_product_loc(xSubC,vec_xI2,xCent2)
          ! distribution at nearby pixels
          ixP1=floor((xCent1-(xI1(1)-half*dxI))/dxI)+1
          ixP2=floor((xCent2-(xI2(1)-half*dxI))/dxI)+1
          ixPmin1=max(1,ixP1-3)
          ixPmax1=min(ixP1+3,numXI1)
          ixPmin2=max(1,ixP2-3)
          ixPmax2=min(ixP2+3,numXI2)
          do ixP1=ixPmin1,ixPmax1
            do ixP2=ixPmin2,ixPmax2
              xerfmin1=((xI1(ixP1)-half*dxI)-xCent1)/(sqrt(2.d0)*sigma0) 
              xerfmax1=((xI1(ixP1)+half*dxI)-xCent1)/(sqrt(2.d0)*sigma0)
              xerfmin2=((xI2(ixP2)-half*dxI)-xCent2)/(sqrt(2.d0)*sigma0)
              xerfmax2=((xI2(ixP2)+half*dxI)-xCent2)/(sqrt(2.d0)*sigma0)
              factor=(erfc(xerfmin1)-erfc(xerfmax1))*(erfc(xerfmin2)-erfc(xerfmax2))/4.d0
              EUV(ixP1,ixP2)=EUV(ixP1,ixP2)+fluxSubC*factor
              Dpl(ixP1,ixP2)=Dpl(ixP1,ixP2)+fluxSubC*factor*v(ix^D)
            enddo !ixP2
          enddo !ixP1
        {enddo\} !iSubC
      {enddo\} !ix

      deallocate(flux,v)
    end subroutine integrate_EUV_inst_resol

    subroutine get_image_data_resol(qunit,datatype,fl)
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
      double precision, allocatable :: SXR(:,:),SXRs(:,:),wI(:,:,:)
      double precision, allocatable :: xI1(:),xI2(:),dxI1(:),dxI2(:),dxIi
      integer :: numXI1,numXI2,numSI,numWI
      double precision :: xI^L
      integer :: iigrid,igrid,i,j
      double precision, allocatable :: xIF1(:),xIF2(:),dxIF1(:),dxIF2(:)
      integer :: nXIF1,nXIF2
      double precision :: xIF^L

      double precision :: unitv,arcsec,RHESSI_rsl
      integer :: strtype^D,nstrb^D,nbb^D,nuni^D,nstr^D,bnx^D
      double precision :: qs^D,dxfirst^D,dxmid^D,lenstr^D

      numX1=domain_nx1*2**(refine_max_level-1)
      numX2=domain_nx2*2**(refine_max_level-1)
      numX3=domain_nx3*2**(refine_max_level-1)

      ! parameters for creating table
      if (LOS_phi==0 .and. LOS_theta==90) then
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

      ! integrate EUV flux and get cell average flux for image
      if (datatype=='image_euv') then
        if (SI_unit) then
          unitv=unit_velocity/1.0e3 ! km/s
        else
          unitv=unit_velocity/1.0e5 ! km/s
        endif
        numWI=2
        allocate(wI(nXIF1,nXIF2,numWI))
        allocate(EUVs(nXIF1,nXIF2),EUV(nXIF1,nXIF2))
        allocate(Dpl(nXIF1,nXIF2),Dpls(nXIF1,nXIF2))
        EUVs=0.0d0
        EUV=0.0d0     
        Dpl=0.d0
        Dpls=0.d0 
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call integrate_EUV_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,EUVs,Dpls)
        enddo
        numSI=nXIF1*nXIF2
        call MPI_ALLREDUCE(EUVs,EUV,numSI,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,icomm,ierrmpi)
        call MPI_ALLREDUCE(Dpls,Dpl,numSI,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,icomm,ierrmpi)
        do ix1=1,nXIF1
          do ix2=1,nXIF2
            if (EUV(ix1,ix2)<smalldouble) EUV(ix1,ix2)=zero
            if(EUV(ix1,ix2)/=0) then
              Dpl(ix1,ix2)=-(Dpl(ix1,ix2)/EUV(ix1,ix2))*unitv
            else
              Dpl(ix1,ix2)=0.d0
            endif
            if (abs(Dpl(ix1,ix2))<smalldouble) Dpl(ix1,ix2)=zero
          enddo
        enddo
        wI(:,:,1)=EUV(:,:)
        wI(:,:,2)=Dpl(:,:)

        call output_data(qunit,xIF1,xIF2,dxIF1,dxIF2,wI,nXIF1,nXIF2,numWI,datatype)
        deallocate(WI,EUV,EUVs,Dpl,Dpls)
      endif

      ! integrate EUV flux and get cell average flux for image
      if (datatype=='image_sxr') then
        if (SI_unit) then
          arcsec=7.25d5
        else
          arcsec=7.25d7
        endif
        RHESSI_rsl=2.3d0
        numWI=1
        allocate(wI(nXIF1,nXIF2,numWI))
        allocate(SXRs(nXIF1,nXIF2),SXR(nXIF1,nXIF2))
        SXRs=0.0d0
        SXR=0.0d0
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call integrate_SXR_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,SXRs)
        enddo
        numSI=nXIF1*nXIF2
        call MPI_ALLREDUCE(SXRs,SXR,numSI,MPI_DOUBLE_PRECISION, &
                           MPI_SUM,icomm,ierrmpi)

        SXR=SXR*(RHESSI_rsl*arcsec)**2
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

    end subroutine get_image_data_resol

    subroutine integrate_EUV_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,EUV,Dpl)
      use mod_global_parameters

      integer, intent(in) :: igrid,nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      double precision, intent(in) :: dxIF1(nXIF1),dxIF2(nXIF2)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: EUV(nXIF1,nXIF2),Dpl(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&),v(:^D&),rho(:^D&)
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

      allocate(flux(ixI^S),v(ixI^S),rho(ixI^S))
      allocate(dxb1(ixI^S),dxb2(ixI^S),dxb3(ixI^S))
      dxb1(ixO^S)=ps(igrid)%dx(ixO^S,1)
      dxb2(ixO^S)=ps(igrid)%dx(ixO^S,2)
      dxb3(ixO^S)=ps(igrid)%dx(ixO^S,3)
      ! get local EUV flux and velocity
      call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,fl,flux)
      call fl%get_rho(ps(igrid)%w,ps(igrid)%x,ixI^L,ixO^L,rho)
      v(ixO^S)=ps(igrid)%w(ixO^S,iw_mom(direction_LOS))/rho(ixO^S)
      deallocate(rho)

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

      deallocate(flux,v,dxb1,dxb2,dxb3,EUVg,Fvg,xg1,xg2,dxg1,dxg2)

    end subroutine integrate_EUV_data_resol

    subroutine integrate_SXR_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,fl,SXR)
      use mod_global_parameters

      integer, intent(in) :: igrid,nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      double precision, intent(in) :: dxIF1(nXIF1),dxIF2(nXIF2)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: SXR(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&)
      double precision, allocatable :: dxb1(:^D&),dxb2(:^D&),dxb3(:^D&)
      double precision, allocatable :: SXRg(:,:),xg1(:),xg2(:),dxg1(:),dxg2(:)
      integer :: levelg,nXg1,nXg2,iXgmin1,iXgmax1,iXgmin2,iXgmax2,rft,iXg^D
      double precision :: SXRt,xc^L,xg^L,r2
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

    end subroutine integrate_SXR_data_resol

    subroutine output_data(qunit,xO1,xO2,dxO1,dxO2,wO,nXO1,nXO2,nWO,datatype)
      ! change the format of data and write data
      use mod_global_parameters

      integer, intent(in) :: qunit,nXO1,nXO2,nWO
      double precision, intent(in) :: dxO1(nxO1),dxO2(nxO2)
      double precision, intent(in) :: xO1(nXO1),xO2(nxO2)
      double precision, intent(in) :: wO(nXO1,nXO2,nWO)
      character(20), intent(in) :: datatype

      integer :: nPiece,nP1,nP2,nC1,nC2,nWC
      integer :: piece_nmax1,piece_nmax2,ix1,ix2,j,ipc,ixc1,ixc2
      double precision, allocatable :: xC(:,:,:,:),wC(:,:,:,:),dxC(:,:,:,:)
      character(len=std_len) :: resolution_type

      ! how many cells in each grid
      if (datatype=='image_euv' .and. resolution_euv=='data') then
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
      else if (datatype=='image_sxr' .and. resolution_sxr=='data') then
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
      else if (datatype=='spectrum_euv' .and. resolution_spectrum=='data') then
        piece_nmax1=16
        if (direction_slit==1) then
          piece_nmax2=block_nx1
        else if (direction_slit==2) then
          piece_nmax2=block_nx2
        else
          piece_nmax2=block_nx3
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

      select case(convert_type)
        case('EIvtuCCmpi','ESvtuCCmpi','SIvtuCCmpi')
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
        case('EIvtiCCmpi','ESvtiCCmpi','SIvtiCCmpi')
          if (convert_type=='EIvtiCCmpi') resolution_type=resolution_euv
          if (convert_type=='ESvtiCCmpi') resolution_type=resolution_spectrum
          if (convert_type=='SIvtiCCmpi') resolution_type=resolution_sxr
          if (sum(stretch_type(:))>0 .and. resolution_type=='data') then
            call mpistop("Error in synthesize emission: vti is not supported for data resolution")
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
      integer :: wholeExtent(1:6),extent(1:6)
      integer :: nP1,nP2,iP1,iP2,iw
      integer :: ixC1,ixC2,ixCmin1,ixCmax1,ixCmin2,ixCmax2

      integer :: filenr
      logical :: fileopen
      character (70) :: subname,wname,vname,nameL,nameS
      character (len=std_len) :: filename
      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: line_center
      double precision :: spatial_rsl,spectral_rsl,sigma_PSF,wslit


      origin(1)=xO1(1)-0.5d0*dxO1(1)
      origin(2)=xO2(1)-0.5d0*dxO2(1)
      origin(3)=zero
      spacing(1)=dxO1(1)
      spacing(2)=dxO2(1)
      spacing(3)=zero
      wholeExtent=zero
      wholeExtent(2)=nXO1
      wholeExtent(4)=nXO2
      nP1=nXO1/nC1
      nP2=nXO2/nC2

      ! get information of emission line
      if (convert_type=='EIvtiCCmpi') then
        call get_line_info(wavelength,ion,mass,logTe,line_center,spatial_rsl,spectral_rsl,sigma_PSF,wslit)
      else if (convert_type=='ESvtiCCmpi') then
        call get_line_info(spectrum_wl,ion,mass,logTe,line_center,spatial_rsl,spectral_rsl,sigma_PSF,wslit)
      endif

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
        if (convert_type=='EIvtiCCmpi' .or. convert_type=='ESvtiCCmpi') then
          write(qunit,'(2a)')'<DataArray type="Float32" Name="logT" ',&
                             'NumberOfTuples="1" format="ascii">'
          write(qunit,*) real(logTe)
          write(qunit,'(a)')'</DataArray>'
        endif
        write(qunit,'(a)')'</FieldData>'
        ! pixel/cell data
        do iP1=1,nP1
          do iP2=1,nP2
            extent=zero
            extent(1)=(iP1-1)*nC1
            extent(2)=iP1*nC1
            extent(3)=(iP2-1)*nC2
            extent(4)=iP2*nC2
            ixCmin1=extent(1)+1
            ixCmax1=extent(2)
            ixCmin2=extent(3)+1
            ixCmax2=extent(4)
            write(qunit,'(a,6(i10),a)') &
                  '<Piece Extent="',extent,'">'
            write(qunit,'(a)')'<CellData>'
            do iw=1,nWO
              ! variable name
              if (convert_type=='EIvtiCCmpi') then
                if (iw==1) write(vname,'(a,i4)') "AIA ",wavelength
                if (iw==2) vname='Doppler velocity'
              else if (convert_type=='SIvtiCCmpi') then
                write(vname,'(a,i2,a,i2,a)') "SXR ",emin_sxr,"-",emax_sxr," keV"
              else if (convert_type=='ESvtiCCmpi') then
                write(vname,'(a,i4)') "spectra ",spectrum_wl
              endif
              write(qunit,'(a,a,a)')&
                '<DataArray type="Float64" Name="',TRIM(vname),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') ((wO(ixC1,ixC2,iw),ixC1=ixCmin1,ixCmax1),ixC2=ixCmin2,ixCmax2)
              write(qunit,'(a)')'</DataArray>'
            enddo
            write(qunit,'(a)')'</CellData>'
            write(qunit,'(a)')'</Piece>'
          enddo
        enddo
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
      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: line_center
      double precision :: spatial_rsl,spectral_rsl,sigma_PSF,wslit

      nP1=nC1+1
      nP2=nC2+1
      np=nP1*nP2
      nc=nC1*nC2
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
      ! get information of emission line
      if (datatype=='image_euv') then
        call get_line_info(wavelength,ion,mass,logTe,line_center,spatial_rsl,spectral_rsl,sigma_PSF,wslit)
      else if (datatype=='spectrum_euv') then
        call get_line_info(spectrum_wl,ion,mass,logTe,line_center,spatial_rsl,spectral_rsl,sigma_PSF,wslit)
      endif

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
        if (datatype=='image_euv' .or. datatype=='spectrum_euv') then
          write(qunit,'(2a)')'<DataArray type="Float32" Name="logT" ',&
                             'NumberOfTuples="1" format="ascii">'
          write(qunit,*) real(logTe)
          write(qunit,'(a)')'</DataArray>'
        endif
        write(qunit,'(a)')'</FieldData>'
        do ixP=1,nPiece
          write(qunit,'(a,i7,a,i7,a)') &
                '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do j=1,nWC
            if (datatype=='image_euv') then
              if (j==1) write(vname,'(a,i4)') "AIA ",wavelength
              if (j==2) vname='Doppler velocity'
            else if (datatype=='image_sxr') then
              write(vname,'(a,i2,a,i2,a)') "SXR ",emin_sxr,"-",emax_sxr," keV"
            else if (datatype=='spectrum_euv') then
              write(vname,'(a,i4)') "spectra ",spectrum_wl
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
              if (datatype=='image_euv' .and. resolution_euv=='data') then
                if (LOS_phi==0 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') 0.d0,xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2)
                else if (LOS_phi==90 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,2),0.d0,xP(ixP,ix1,ix2,1)
                else
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2),0.d0
                endif
              else if (datatype=='image_sxr' .and. resolution_sxr=='data') then
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

    subroutine dot_product_loc(vec_in1,vec_in2,res)
      double precision, intent(in) :: vec_in1(1:3),vec_in2(1:3)
      double precision, intent(out) :: res
      integer :: j

      res=zero
      do j=1,3
        res=res+vec_in1(j)*vec_in2(j)
      enddo

    end subroutine dot_product_loc

    subroutine cross_product_loc(vec_in1,vec_in2,vec_out)
      double precision, intent(in) :: vec_in1(1:3),vec_in2(1:3)
      double precision, intent(out) :: vec_out(1:3)

      vec_out(1)=vec_in1(2)*vec_in2(3)-vec_in1(3)*vec_in2(2)
      vec_out(2)=vec_in1(3)*vec_in2(1)-vec_in1(1)*vec_in2(3)
      vec_out(3)=vec_in1(1)*vec_in2(2)-vec_in1(2)*vec_in2(1)

    end subroutine cross_product_loc

    subroutine init_vectors()
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
        vec_xI1=zero
        vec_xI2=zero
        vec_xI1(1)=1.d0
        vec_xI2(2)=1.d0
      else
        call cross_product_loc(vec_LOS,vec_z,vec_xI1)
        call cross_product_loc(vec_xI1,vec_LOS,vec_xI2)
      endif
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

    end subroutine init_vectors

end module mod_thermal_emission

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

  contains

    subroutine get_line_info(ion,mass,logTe,line_center,spatial_rsl,spectral_rsl,sigma_PSF)
      ! get information of the spectral line
      ! mass: ion mass, unit -- proton mass
      ! logTe: peak temperature of emission line in logarithm
      ! line_center: center wavelength of emission line, unit -- Angstrom (0.1 nm) 
      ! spatial_rsl: spatial resolution of instrument (for image), unit -- arcsec
      ! spectral_rsl: spectral resolution in wagelength of instrument (for spectrum), unit -- Angstrom
      ! sigma_PSF: width of point spread function core (for instrument), unit -- pixel (spatial_rsl)
      use mod_global_parameters

      integer, intent(out) :: mass
      character(len=30), intent(out) :: ion
      double precision, intent(out) :: logTe
      double precision, intent(out) :: line_center
      double precision, intent(out) :: spatial_rsl,spectral_rsl
      double precision, intent(out) :: sigma_PSF

      select case(wavelength)
      case(304)
        ion='He II'
        mass=4
        logTe=4.7
        line_center=303.8
        spatial_rsl=0.6
        spectral_rsl=0.2
        sigma_PSF=0.895
      case(171)
        ion='Fe IX'
        mass=56
        logTe=5.8
        line_center=171.1
        spatial_rsl=0.6
        spectral_rsl=0.2 
        sigma_PSF=1.019
      case(193)
        ion='Fe XXIV'
        mass=56
        logTe=7.3
        line_center=193.5
        spatial_rsl=0.6
        spectral_rsl=0.2
        sigma_PSF=0.813
      case(211)
        ion='Fe XIV'
        mass=56
        logTe=6.3
        line_center=211.3
        spatial_rsl=0.6
        spectral_rsl=0.2
        sigma_PSF=0.913
      case(335)
        ion='Fe XVI'
        mass=56
        logTe=6.4
        line_center=335.4
        spatial_rsl=0.6
        spectral_rsl=0.2
        sigma_PSF=1.019
      case(94)
        ion='Fe XVIII'
        mass=56
        logTe=6.8
        line_center=93.9
        spatial_rsl=0.6
        spectral_rsl=0.2
        sigma_PSF=1.025
      case(131)
        ion='Fe XXI'
        mass=56
        logTe=7.0
        line_center=131.0
        spatial_rsl=0.6
        spectral_rsl=0.2
        sigma_PSF=0.984
      case default
        call mpistop("No information about this line")
      end select
    end subroutine get_line_info
    
    subroutine get_EUV(wl,ixI^L,ixO^L,w,x,flux)
      ! calculate the local emission intensity of given EUV line (optically thin)
      ! wavelength is the wave length of the emission line
      ! unit [DN cm^-1 s^-1 pixel^-1]
      use mod_global_parameters

      integer, intent(in) :: wl
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(out) :: flux(ixI^S)

      integer :: n_table
      double precision, allocatable :: t_table(:),f_table(:)
      integer :: ix^D,iTt,i
      double precision :: pth(ixI^S),Te(ixI^S),Ne(ixI^S)
      double precision :: logT,logGT

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
      case default
        allocate(t_table(1))
        allocate(f_table(1))
        call mpistop("This wavelength is unknown")
      end select
      do i=1,n_table
        if(f_table(i) .gt. 1.d-99) then
          f_table(i)=log10(f_table(i))
        else
          f_table(i)=-99.d0
        endif
      enddo 
      call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)
      Te(ixO^S)=pth(ixO^S)/w(ixO^S,iw_rho)*unit_temperature
      Ne(ixO^S)=w(ixO^S,iw_rho)*unit_numberdensity
      flux(ixO^S)=Ne(ixO^S)**2
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
      deallocate(t_table,f_table)
    end subroutine get_EUV

    subroutine get_SXR(ixI^L,ixO^L,w,x,flux,El,Eu)
      !synthesize thermal SXR from El keV to Eu keV
      use mod_global_parameters

      integer, intent(in)           :: ixI^L,ixO^L
      integer, intent(in)           :: El,Eu
      double precision, intent(in)  :: x(ixI^S,1:ndim)
      double precision, intent(in)  :: w(ixI^S,nw)
      double precision, intent(out) :: flux(ixI^S)

      integer :: ix^D,ixO^D
      integer :: iE,numE
      double precision :: I0,kb,keV,dE,Ei
      double precision :: pth(ixI^S),Te(ixI^S),kbT(ixI^S)
      double precision :: Ne(ixI^S),gff(ixI^S),fi(ixI^S)
      double precision :: EM(ixI^S)

      I0=1.07d-42    ! photon flux index at 1AU [cm^-2 s^-1 keV^-1]
      kb=const_kb
      keV=1.0d3*const_ev
      dE=0.1
      numE=floor((Eu-El)/dE)
      call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)
      Te(ixO^S)=pth(ixO^S)/w(ixO^S,iw_rho)*unit_temperature
      Ne(ixO^S)=w(ixO^S,iw_rho)*unit_numberdensity
      kbT(ixO^S)=kb*Te(ixO^S)/keV
      flux(ixO^S)=0.0d0
      EM(ixO^S)=(Ne(ixO^S))**2
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

    subroutine get_GOES_SXR_flux(xbox^L,eflux)
      !get GOES SXR 1-8A flux observed at 1AU from given box [w/m^2]
      use mod_global_parameters

      double precision, intent(in) :: xbox^L
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
        call get_GOES_flux_grid(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,ps(igrid)%dvolume(ixI^S),xbox^L,xb^L,eflux_grid)
        eflux_pe=eflux_pe+eflux_grid
      enddo
      call MPI_ALLREDUCE(eflux_pe,eflux,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    end subroutine get_GOES_SXR_flux

    subroutine get_GOES_flux_grid(ixI^L,ixO^L,w,x,dV,xbox^L,xb^L,eflux_grid)
      use mod_global_parameters

      integer, intent(in)           :: ixI^L,ixO^L
      double precision, intent(in)  :: x(ixI^S,1:ndim),dV(ixI^S)
      double precision, intent(in)  :: w(ixI^S,nw)
      double precision, intent(in)  :: xbox^L,xb^L
      double precision, intent(out) :: eflux_grid

      integer :: ix^D,ixO^D,ixb^L
      integer :: iE,numE,j,inbox
      double precision :: I0,kb,keV,dE,Ei,El,Eu,A_cgs
      double precision :: pth(ixI^S),Te(ixI^S),kbT(ixI^S)
      double precision :: Ne(ixI^S),EM(ixI^S)
      double precision :: gff,fi,erg_SI

      ! check whether the grid is inside given box
      inbox=zero
      {if (xbmin^D<xboxmax^D .and. xbmax^D>xboxmin^D) inbox=inbox+1\}

      if (inbox==ndim) then
        ! indexes for cells inside given box
        ^D&ixbmin^D=ixOmin^D;
        ^D&ixbmax^D=ixOmax^D;
        {if (xbmax^D>xboxmax^D) ixbmax^D=ixOmax^D-ceiling((xbmax^D-xboxmax^D)/dxlevel(^D))\}
        {if (xbmin^D<xboxmin^D) ixbmin^D=ceiling((xboxmin^D-xbmin^D)/dxlevel(^D))+ixOmin^D\}

        I0=1.07d-38 ! photon flux index at 1AU [m^-2 s^-1 keV^-1]
        kb=const_kb
        keV=1.0d3*const_ev
        erg_SI=1.d-7
        A_cgs=1.d-8 ! Angstrom
        El=const_h*const_c/(8.d0*A_cgs)/keV ! 8 A
        Eu=const_h*const_c/(1.d0*A_cgs)/keV ! 1 A
        dE=0.1  ! keV
        numE=floor((Eu-El)/dE)
        call phys_get_pthermal(w,x,ixI^L,ixb^L,pth)
        Te(ixb^S)=pth(ixb^S)/w(ixb^S,iw_rho)*unit_temperature
        Ne(ixb^S)=w(ixb^S,iw_rho)*unit_numberdensity
        kbT(ixb^S)=kb*Te(ixb^S)/keV
        EM(ixb^S)=(I0*(Ne(ixb^S))**2)*dV(ixb^S)*unit_length**3
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
    subroutine get_EUV_image(qunit)
      use mod_global_parameters

      integer, intent(in) :: qunit
      character(20) :: datatype

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl

      if (mype==0) print *, 'Systhesizing EUV image (observed by SDO/AIA).'
      if (mype==0) write(*,'(a,i3,a)') ' Wavelength: ',wavelength,' Angstrom'
      if (mype==0) print *, 'Unit of EUV flux: DN s^-1 pixel^-1'
      call get_line_info(ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF)
      if (mype==0) write(*,'(a,f5.1,a,f5.1,a)') ' Pixel: ',spaceRsl*725.0,' km x ',spaceRsl*725.0, ' km'

      datatype='image_euv'

      if (image_euv) then
        if (resolution_euv=='data') then
          if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d5,' km'
          if (LOS_theta==0 .and. LOS_phi==90) then
            call get_image_data_resol(qunit,datatype)
          else if (LOS_theta==90 .and. LOS_phi==90) then
            call get_image_data_resol(qunit,datatype)
          else if (LOS_phi==0) then
            call get_image_data_resol(qunit,datatype)
          else
            call MPISTOP('ERROR: Wrong LOS for synthesizing emission!')
          endif
        else if (resolution_euv=='instrument') then
          if (mype==0) print *, 'Unit of length: arcsec (~725 km)'
          call get_image_inst_resol(qunit,datatype)
        else
          call MPISTOP('ERROR: Wrong resolution type')
        endif
      endif
      
    end subroutine get_EUV_image

    subroutine get_SXR_image(qunit)
      use mod_global_parameters

      integer, intent(in) :: qunit
      character(20) :: datatype

      if (mype==0) print *, 'Systhesizing SXR image (observed at 1 AU).'
      if (mype==0) write(*,'(a,i2,a,i2,a)') ' Passband: ',emin_sxr,' - ',emax_sxr,' keV'
      if (mype==0) print *, 'Unit of SXR flux: photons cm^-2 s^-1 pixel^-1'
      if (mype==0) write(*,'(a,f6.1,a,f6.1,a)') ' Pixel: ',2.3*725.0,' km x ',2.3*725.0, ' km'

      datatype='image_sxr'

      if (image_sxr) then
        if (resolution_sxr=='data') then
          !if (mype==0) print *, 'Unit of length: ',unit_length, ' cm'
          if (mype==0) write(*,'(a,f8.1,a)') ' Unit of length: ',unit_length/1.d5,' km'
          if (LOS_theta==0 .and. LOS_phi==90) then
            call get_image_data_resol(qunit,datatype)
          else if (LOS_theta==90 .and. LOS_phi==90) then
            call get_image_data_resol(qunit,datatype)
          else if (LOS_phi==0) then
            call get_image_data_resol(qunit,datatype)
          else
            call MPISTOP('ERROR: Wrong LOS for synthesizing emission!')
          endif
        else if (resolution_sxr=='instrument') then
          if (mype==0) print *, 'Unit of length: arcsec (~725 km)'
          call get_image_inst_resol(qunit,datatype)
        else
          call MPISTOP('ERROR: Wrong resolution type')
        endif
      endif

    end subroutine get_SXR_image

    subroutine get_image_inst_resol(qunit,datatype)
      ! integrate emission flux along line of sight (LOS) 
      ! in a 3D simulation box and get a 2D EUV image
      use mod_global_parameters
      use mod_constants

      integer, intent(in) :: qunit
      character(20), intent(in) :: datatype

      integer :: ix^D,numXI1,numXI2,numWI
      double precision :: xImin1,xImax1,xImin2,xImax2,xIcent1,xIcent2,dxI
      double precision, allocatable :: xI1(:),xI2(:),dxI1(:),dxI2(:),wI(:,:,:)
      double precision, allocatable :: EUVs(:,:),EUV(:,:),Dpls(:,:),Dpl(:,:)
      double precision, allocatable :: SXRs(:,:),SXR(:,:)
      double precision :: vec_LOS(1:3),vec_xI1(1:3),vec_xI2(1:3)
      double precision :: vec_temp1(1:3),vec_temp2(1:3)
      double precision :: vec_z(1:3),vec_cor(1:3),xI_cor(1:2)
      double precision :: res,LOS_psi

      integer :: mass
      character (30) :: ion
      double precision :: logTe,lineCent,sigma_PSF,spaceRsl,wlRsl
      double precision :: unitv,arcsec,RHESSI_rsl,pixel
      integer :: iigrid,igrid,i,j,numSI

      ! vectors for image coordinate
      vec_LOS(1)=-cos(dpi*LOS_theta/180.d0)*sin(dpi*LOS_phi/180.d0)
      vec_LOS(2)=-sin(dpi*LOS_theta/180.d0)*sin(dpi*LOS_phi/180.d0)
      vec_LOS(3)=-cos(dpi*LOS_phi/180.d0)
      vec_z(:)=zero
      vec_z(3)=1.d0
      if (LOS_phi==zero) then
        vec_xI1=zero
        vec_xI2=zero
        vec_xI1(1)=1.d0
        vec_xI2(2)=1.d0
      else
        call cross_product_loc(vec_z,vec_LOS,vec_xI1)
        call cross_product_loc(vec_LOS,vec_xI1,vec_xI2)
      endif
      vec_temp1=vec_xI1/sqrt(vec_xI1(1)**2+vec_xI1(2)**2+vec_xI1(3)**2)
      vec_temp2=vec_xI2/sqrt(vec_xI2(1)**2+vec_xI2(2)**2+vec_xI2(3)**2)
      LOS_psi=dpi*image_rotate/180.d0
      vec_xI1=vec_temp1*cos(LOS_psi)-vec_temp2*sin(LOS_psi)
      vec_xI2=vec_temp2*cos(LOS_psi)+vec_temp1*sin(LOS_psi)

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
      unitv=unit_velocity/1.0e5 ! km/s
      arcsec=7.25d7/unit_length
      if (datatype=='image_euv') then
        call get_line_info(ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF)
        dxI=spaceRsl*arcsec  ! intrument resolution of image
      else
        RHESSI_rsl=2.3d0
        dxI=RHESSI_rsl*arcsec
      endif
      numXI1=ceiling((xImax1-xIcent1)/dxI)
      xImin1=xIcent1-numXI1*dxI
      xImax1=xIcent1+numXI1*dxI
      numXI1=numXI1*2
      numXI2=ceiling((xImax2-xIcent2)/dxI)
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
          call integrate_EUV_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,EUVs,Dpls,vec_LOS,vec_xI1,vec_xI2)
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
          call integrate_SXR_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,SXRs,vec_LOS,vec_xI1,vec_xI2)
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

    subroutine integrate_SXR_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,SXR,vec_LOS,vec_xI1,vec_xI2)
      integer, intent(in) :: igrid,numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2)
      double precision, intent(in) :: dxI
      double precision, intent(out) :: SXR(numXI1,numXI2)
      double precision, intent(in) ::  vec_LOS(1:3),vec_xI1(1:3),vec_xI2(1:3)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&)
      double precision :: vloc(1:3),res
      integer :: ixP^L,ixP^D,ixIF^D
      integer :: nSubP,iSubP1,iSubP2,nSubC^D,iSubC^D
      double precision :: xSubP1,xSubP2,dxSubP,SXRsub
      double precision :: xSubC(1:3),dxSubC^D,xCent1,xCent2,r2

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
      call get_SXR(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,flux,emin_sxr,emax_sxr)

      ! integrate emission
      arcsec=7.25d7/unit_length
      RHESSI_rsl=2.3d0      
      sigma_PSF=1.d0
      pixel=RHESSI_rsl*arcsec
      sigma0=sigma_PSF*pixel
      nSubP=4   ! dividing one pixel to nSubP*nSubP parts in calculation
      dxSubP=dxI/nSubP
      {do ix^D=ixOmin^D,ixOmax^D\}
        ^D&nSubC^D=1;
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI1(^D))/dxSubP));
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI2(^D))/dxSubP));
        ^D&dxSubC^D=ps(igrid)%dx(ix^DD,^D)/nSubC^D;
        ! dividing a cell to several parts to get more accurate integrating values
        {do iSubC^D=1,nSubC^D\}
          ^D&xSubC(^D)=ps(igrid)%x(ix^DD,^D)-half*ps(igrid)%dx(ix^DD,^D)+(iSubC^D-half)*dxSubC^D;
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
              ! dividing a pixel to several parts to get more accurate integrating values
              do iSubP1=1,nSubP
                do iSubP2=1,nSubP
                  xSubP1=xI1(ixP1)-half*dxI+(iSubP1-half)*dxSubP
                  xSubP2=xI2(ixP2)-half*dxI+(iSubP2-half)*dxSubP
                  r2=(xCent1-xSubP1)**2+(xCent2-xSubP2)**2
                  factor=dxSubC1*dxSubC2*dxSubC3*unit_length**3
                  factor=factor*exp(-r2/(2*sigma0**2))/6.25
                  factor=factor*dxSubP*dxSubP/dxI/dxI
                  SXRsub=flux(ix^D)*factor
                  SXR(ixP1,ixP2)=SXR(ixP1,ixP2)+SXRsub
                enddo !iSubP2
              enddo !iSubP1
            enddo !ixP2
          enddo !ixP1
        {enddo\} !iSubC
      {enddo\} !ix

      deallocate(flux)
    end subroutine integrate_SXR_inst_resol

    subroutine integrate_EUV_inst_resol(igrid,numXI1,numXI2,xI1,xI2,dxI,EUV,Dpl,vec_LOS,vec_xI1,vec_xI2)
      integer, intent(in) :: igrid,numXI1,numXI2
      double precision, intent(in) :: xI1(numXI1),xI2(numXI2)
      double precision, intent(in) :: dxI
      double precision, intent(out) :: EUV(numXI1,numXI2),Dpl(numXI1,numXI2)
      double precision, intent(in) ::  vec_LOS(1:3),vec_xI1(1:3),vec_xI2(1:3)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&),v(:^D&)
      double precision :: vloc(1:3),res
      integer :: ixP^L,ixP^D,ixIF^D
      integer :: nSubP,iSubP1,iSubP2,nSubC^D,iSubC^D
      double precision :: xSubP1,xSubP2,dxSubP,EUVsub,Dplsub
      double precision :: xSubC(1:3),dxSubC^D,xCent1,xCent2,r2

      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: lineCent
      double precision :: sigma_PSF,spaceRsl,wlRsl,sigma0,factor
      double precision :: unitv,arcsec,pixel


      ^D&ixOmin^D=ixmlo^D\
      ^D&ixOmax^D=ixmhi^D\
      ^D&ixImin^D=ixglo^D\
      ^D&ixImax^D=ixghi^D\
      ^D&xbmin^D=rnode(rpxmin^D_,igrid)\
      ^D&xbmax^D=rnode(rpxmax^D_,igrid)\

      allocate(flux(ixI^S),v(ixI^S))
      ! get local EUV flux and velocity
      call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,flux)
      {do ix^D=ixOmin^D,ixOmax^D\}
        do j=1,3
          vloc(j)=ps(igrid)%w(ix^D,iw_mom(j))/ps(igrid)%w(ix^D,iw_rho)
        enddo
        call dot_product_loc(vloc,vec_LOS,res)
        v(ix^D)=res
      {enddo\}

      ! integrate emission
      call get_line_info(ion,mass,logTe,lineCent,spaceRsl,wlRsl,sigma_PSF)
      arcsec=7.25d7/unit_length
      pixel=spaceRsl*arcsec
      sigma0=sigma_PSF*pixel
      nSubP=2   ! dividing one pixel to nSubP*nSubP parts in calculation
      dxSubP=dxI/nSubP
      {do ix^D=ixOmin^D,ixOmax^D\}
        ^D&nSubC^D=1;
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI1(^D))/dxSubP));
        ^D&nSubC^D=max(nSubC^D,ceiling(ps(igrid)%dx(ix^DD,^D)*abs(vec_xI2(^D))/dxSubP));
        ^D&dxSubC^D=ps(igrid)%dx(ix^DD,^D)/nSubC^D;
        ! dividing a cell to several parts to get more accurate integrating values
        {do iSubC^D=1,nSubC^D\}
          ^D&xSubC(^D)=ps(igrid)%x(ix^DD,^D)-half*ps(igrid)%dx(ix^DD,^D)+(iSubC^D-half)*dxSubC^D;
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
              ! dividing a pixel to several parts to get more accurate integrating values
              do iSubP1=1,nSubP
                do iSubP2=1,nSubP
                  xSubP1=xI1(ixP1)-half*dxI+(iSubP1-half)*dxSubP
                  xSubP2=xI2(ixP2)-half*dxI+(iSubP2-half)*dxSubP
                  r2=(xCent1-xSubP1)**2+(xCent2-xSubP2)**2
                  factor=dxSubC1*dxSubC2*dxSubC3*unit_length
                  factor=factor*exp(-r2/(2*sigma0**2))/6.25
                  factor=factor*dxSubP*dxSubP/dxI/dxI
                  factor=factor/dxI/dxI
                  EUVsub=flux(ix^D)*factor
                  Dplsub=flux(ix^D)*v(ix^D)*factor
                  EUV(ixP1,ixP2)=EUV(ixP1,ixP2)+EUVsub
                  Dpl(ixP1,ixP2)=Dpl(ixP1,ixP2)+Dplsub
                enddo !iSubP2
              enddo !iSubP1
            enddo !ixP2
          enddo !ixP1
        {enddo\} !iSubC
      {enddo\} !ix

      deallocate(flux,v)
    end subroutine integrate_EUV_inst_resol

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

    subroutine get_image_data_resol(qunit,datatype)
      ! integrate emission flux along line of sight (LOS) 
      ! in a 3D simulation box and get a 2D EUV image
      use mod_global_parameters
      use mod_constants

      integer, intent(in) :: qunit
      character(20), intent(in) :: datatype

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

      ! finest table
      if (LOS_theta==0 .and. LOS_phi==90) then
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
      else if (LOS_theta==90 .and. LOS_phi==90) then
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

      ! initialize image coordinate
      allocate(xIF1(nXIF1),xIF2(nXIF2),dxIF1(nXIF1),dxIF2(nXIF2))

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

      if (datatype=='image_euv') then
        unitv=unit_velocity/1.0e5 ! km/s
        ! integrate flux and get cell center flux for image
        numWI=2
        allocate(wI(nXIF1,nXIF2,numWI))
        allocate(EUVs(nXIF1,nXIF2),EUV(nXIF1,nXIF2))
        allocate(Dpl(nXIF1,nXIF2),Dpls(nXIF1,nXIF2))
        EUVs=0.0d0
        EUV=0.0d0     
        Dpl=0.d0
        Dpls=0.d0 
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call integrate_EUV_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,EUVs,Dpls)
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

      if (datatype=='image_sxr') then
        arcsec=7.25d7
        RHESSI_rsl=2.3d0
        ! integrate flux and get cell center flux for image
        numWI=1
        allocate(wI(nXIF1,nXIF2,numWI))
        allocate(SXRs(nXIF1,nXIF2),SXR(nXIF1,nXIF2))
        SXRs=0.0d0
        SXR=0.0d0
        do iigrid=1,igridstail; igrid=igrids(iigrid);
          call integrate_SXR_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,SXRs)
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

    subroutine integrate_EUV_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,EUV,Dpl)
      use mod_global_parameters

      integer, intent(in) :: igrid,nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      double precision, intent(in) :: dxIF1(nXIF1),dxIF2(nXIF2)
      double precision, intent(out) :: EUV(nXIF1,nXIF2),Dpl(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&),v(:^D&)
      double precision, allocatable :: dxb1(:^D&),dxb2(:^D&),dxb3(:^D&)
      double precision, allocatable :: EUVg(:,:),Fvg(:,:),xg1(:),xg2(:),dxg1(:),dxg2(:)
      integer :: levelg,nXg1,nXg2,iXgmin1,iXgmax1,iXgmin2,iXgmax2,rft,iXg^D
      double precision :: EUVt,Fvt,xc^L,xg^L,r2
      integer :: ixP^L,ixP^D,ixIF^D
      integer :: direction_LOS

      if (LOS_theta==0 .and. LOS_phi==90) then
        direction_LOS=1
      else if (LOS_theta==90 .and. LOS_phi==90) then
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

      allocate(flux(ixI^S),v(ixI^S))
      allocate(dxb1(ixI^S),dxb2(ixI^S),dxb3(ixI^S))
      dxb1(ixO^S)=ps(igrid)%dx(ixO^S,1)
      dxb2(ixO^S)=ps(igrid)%dx(ixO^S,2)
      dxb3(ixO^S)=ps(igrid)%dx(ixO^S,3)
      ! get local EUV flux and velocity
      call get_EUV(wavelength,ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,flux)
      v(ixO^S)=ps(igrid)%w(ixO^S,iw_mom(direction_LOS))/ps(igrid)%w(ixO^S,iw_rho)

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

    subroutine integrate_SXR_data_resol(igrid,nXIF1,nXIF2,xIF1,xIF2,dxIF1,dxIF2,SXR)
      use mod_global_parameters

      integer, intent(in) :: igrid,nXIF1,nXIF2
      double precision, intent(in) :: xIF1(nXIF1),xIF2(nXIF2)
      double precision, intent(in) :: dxIF1(nXIF1),dxIF2(nXIF2)
      double precision, intent(out) :: SXR(nXIF1,nXIF2)

      integer :: ixO^L,ixO^D,ixI^L,ix^D,i,j
      double precision :: xb^L,xd^D
      double precision, allocatable :: flux(:^D&)
      double precision, allocatable :: dxb1(:^D&),dxb2(:^D&),dxb3(:^D&)
      double precision, allocatable :: SXRg(:,:),xg1(:),xg2(:),dxg1(:),dxg2(:)
      integer :: levelg,nXg1,nXg2,iXgmin1,iXgmax1,iXgmin2,iXgmax2,rft,iXg^D
      double precision :: SXRt,xc^L,xg^L,r2
      integer :: ixP^L,ixP^D,ixIF^D
      integer :: direction_LOS

      if (LOS_theta==0 .and. LOS_phi==90) then
        direction_LOS=1
      else if (LOS_theta==90 .and. LOS_phi==90) then
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
      call get_SXR(ixI^L,ixO^L,ps(igrid)%w,ps(igrid)%x,flux,emin_sxr,emax_sxr)

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

      ! how many cells in each grid
      if (datatype=='image_euv' .and. resolution_euv=='data') then
        if (LOS_theta==0 .and. LOS_phi==90) then
          piece_nmax1=block_nx2
          piece_nmax2=block_nx3
        else if (LOS_theta==90 .and. LOS_phi==90) then
          piece_nmax1=block_nx3
          piece_nmax2=block_nx1
        else
          piece_nmax1=block_nx1
          piece_nmax2=block_nx2
        endif
      else if (datatype=='image_sxr' .and. resolution_sxr=='data') then
        if (LOS_theta==0 .and. LOS_phi==90) then
          piece_nmax1=block_nx2
          piece_nmax2=block_nx3
        else if (LOS_theta==90 .and. LOS_phi==90) then
          piece_nmax1=block_nx3
          piece_nmax2=block_nx1
        else
          piece_nmax1=block_nx1
          piece_nmax2=block_nx2
        endif
      else
        piece_nmax1=16
        piece_nmax2=16
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
      call write_image(qunit,xC,wC,dxC,nPiece,nC1,nC2,nWC,datatype)
      deallocate(xC,wC)
    end subroutine output_data

    subroutine write_image(qunit,xC,wC,dxC,nPiece,nC1,nC2,nWC,datatype)
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
      double precision :: spatial_rsl,spectral_rsl,sigma_PSF

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
        call get_line_info(ion,mass,logTe,line_center,spatial_rsl,spectral_rsl,sigma_PSF)
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
        if (datatype=='image_euv') then
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
            else
              write(vname,'(a,i2,a,i2,a)') "SXR ",emin_sxr,"-",emax_sxr," keV"
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
              if (datatype=='image_euv' .and. resolution_euv=='instrument') then
                write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2),0.d0
              else if (datatype=='image_sxr' .and. resolution_sxr=='instrument') then
                write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2),0.d0
              else
                if (LOS_theta==0 .and. LOS_phi==90) then
                  write(qunit,'(3(1pe14.6))') 0.d0,xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2)
                else if (LOS_theta==90 .and. LOS_phi==90) then
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,2),0.d0,xP(ixP,ix1,ix2,1)
                else
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,2),0.d0
                endif
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
    end subroutine write_image
  }
end module mod_thermal_emission

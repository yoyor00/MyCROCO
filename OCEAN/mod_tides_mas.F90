!============================================================================================
! Reference : Coastal Tides book by Simon (SHOM)
! https://diffusion.shom.fr/coastal-tides-version-anglaise-de-la-maree-oceanique-cotiere.html
!=============================================================================================
#include "cppdefs.h"
module mod_tides_mas
#ifdef TIDES_MAS
  use param
  integer, parameter :: nmarmax=144
  integer, dimension(nmarmax)            :: nu
  integer, dimension(30,8)               :: nomd
  real, dimension(30,8)                  :: fr1,fr2
  character(len=100), dimension(nmarmax) :: nommar
  logical, dimension(nmarmax)            :: l_presence
  integer, dimension(nmarmax)            :: doodson
  real, dimension(nmarmax)               :: omega
  real, dimension(nmarmax)               :: frequence,equitide

  integer :: inivm,isa,issa,imsm,imm,imsf,imf,imfm,i2q1,isigma1
  integer :: iq1,irho1,io1
  integer :: ims1,imp1,ia19,im1,ikhi1,ipi1,i2ns2,i3m2s2,ioq2mnk2,i33,imns2
  integer :: ip1,is1,ik1,ipsi1,iphi1,itheta1,ij1,iso1,ioo1,ikq1,i2mn2s2
  integer :: imnus2,i2mk2,i2n22nm2,imu2ms2,i39,in2,inu2,iop2msk2
  integer :: i43,imsk2,im2,imks2,imks22,ilambda2,i2mk3,im3,iso3
  integer :: ims3,imk3,ia87
  integer :: il22mn2,inkm2,it2,is2,ir22,ik2,imsn2,ikj2mkn2,i2sm2
  integer :: iskm2,imq3
  integer :: isp3,is3,isk3,ik3,i2nms4,i2mmus4,i2mns4,i2mnus4,i3mk4,in4
  integer :: i3ms4,imn4,imnu4,i2msk4,im4,i2msn4,i2mkn4,is4,isk4,i3mnk6
  integer :: i2mks4,isn4,i3mn4ml4,ink4,i2smk4,imt4,ims4,imk4,i2snm4
  integer :: i3mns6,i3mnus6,i4mk6,i2nm6,i4ms6,i2mn6,i2mnu6,i3msk6,im6
  integer :: i3mks6,imsn6,i4mn2ml6,imnk6,i2mt6,i3mnu8,i2msk8,im8,i4mks8
  integer :: i2ms6,i2mk6,i2sn6,i3msn6,i3mkn6,i2sm6,imsk6,i2mns8,i2mn8,i3mn8
  integer :: i2msn8,i3ml8,i2mnk8,i3ms8,i3mk8,i2smn8,i4msn8
  integer :: imsnk8,i4mnk8,i2ms8
  integer :: i4mn10,im10,i3msn10,i4ms10,i4mk10,i5msn10,i2msnk10,i3m2s10
  integer :: i2msk82,i5mns10,i3m2n10

  integer, dimension(0:10), parameter :: esp_corres= &
    (/1,2,3,4,5,0,6,0,7,0,8/)

!tidepredic
  integer :: jourmem=-1
  integer :: mas_ijour=-1,mas_imois=-1,mas_ian=-1
  real    :: mas_tempis=0.0
  logical :: mas_is_new_day=.false.
  integer, dimension(8), save :: num
  integer, dimension(144), parameter :: ide=(/ &
      55555,   56555,   57555,   63655,   65455,   73555,   75555,  125755, &
     127555,  135655,  137455,  145555,  146555,  147555,  153655,  155655, &
     157455,  162556,  163555,  164555,  165555,  166554,  167555,  173655, &
     175455,  183555,  185555,  195455,  209655,  217755,  219555,  225655, &
     225855,  227655,  229455,  235555,  235755,  237555,  238554,  245655, &
     247455,  253555,  253755,  254556,  255555,  256554,  257555,  263655, &
     265455,  265655,  272556,  273555,  274554,  275555,  283455,  285455, &
     291555,  293555,  335655,  345555,  355555,  363555,  364555,  365555, &
     375555,  381555,  382555,  383555,  385555,  417755,  419555,  427655, &
     429455,  435555,  435755,  437555,  445655,  447455,  453555,  455555, &
     457555,  463655,  465455,  465655,  471555,  472556,  473555,  475555, &
     481655,  483455,  485455,  491555,  493555,  625655,  627655,  629455, &
     635555,  635755,  637555,  645655,  647455,  653555,  655555,  657555, &
     663655,  665455,  665655,  672556,  673555,  675555,  681655,  683455, &
     685455,  691555,  693555,  827655,  835755,  845655,  847455,  853555, &
     855555,  857555,  863655,  865455,  865655,  873555,  875555,  881655, &
     883455,  883655,  885455,  891555,  893555, 1027655, 1035755, 1045655, &
    1055555, 1063655, 1073555, 1075555, 1083455, 1083655, 1091555,   85455 /)
  integer, dimension(144) :: rr1=(/ &
       0,    0,    0,  -72,  -66,   72,    0,  189,  189,  189,  188,  189, &
       0,    0,  227,    0,    0,    0,    0,    0,   20,    0,    0,    0, &
       0,    0,    0,    0, -112,  -75, -112,  279,    0,  -75,  -75,  279, &
     -37,  -37,    0,  -37,  -37,  260,    0,    0,  -37,    0,  -37,  -37, &
     -37,  -37,    0,    2,    0,    0,  -37,  -37,    0,    0,  152,   61, &
     -56,  189,    0,    0,    0,    0,    0,   20,   20, -112, -112, -112, &
    -112,  186,  -75, -112,  -75,  -75,  223,  -75,  -75,  -37, -112,  -37, &
     261,  -37,  -37,  -37,  -37,  -75,  -75,    0,    0,  148, -149, -149, &
     148, -112, -149, -112, -112,  186, -112, -112,  -75, -149,  -75,  -75, &
     -75,  -75,  -37, -112, -112,  -37,  -37, -112, -149, -149, -149,  148, &
    -149, -149, -112, -149, -112, -112, -112,  -75, -149,  -75, -149,  -75, &
     -75, -224, -186, -186, -186, -149, -149, -149, -186, -112, -112,   41 /)
  integer, dimension(144) :: rr2=(/ &
       0,    0,    0,  -65,  -65,  -64,  414,    0,    0,    0,    0,    0, &
       0,  218,    0, -200, -219,    0,    0,    0, -135,    0,    0,    0, &
    -198, -195, -640, -640,    0,    0,    0,    0,    0,    0,    0,    0, &
       0,    0,    0,    0,    0,    0,    0,    0,    0, -477,  298,    0, &
       0,  260,    0,    0,    0,  298,  -37,  260,  -37,  260,    0,   20, &
       0,    0,    0,  135,  439,    0,    0, -135, -433,    0,    0,    0, &
       0,    0,    0,    0,    0,    0,    0,    0,  298,    0,  -37,  298, &
       0,    0,    0,  298,  -37,  -37,  261,    0,  298,    0,    0,    0, &
       0,    0,    0,    0,    0,    0,    0,  297,    0,  -37,  297,    0, &
       0,  297,    0,  -37,  260,    0,  297,    0,    0,    0,    0,    0, &
       0,  297,    0,    0,  297,    0,  297,    0,  -37,  297,  260,    0, &
     297,    0,    0,    0,    0,    0,    0,  297,  -37,  297,    0,    0 /)
  integer, dimension(144) :: rr0=(/ &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000, -1000, -1000, -1000, -1000,  1000,  1000,  1000, &
    -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,  1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000, -1000, -1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000, -1000, -1000,  1000,  1000, -1000, -1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000, -1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000,  1000, &
     1000,  1000,  1000,  1000 /)

  integer, parameter :: rsh=4
  real, parameter    :: pi=3.141592654
  real, parameter    :: rad2deg=180.0/pi
  real(kind=rsh), parameter :: valmanq=999.
!
! pi/rad2deg are private: callers already get identically-named
! PARAMETERs from scalars.h, and USE-associating both would conflict.
  private :: pi,rad2deg

contains

  subroutine mas_init_waves(Ntides,k_of,nn_of)
!
!  Called once (serial, before the OMP tile loop) from get_tides:
!  ir0/ir1/ir2/nol are local here (not the module-level scratch
!  scalars of the same name) so this stays safe even if ever called
!  from more than one place.
    implicit none
    integer, intent(in)  :: Ntides
    integer, intent(out) :: k_of(:),nn_of(:)
    integer :: ht,k,nn,nol,ir0,ir1,ir2

    ! Find the active waves and store their information
    num(:) = 0
    nomd(:,:) = 0
    fr1(:,:) = 0.0
    fr2(:,:) = 0.0

    do ht=1,Ntides
      if (l_presence(nu(ht))) then
        nol=ide(nu(ht))
        ir0=rr0(nu(ht))
        ir1=rr1(nu(ht))
        ir2=rr2(nu(ht))
        k=esp_corres(int(nol/100000))
        num(k)=num(k)+1
        nn=num(k)

        k_of(ht) = k
        nn_of(ht) = nn

        nomd(nn,k)=nol
        fr1(nn,k)=ir1/float(ir0)
        fr2(nn,k)=ir2/float(ir0)
      endif
    enddo
  end subroutine mas_init_waves


  subroutine tide_data
!
!  Populate the module-level wave catalog (nommar/doodson/frequence/
!  equitide) and each named species index (inivm, isa, ..., i3m2n10)
!  for the 144 MAS harmonics, then derive frequence/omega for all of
!  them. Call once at initialization, before mas_init_waves.
!
    implicit none
    integer              :: i
    integer, dimension(6) :: nd_tmp
    real :: asfo,vathau,valonmlun,valonmsol,valonperl,valonascl
    real :: valonpers

    inivm    =   1 ; call set_wave(inivm    ,'niveau_moyen' ,   55555,0.00000       ,0.0)
    isa      =   2 ; call set_wave(isa      ,'sa'           ,   56555,0.041068639   ,0.0)
    issa     =   3 ; call set_wave(issa     ,'ssa'          ,   57555,0.082137278   ,0.0)
    imsm     =   4 ; call set_wave(imsm     ,'msm'          ,   63655,0.471521089   ,0.0)
    imm      =   5 ; call set_wave(imm      ,'mm'           ,   65455,0.544374694   ,0.0)
    imsf     =   6 ; call set_wave(imsf     ,'msf'          ,   73555,1.015895783   ,0.0)
    imf      =   7 ; call set_wave(imf      ,'mf'           ,   75555,1.098033061   ,0.0)
    imfm     = 144 ; call set_wave(imfm     ,'mfm'          ,   85455,1.62240772    ,0.0)
    i2q1     =   8 ; call set_wave(i2q1     ,'2Q1'          ,  125755,12.854286191  ,0.0)
    isigma1  =   9 ; call set_wave(isigma1  ,'sigma1'       ,  127555,12.927139796  ,0.0)
    iq1      =  10 ; call set_wave(iq1      ,'Q1'           ,  135655,13.398660885  ,0.019256)
    irho1    =  11 ; call set_wave(irho1    ,'rho1'         ,  137455,13.471514490  ,0.0)
    io1      =  12 ; call set_wave(io1      ,'O1'           ,  145555,13.943035578  ,0.100574)
    ims1     =  13 ; call set_wave(ims1     ,'MS1'          ,  146555,13.984104217  ,0.0)
    imp1     =  14 ; call set_wave(imp1     ,'MP1'          ,  147555,14.025172857  ,0.0)
    ia19     =  15 ; call set_wave(ia19     ,'a19'          ,  153655,14.414556667  ,0.0)
    im1      =  16 ; call set_wave(im1      ,'M1'           ,  155655,14.496693945  ,0.0)
    ikhi1    =  17 ; call set_wave(ikhi1    ,'khi1'         ,  157455,14.569547550  ,0.0)
    ipi1     =  18 ; call set_wave(ipi1     ,'pi1'          ,  162556,14.917864683  ,0.0)
    ip1      =  19 ; call set_wave(ip1      ,'p1'           ,  163555,14.958931361  ,0.046843)
    is1      =  20 ; call set_wave(is1      ,'s1'           ,  164555,15.0          ,0.0)
    ik1      =  21 ; call set_wave(ik1      ,'k1'           ,  165555,15.041068639  ,0.141565)
    ipsi1    =  22 ; call set_wave(ipsi1    ,'psi1'         ,  166554,15.082135317  ,0.0)
    iphi1    =  23 ; call set_wave(iphi1    ,'phi1'         ,  167555,15.123205917  ,0.0)
    itheta1  =  24 ; call set_wave(itheta1  ,'theta1'       ,  173655,15.512589728  ,0.0)
    ij1      =  25 ; call set_wave(ij1      ,'J1'           ,  175455,15.585443333  ,0.0)
    iso1     =  26 ; call set_wave(iso1     ,'SO1'          ,  183555,16.056964422  ,0.0)
    ioo1     =  27 ; call set_wave(ioo1     ,'OO1'          ,  185555,16.1391017    ,0.0)
    ikq1     =  28 ; call set_wave(ikq1     ,'KQ1'          ,  195455,16.683476394  ,0.0)
    i2mn2s2  =  29 ; call set_wave(i2mn2s2  ,'2MN2S2'       ,  209655,26.407937959  ,0.0)
    i2ns2    =  30 ; call set_wave(i2ns2    ,'2NS2'         ,  217755,26.879459048  ,0.0)
    i3m2s2   =  31 ; call set_wave(i3m2s2   ,'3M2S2'        ,  219555,26.952312652  ,0.0)
    ioq2mnk2 =  32 ; call set_wave(ioq2mnk2 ,'OQ2MNK2'      ,  225655,27.341696463  ,0.0)
    i33      =  33 ; call set_wave(i33      ,'33'           ,  225855,27.350980136  ,0.0)
    imns2    =  34 ; call set_wave(imns2    ,'MNS2'         ,  227655,27.423833741  ,0.0)
    imnus2   =  35 ; call set_wave(imnus2   ,'MNUS2'        ,  229455,27.496687346  ,0.0)
    i2mk2    =  36 ; call set_wave(i2mk2    ,'2MK2'         ,  235555,27.886071157  ,0.0)
    i2n22nm2 =  37 ; call set_wave(i2n22nm2 ,'2N22NM2'      ,  235755,27.89535483   ,0.0)
    imu2ms2  =  38 ; call set_wave(imu2ms2  ,'MU22MS2'      ,  237555,27.968208435  ,0.0)
    i39      =  39 ; call set_wave(i39      ,'39'           ,  238554,28.009275113  ,0.0)
    in2      =  40 ; call set_wave(in2      ,'N2'           ,  245655,28.439729524  ,0.046398)
    inu2     =  41 ; call set_wave(inu2     ,'nu2'          ,  247455,28.512583129  ,0.0)
    iop2msk2 =  42 ; call set_wave(iop2msk2 ,'OP2MSK2'      ,  253555,28.901966939  ,0.0)
    i43      =  43 ; call set_wave(i43      ,'43'           ,  253755,28.911250612  ,0.0)
    imsk2    =  44 ; call set_wave(imsk2    ,'M_SK_2'       ,  254556,28.94303754   ,0.0)
    im2      =  45 ; call set_wave(im2      ,'M2'           ,  255555,28.984104218  ,0.24230)
    imks2    =  46 ; call set_wave(imks2    ,'M_KS_2'       ,  256554,29.025170895  ,0.0)
    imks22   =  47 ; call set_wave(imks22   ,'MKS2'         ,  257555,29.066241496  ,0.0)
    ilambda2 =  48 ; call set_wave(ilambda2 ,'lambda2'      ,  263655,29.455625306  ,0.0)
    il22mn2  =  49 ; call set_wave(il22mn2  ,'L2_2MN2'      ,  265455,29.528478911  ,0.0)
    inkm2    =  50 ; call set_wave(inkm2    ,'NKM2'         ,  265655,29.537762585  ,0.0)
    it2      =  51 ; call set_wave(it2      ,'T2'           ,  272556,29.958933322  ,0.0)
    is2      =  52 ; call set_wave(is2      ,'S2'           ,  273555,30.0          ,0.112841)
    ir22     =  53 ; call set_wave(ir22     ,'R2'           ,  274554,30.041066678  ,0.0)
    ik2      =  54 ; call set_wave(ik2      ,'K2'           ,  275555,30.082137278  ,0.030704)
    imsn2    =  55 ; call set_wave(imsn2    ,'MSN2'         ,  283455,30.544374694  ,0.0)
    ikj2mkn2 =  56 ; call set_wave(ikj2mkn2 ,'KJ2MKN2'      ,  285455,30.626511972  ,0.0)
    i2sm2    =  57 ; call set_wave(i2sm2    ,'2SM2'         ,  291555,31.015895782  ,0.0)
    iskm2    =  58 ; call set_wave(iskm2    ,'SKM2'         ,  293555,31.098033061  ,0.0)
    imq3     =  59 ; call set_wave(imq3     ,'MQ3'          ,  335655,42.382765102  ,0.0)
    i2mk3    =  60 ; call set_wave(i2mk3    ,'2MK3'         ,  345555,42.927139796  ,0.0)
    im3      =  61 ; call set_wave(im3      ,'M3'           ,  355555,43.476156326  ,0.0)
    iso3     =  62 ; call set_wave(iso3     ,'SO3'          ,  363555,43.943035578  ,0.0)
    ims3     =  63 ; call set_wave(ims3     ,'MS3'          ,  364555,43.984104217  ,0.0)
    imk3     =  64 ; call set_wave(imk3     ,'MK3'          ,  365555,44.025172857  ,0.0)
    ia87     =  65 ; call set_wave(ia87     ,'A87'          ,  375555,44.574189387  ,0.0)
    isp3     =  66 ; call set_wave(isp3     ,'SP3'          ,  381555,44.958931361  ,0.0)
    is3      =  67 ; call set_wave(is3      ,'s3'           ,  382555,45.0          ,0.0)
    isk3     =  68 ; call set_wave(isk3     ,'sk3'          ,  383555,45.041068639  ,0.0)
    ik3      =  69 ; call set_wave(ik3      ,'k3'           ,  385555,45.123205918  ,0.0)
    i2nms4   =  70 ; call set_wave(i2nms4   ,'2NMS4'        ,  417755,55.863563265  ,0.0)
    i2mmus4  =  71 ; call set_wave(i2mmus4  ,'2MMUS4'       ,  419555,55.936416870  ,0.0)
    i2mns4   =  72 ; call set_wave(i2mns4   ,'2MNS4'        ,  427655,56.407937959  ,0.0)
    i2mnus4  =  73 ; call set_wave(i2mnus4  ,'2MNUS4'       ,  429455,56.480791564  ,0.0)
    i3mk4    =  74 ; call set_wave(i3mk4    ,'3MK4'         ,  435555,56.870175374  ,0.0)
    in4      =  75 ; call set_wave(in4      ,'N4'           ,  435755,56.879459047  ,0.0)
    i3ms4    =  76 ; call set_wave(i3ms4    ,'3MS4'         ,  437555,56.952312653  ,0.0)
    imn4     =  77 ; call set_wave(imn4     ,'MN4'          ,  445655,57.423833741  ,0.0)
    imnu4    =  78 ; call set_wave(imnu4    ,'MNU4'         ,  447455,57.496687346  ,0.0)
    i2msk4   =  79 ; call set_wave(i2msk4   ,'2MSK4'        ,  453555,57.886071157  ,0.0)
    im4      =  80 ; call set_wave(im4      ,'M4'           ,  455555,57.968208435  ,0.0)
    i2mks4   =  81 ; call set_wave(i2mks4   ,'2MKS4'        ,  457555,57.968208435  ,0.0)
    isn4     =  82 ; call set_wave(isn4     ,'SN4'          ,  463655,57.968208435  ,0.0)
    i3mn4ml4 =  83 ; call set_wave(i3mn4ml4 ,'3MN4ML4'      ,  465455,57.968208435  ,0.0)
    ink4     =  84 ; call set_wave(ink4     ,'NK4'          ,  465655,57.968208435  ,0.0)
    i2smk4   =  85 ; call set_wave(i2smk4   ,'2SMK4'        ,  471555,57.968208435  ,0.0)
    imt4     =  86 ; call set_wave(imt4     ,'MT4'          ,  472556,57.968208435  ,0.0)
    ims4     =  87 ; call set_wave(ims4     ,'MS4'          ,  473555,57.968208435  ,0.0)
    imk4     =  88 ; call set_wave(imk4     ,'MK4'          ,  475555,57.968208435  ,0.0)
    i2snm4   =  89 ; call set_wave(i2snm4   ,'2SNM4'        ,  481655,57.968208435  ,0.0)
    i2msn4   =  90 ; call set_wave(i2msn4   ,'2MSN4'        ,  483455,57.968208435  ,0.0)
    i2mkn4   =  91 ; call set_wave(i2mkn4   ,'2MKN4'        ,  485455,57.968208435  ,0.0)
    is4      =  92 ; call set_wave(is4      ,'S4'           ,  491555,57.968208435  ,0.0)
    isk4     =  93 ; call set_wave(isk4     ,'SK4'          ,  493555,57.968208435  ,0.0)
    i3mnk6   =  94 ; call set_wave(i3mnk6   ,'3MNK6'        ,  625655,57.968208435  ,0.0)
    i3mns6   =  95 ; call set_wave(i3mns6   ,'3MNS6'        ,  627655,57.968208435  ,0.0)
    i3mnus6  =  96 ; call set_wave(i3mnus6  ,'3MNUS6'       ,  629455,57.968208435  ,0.0)
    i4mk6    =  97 ; call set_wave(i4mk6    ,'4MK6'         ,  635555,57.968208435  ,0.0)
    i2nm6    =  98 ; call set_wave(i2nm6    ,'2NM6'         ,  635755,57.968208435  ,0.0)
    i4ms6    =  99 ; call set_wave(i4ms6    ,'4MS6'         ,  637555,57.968208435  ,0.0)
    i2mn6    = 100 ; call set_wave(i2mn6    ,'2MN6'         ,  645655,57.968208435  ,0.0)
    i2mnu6   = 101 ; call set_wave(i2mnu6   ,'2MNU6'        ,  647455,57.968208435  ,0.0)
    i3msk6   = 102 ; call set_wave(i3msk6   ,'3MSK6'        ,  653555,57.968208435  ,0.0)
    im6      = 103 ; call set_wave(im6      ,'M6'           ,  655555,57.968208435  ,0.0)
    i3mks6   = 104 ; call set_wave(i3mks6   ,'3MKS6'        ,  657555,57.968208435  ,0.0)
    imsn6    = 105 ; call set_wave(imsn6    ,'MSN6'         ,  663655,57.968208435  ,0.0)
    i4mn2ml6 = 106 ; call set_wave(i4mn2ml6 ,'4MN.2ML6'     ,  665455,57.968208435  ,0.0)
    imnk6    = 107 ; call set_wave(imnk6    ,'MNK6'         ,  665655,57.968208435  ,0.0)
    i2mt6    = 108 ; call set_wave(i2mt6    ,'2MT6'         ,  672556,57.968208435  ,0.0)
    i2ms6    = 109 ; call set_wave(i2ms6    ,'2MS6'         ,  673555,57.968208435  ,0.0)
    i2mk6    = 110 ; call set_wave(i2mk6    ,'2MK6'         ,  675555,57.968208435  ,0.0)
    i2sn6    = 111 ; call set_wave(i2sn6    ,'2SN6'         ,  681655,57.968208435  ,0.0)
    i3msn6   = 112 ; call set_wave(i3msn6   ,'3MSN6'        ,  683455,57.968208435  ,0.0)
    i3mkn6   = 113 ; call set_wave(i3mkn6   ,'3MKN6'        ,  685455,57.968208435  ,0.0)
    i2sm6    = 114 ; call set_wave(i2sm6    ,'2SM6'         ,  691555,57.968208435  ,0.0)
    imsk6    = 115 ; call set_wave(imsk6    ,'MSK6'         ,  693555,57.968208435  ,0.0)
    i2mns8   = 116 ; call set_wave(i2mns8   ,'2MNS8'        ,  827655,57.968208435  ,0.0)
    i2mn8    = 117 ; call set_wave(i2mn8    ,'2_MN_8'       ,  835755,57.968208435  ,0.0)
    i3mn8    = 118 ; call set_wave(i3mn8    ,'3MN8'         ,  845655,57.968208435  ,0.0)
    i3mnu8   = 119 ; call set_wave(i3mnu8   ,'3MNU8'        ,  847455,57.968208435  ,0.0)
    i2msk8   = 120 ; call set_wave(i2msk8   ,'2MSK8'        ,  853555,57.968208435  ,0.0)
    im8      = 121 ; call set_wave(im8      ,'M8'           ,  855555,57.968208435  ,0.0)
    i4mks8   = 122 ; call set_wave(i4mks8   ,'4MKS8'        ,  857555,57.968208435  ,0.0)
    i2msn8   = 123 ; call set_wave(i2msn8   ,'2MSN8'        ,  863655,57.968208435  ,0.0)
    i3ml8    = 124 ; call set_wave(i3ml8    ,'3ML8'         ,  865455,57.968208435  ,0.0)
    i2mnk8   = 125 ; call set_wave(i2mnk8   ,'2MNK8'        ,  865655,57.968208435  ,0.0)
    i3ms8    = 126 ; call set_wave(i3ms8    ,'3MS8'         ,  873555,57.968208435  ,0.0)
    i3mk8    = 127 ; call set_wave(i3mk8    ,'3MK8'         ,  875555,57.968208435  ,0.0)
    i2smn8   = 128 ; call set_wave(i2smn8   ,'2SMN8'        ,  881655,57.968208435  ,0.0)
    i4msn8   = 129 ; call set_wave(i4msn8   ,'4MSN8'        ,  883455,57.968208435  ,0.0)
    imsnk8   = 130 ; call set_wave(imsnk8   ,'MSNK8'        ,  883655,57.968208435  ,0.0)
    i4mnk8   = 131 ; call set_wave(i4mnk8   ,'4MNK8'        ,  885455,57.968208435  ,0.0)
    i2ms8    = 132 ; call set_wave(i2ms8    ,'2_MS_8'       ,  891555,57.968208435  ,0.0)
    i2msk82  = 133 ; call set_wave(i2msk82  ,'2MSK8'        ,  893555,57.968208435  ,0.0)
    i5mns10  = 134 ; call set_wave(i5mns10  ,'5MNS10'       , 1027655,57.968208435  ,0.0)
    i3m2n10  = 135 ; call set_wave(i3m2n10  ,'3M2N10'       , 1035755,57.968208435  ,0.0)
    i4mn10   = 136 ; call set_wave(i4mn10   ,'4MN10'        , 1045655,57.968208435  ,0.0)
    im10     = 137 ; call set_wave(im10     ,'M10'          , 1055555,57.968208435  ,0.0)
    i3msn10  = 138 ; call set_wave(i3msn10  ,'3MSN10'       , 1063655,57.968208435  ,0.0)
    i4ms10   = 139 ; call set_wave(i4ms10   ,'4MS10'        , 1073555,57.968208435  ,0.0)
    i4mk10   = 140 ; call set_wave(i4mk10   ,'4MK10'        , 1075555,57.968208435  ,0.0)
    i5msn10  = 141 ; call set_wave(i5msn10  ,'5MSN10'       , 1083455,57.968208435  ,0.0)
    i2msnk10 = 142 ; call set_wave(i2msnk10 ,'2MSNK10'      , 1083655,57.968208435  ,0.0)
    i3m2s10  = 143 ; call set_wave(i3m2s10  ,'3M2S10'       , 1091555,57.968208435  ,0.0)
!
!-------------------------------------------------
!
    vathau=347.80925080
!   vathau=347.8092506160
!   vathau=360./1.03505
    valonmlun=13.17639673
!   valonmlun=13.1763967440
!   valonmlun=360/27.321582
    valonmsol=0.98564734
!   valonmsol=0.985647359999986
!   valonmsol=0.98564736
!   valonmsol=360./365.242199
    valonperl=0.11140408
!   valonperl=360./8.847309/365.242199
    valonascl=0.05295392
!   valonascl=360./18.612904/365.242199
    valonpers=0.00004707
!   valonpers=360./209.4021/36524.2199
!
    do i=1,144
      nd_tmp(6)=doodson(i)
      nd_tmp(5)=doodson(i)/10
      nd_tmp(4)=doodson(i)/100
      nd_tmp(3)=doodson(i)/1000
      nd_tmp(2)=doodson(i)/10000
      nd_tmp(1)=doodson(i)/100000
      nd_tmp(6)=nd_tmp(6)-nd_tmp(5)*10
      nd_tmp(5)=nd_tmp(5)-nd_tmp(4)*10
      nd_tmp(4)=nd_tmp(4)-nd_tmp(3)*10
      nd_tmp(3)=nd_tmp(3)-nd_tmp(2)*10
      nd_tmp(2)=nd_tmp(2)-nd_tmp(1)*10
      asfo=vathau*nd_tmp(1)+valonmlun*(nd_tmp(2)-5)+           &
           valonmsol*(nd_tmp(3)-5)+valonperl*(nd_tmp(4)-5)     &
           +valonascl*(nd_tmp(5)-5)+valonpers*(nd_tmp(6)-5)
      frequence(i)=max(mod(asfo/24.0,360.0),1e-10)
      omega(i)=2*pi*frequence(i)*(pi/180.0/3600.0)
    end do
  end subroutine tide_data


  subroutine set_wave(iw,name,dood,freq,equi)
!
!  Store one entry of the wave catalog (used by tide_data).
!
    implicit none
    integer, intent(in)          :: iw,dood
    character(len=*), intent(in) :: name
    real, intent(in)             :: freq,equi

    nommar(iw)    = name
    doodson(iw)   = dood
    frequence(iw) = freq
    equitide(iw)  = equi
  end subroutine set_wave


  subroutine mas_init_ssh(Istr,Iend,Jstr,Jend,IstrR,IendR,       &
                          JstrR,JendR,Ntides,k_of,nn_of,         &
                          SSH_Tamp,SSH_Tphase,mrssh,mgssh,mhssh)
!
!  Convert the standard SSH_Tamp/SSH_Tphase harmonic amplitude/phase
!  into SHOM reduced amplitude/phase tables, at every grid point.
!
    implicit none
    integer, intent(in) :: Istr,Iend,Jstr,Jend,IstrR,IendR,       &
                           JstrR,JendR,Ntides
    integer, intent(in) :: k_of(:),nn_of(:)
    real, intent(in)    :: SSH_Tamp(GLOBAL_2D_ARRAY,Ntides)
    real, intent(in)    :: SSH_Tphase(GLOBAL_2D_ARRAY,Ntides)
    real, intent(inout) :: mrssh(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mgssh(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mhssh(130,GLOBAL_2D_ARRAY)
    integer :: i,j,ht,k,nn,ir0

    mhssh(129,:,:) = -999.0
    do ht=1,Ntides
      if (l_presence(nu(ht))) then
        ir0=rr0(nu(ht))
        k =k_of(ht)
        nn=nn_of(ht)
        do j=JstrR,JendR
          do i=IstrR,IendR
            mrssh(nn,k,i,j)=SSH_Tamp(i,j,ht)*100.0
            mgssh(nn,k,i,j)=SSH_Tphase(i,j,ht)
            if (ir0 < 0) mgssh(nn,k,i,j)=mgssh(nn,k,i,j)+180.0
          enddo
        enddo
      endif
    enddo
  end subroutine mas_init_ssh


  subroutine mas_init_pot(Istr,Iend,Jstr,Jend,IstrR,IendR,       &
                          JstrR,JendR,Ntides,k_of,nn_of,         &
                          POT_Tamp,POT_Tphase,mrpot,mgpot,mhpot)
!
    implicit none
    integer, intent(in) :: Istr,Iend,Jstr,Jend,IstrR,IendR,       &
                           JstrR,JendR,Ntides
    integer, intent(in) :: k_of(:),nn_of(:)
    real, intent(in)    :: POT_Tamp(GLOBAL_2D_ARRAY,Ntides)
    real, intent(in)    :: POT_Tphase(GLOBAL_2D_ARRAY,Ntides)
    real, intent(inout) :: mrpot(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mgpot(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mhpot(130,GLOBAL_2D_ARRAY)

    call mas_init_ssh(Istr,Iend,Jstr,Jend,IstrR,IendR,            &
                      JstrR,JendR,Ntides,k_of,nn_of,              &
                      POT_Tamp,POT_Tphase,mrpot,mgpot,mhpot)
  end subroutine mas_init_pot


# if defined UV_TIDES
  subroutine mas_init_uv(Istr,Iend,Jstr,Jend,IstrR,IendR,        &
                         JstrR,JendR,Ntides,k_of,nn_of,          &
                         UV_Tmajor,UV_Tminor,UV_Tangle,UV_Tphase, &
                         angler,mru,mgu,mrv,mgv,                 &
                         mhu,mhv)
!
!  Decompose the U/V tidal ellipses into complex U/V harmonics
!  (ep2ap), then convert their rho-to-u/v averages into SHOM
!  reduced amplitude/phase tables, at every u- and v-point.
!
    implicit none
    integer, intent(in) :: Istr,Iend,Jstr,Jend,IstrR,IendR,       &
                           JstrR,JendR,Ntides
    integer, intent(in) :: k_of(:),nn_of(:)
    real, intent(in)    :: UV_Tmajor(GLOBAL_2D_ARRAY,Ntides)
    real, intent(in)    :: UV_Tminor(GLOBAL_2D_ARRAY,Ntides)
    real, intent(in)    :: UV_Tangle(GLOBAL_2D_ARRAY,Ntides)
    real, intent(in)    :: UV_Tphase(GLOBAL_2D_ARRAY,Ntides)
    real, intent(in)    :: angler(GLOBAL_2D_ARRAY)
    real, intent(inout) :: mru(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mgu(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mrv(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mgv(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mhu(130,GLOBAL_2D_ARRAY)
    real, intent(inout) :: mhv(130,GLOBAL_2D_ARRAY)
    complex :: u_cmplx(GLOBAL_2D_ARRAY)
    complex :: v_cmplx(GLOBAL_2D_ARRAY)
    integer :: i,j,ht,k,nn,ir0
    real    :: angle
    complex :: U_interp,V_interp
!  Shadows the module's own (deliberately truncated, to match fft's
!  reference behavior) private rad2deg: this atan2-to-degrees step
!  used to run in the caller against scalars.h's full-precision
!  rad2deg before this code moved into the module, and needs that
!  same precision to avoid a last-digit drift in mgu/mgv.
    real(8), parameter :: rad2deg=180.0d0/3.14159265358979323846d0

    mhu(129,Istr:,:) = -999.0
    mhv(129,:,Jstr:) = -999.0

    do ht=1,Ntides
      if (l_presence(nu(ht))) then
        ir0=rr0(nu(ht))
        k =k_of(ht)
        nn=nn_of(ht)

        do j=Jstr-1,JendR
          do i=Istr-1,IendR
            angle=UV_Tangle(i,j,ht)
#  ifdef CURVGRID
            angle=angle-angler(i,j)*rad2deg
#  endif
            call ep2ap(UV_Tmajor(i,j,ht),UV_Tminor(i,j,ht),       &
                       angle,UV_Tphase(i,j,ht),                   &
                       u_cmplx(i,j),v_cmplx(i,j))
          enddo
        enddo
        do j=JstrR,JendR
          do i=Istr,IendR
            U_interp = 0.5*(u_cmplx(i,j)+u_cmplx(i-1,j))
            mru(nn,k,i,j) = abs(U_interp)*100.0
            mgu(nn,k,i,j) = atan2(aimag(U_interp),real(U_interp))  &
                            *rad2deg
            if (mgu(nn,k,i,j) < 0.0) mgu(nn,k,i,j)=mgu(nn,k,i,j)+360.0
            if (ir0 < 0) mgu(nn,k,i,j)=mgu(nn,k,i,j)+180.0
          enddo
        enddo
        do j=Jstr,JendR
          do i=IstrR,IendR
            V_interp = 0.5*(v_cmplx(i,j)+v_cmplx(i,j-1))
            mrv(nn,k,i,j) = abs(V_interp)*100.0
            mgv(nn,k,i,j) = atan2(aimag(V_interp),real(V_interp))  &
                            *rad2deg
            if (mgv(nn,k,i,j) < 0.0) mgv(nn,k,i,j)=mgv(nn,k,i,j)+360.0
            if (ir0 < 0) mgv(nn,k,i,j)=mgv(nn,k,i,j)+180.0
          enddo
        enddo
      endif
    enddo
  end subroutine mas_init_uv
# endif


  subroutine mas_update_calendar(time)
!
!  Decomposes time into a calendar date and flags whether a new
!  day has started since the last call, publishing the result via
!  the module mas_ijour/mas_imois/mas_ian/mas_tempis/mas_is_new_day.
!  Must be called exactly once per time step, before any tile reads
!  those values (e.g. under OMP MASTER + BARRIER, ahead of the tile
!  loop that calls clm_tides): this is the only place jourmem is
!  touched, so calling it per-tile instead would race across threads.
!
    implicit none
    real, intent(in) :: time
    character(len=19) :: cdate,tool_sectodat
    integer :: iheure,iminu,isec

    cdate = tool_sectodat(time)
    call tool_decompdate(cdate,mas_ijour,mas_imois,mas_ian,       &
                         iheure,iminu,isec)
    mas_tempis = iheure*3600.0+iminu*60.0+isec
    mas_is_new_day = (mas_ijour /= jourmem)
    if (mas_is_new_day) jourmem = mas_ijour
  end subroutine mas_update_calendar


  subroutine mas_daily_update(imin,imax,jmin,jmax,                &
                              new_day,ijour,imois,ian,r,g,h)
!
!  Once per new day, rebuild the 128-heights/day reduced-tide table
!  (ma1) at every grid point, for a single (r,g,h) field family
!  (e.g. mrssh/mgssh/mhssh, or mru/mgu/mhu). Call once per active field.
!
    implicit none
    integer, intent(in) :: imin,imax,jmin,jmax
    logical, intent(in) :: new_day
    integer, intent(in) :: ijour,imois,ian
    real, intent(in)    :: r(30,8,GLOBAL_2D_ARRAY)
    real, intent(in)    :: g(30,8,GLOBAL_2D_ARRAY)
    real, intent(inout) :: h(130,GLOBAL_2D_ARRAY)
    integer :: i,j

    if (new_day) then
      do j=jmin,jmax
        do i=imin,imax
          call ma1(ijour,imois,ian,h(:,i,j),r(:,:,i,j),g(:,:,i,j))
        enddo
      enddo
    endif
  end subroutine mas_daily_update


  subroutine mas_interpolate(imin,imax,jmin,jmax,                 &
                             Istr,Iend,Jstr,Jend,tempis,h,X)
!
!  Interpolate the current 128-heights/day table to the current time
!  of day (ma7), at every grid point, for a single h/X field pair
!  (e.g. mhssh/Etide, or mhu/Utide). Call once per active field.
!
!  X is a per-tile private scratch array (Istr-2:Iend+2,Jstr-2:
!  Jend+2), unlike h which is a persistent global array: its
!  bounds are a runtime offset of the caller own tile, so Istr/
!  Iend/Jstr/Jend (not just imin/imax/jmin/jmax) must be passed in
!  to declare it with matching bounds.
!
    implicit none
    integer :: i,j
    integer, intent(in) :: imin,imax,jmin,jmax
    integer, intent(in) :: Istr,Iend,Jstr,Jend
    real, intent(in)    :: tempis
    real, intent(in)    :: h(130,GLOBAL_2D_ARRAY)
    real, intent(out)   :: X(PRIVATE_2D_SCRATCH_ARRAY)

    do j=jmin,jmax
      do i=imin,imax
        call ma7(tempis,h(:,i,j),X(i,j))
      enddo
    enddo
  end subroutine mas_interpolate


# if defined UV_TIDES
  subroutine ep2ap(semi_major, semi_minor, inclination, phase,    &
                   u_complex, v_complex)
!
! Compute U/V amplitude and phase from ellipse parameters
!
    implicit none
    real(8), intent(in)  :: semi_major, semi_minor
    real(8), intent(in)  :: inclination, phase
    complex, intent(out) :: u_complex, v_complex

    complex(8) :: wp, wm
    real(8) :: incl_rad, phase_rad

    real, parameter :: pi=3.14159265358979323846
    real, parameter :: deg2rad=pi/180.
    real, parameter :: rad2deg=180./pi

    ! Convert to radians
    incl_rad = inclination * deg2rad
    phase_rad = phase * deg2rad

    ! Decompose ellipse into counter-rotating circles
    wp = (semi_major + semi_minor) / 2.0d0 *                      &
         exp(cmplx(0.0d0,-1.D0)*(incl_rad - phase_rad))
    wm = (semi_major - semi_minor) / 2.0d0 *                      &
         exp(cmplx(0.0d0,-1.D0)*(incl_rad + phase_rad))

    ! Reconstruct U and V components
    u_complex = wp + conjg(wm)
    v_complex = cmplx(0.0d0, 1.0d0) * (wp - conjg(wm))

  end subroutine ep2ap
# endif
!******************************************************
!
! Interpolation between the 128 precomputed heights
!
  subroutine ma7(seci,hn,haut)
    implicit none
    integer               :: idk,l,k
    real                  :: d1,haut,fk,d,hh,d2
    real, dimension(130)  :: hn
    real                  :: seci

    fk=seci/60.0/11.25+2.0
    k=fk+0.5
    d=fk-k
    if (k >= 130) then
      idk=k-129
      k=k-idk
      d=d+idk
    endif
    l=k
    if (l == 130) l=129
    hh=hn(l)+hn(l)
    d1=hn(l+1)-hn(l-1)
    d2=hn(l+1)+hn(l-1)-hh

    haut=(hh+d*(d*d2+d1))/800.0
  end subroutine ma7



  subroutine ma1(jour,mois,ia,hn,rr,gg)
!
!*****  compute hn: 128 reduced tide heights per day
!
    implicit none
    integer :: l,nj,i,nd0,n,m,jour,mois,ia,kg,ng,j,k
    integer :: ll,kk
    integer, dimension(30) :: noa
    real                   :: fq,co,so
    double precision       :: fnjd
    real, dimension(30)    :: f,q,q0,r1,r2,v0
    real, dimension(31)    :: g,r
    real, dimension(130)   :: hn
    real, dimension(91), parameter :: cs=(/ &
      1.000000, 0.999848, 0.999391, 0.998630, 0.997564, 0.996195, 0.994522, &
      0.992546, 0.990268, 0.987688, 0.984808, 0.981627, 0.978148, 0.974370, &
      0.970296, 0.965926, 0.961262, 0.956305, 0.951057, 0.945519, 0.939693, &
      0.933580, 0.927184, 0.920505, 0.913545, 0.906308, 0.898794, 0.891007, &
      0.882948, 0.874620, 0.866025, 0.857167, 0.848048, 0.838671, 0.829038, &
      0.819152, 0.809017, 0.798636, 0.788011, 0.777146, 0.766044, 0.754710, &
      0.743145, 0.731354, 0.719340, 0.707107, 0.694658, 0.681998, 0.669131, &
      0.656059, 0.642788, 0.629320, 0.615661, 0.601815, 0.587785, 0.573576, &
      0.559193, 0.544639, 0.529919, 0.515038, 0.500000, 0.484810, 0.469472, &
      0.453991, 0.438371, 0.422618, 0.406737, 0.390731, 0.374607, 0.358368, &
      0.342020, 0.325568, 0.309017, 0.292372, 0.275637, 0.258819, 0.241922, &
      0.224951, 0.207912, 0.190809, 0.173648, 0.156434, 0.139173, 0.121869, &
      0.104528, 0.087156, 0.069756, 0.052336, 0.034899, 0.017452, 0.000000 /)
!  c(k) = cos((k-1) deg), k=1..362, built from the 0..90 deg table cs
    real, dimension(362), parameter :: c=(/ cs(1:91),     & ! c(1:91)    :   0..90  deg
                                            -cs(90:1:-1), & ! c(92:181)  :  91..180 deg
                                            -cs(2:90),    & ! c(182:270) : 181..269 deg
                                            0.0,          & ! c(271)     : 270      deg
                                            cs(90:1:-1),  & ! c(272:361) : 271..360 deg
                                            cs(2) /)        ! c(362)     : 361      deg
    real, dimension(30,8)  :: rr,gg
    real, dimension(11,3)  :: x,y

    l=1

    x(:,:)=0.0
    y(:,:)=0.0

    nj=julien2(ia,jour,mois)-2415021
    fnjd=real(nj+l,8)-dble(1.4921875)
    nj=mod(nj,7)+1
    if (nj <= 0) nj=nj+7
    nd0=-1
    n=0
    m=1

    do kg=1,8
      n=num(kg)
      if (n /= 0) then
        do i=1,n
          r(i)=rr(i,kg)
          g(i)=gg(i,kg)
          noa(i)=nomd(i,kg)
          r1(i)=fr1(i,kg)
          r2(i)=fr2(i,kg)
        end do
        ! mas3 computes nodal correction
        call ma3(f,v0,q,n,r,g,fnjd,ng,r1,r2,noa)

        do i=1,n
          r(i)=r(i)*f(i)
          q0(i)=v0(i)-g(i)+360.0
          q(i)=mod(q(i),360.0)
          q0(i)=mod(q0(i),360.0)
          if (q0(i) < 0.0) q0(i)=q0(i)+360.0
        end do
        do i=l,3
          do j=1,n
            fq=q0(j)
            k=fq+1
            fq=fq-k+1
            co=((c(k+1)-c(k))*fq+c(k))*r(j)
            k=k-90
            if (k <= 0) k=k+360
            so=((c(k+1)-c(k))*fq+c(k))*r(j)
            x(ng,i)=x(ng,i)+co
            y(ng,i)=y(ng,i)+so
            q0(j)=q0(j)+q(j)
            if (q0(j) > 360.0) q0(j)=q0(j)-360.0
          end do
        end do
      end if
    end do
    call ma4(hn,l,x,y)
  end subroutine ma1


  subroutine ma4(hn,l,hm1,hm2)
!  Compute the reduced heights h for days day-1, day and day+1
!  and call the subroutine computing the 128 heights per day (ma5)
    implicit none
    integer                 :: k,l,i
    real                    :: fnm
    real, dimension(130)    :: hn
    real, dimension(11,3)   :: hm1,hm2
    real, dimension(128,3)  :: h
    complex, dimension(128) :: x

    do k=l,3
      h(:,k)=0.0
      x(:)=0.0
      x(1)=hm1(1,k)
      fnm=x(1)*2.0
      hm2(1,k)=0.0
      hm1(1,k)=0.0
      do i=2,11
        x(i)=cmplx(hm1(i,k),hm2(i,k))
        hm1(i,k)=0.0
        hm2(i,k)=0.0
        x(130-i)=conjg(x(i))
      end do
      call fft(x)
      h(:,k)=real(x(:))
    end do
    call ma5(h,hn)
  end subroutine ma4


  subroutine ma5(h,x)
!
!  Compute the 128 heights by interpolation between the reduced heights
!
    implicit none
    integer                :: k
    real                   :: b,hh
    real, dimension(128,3) :: h
    real, dimension(130)   :: x

    b=-0.5078125

    x(1)=x(129)
    x(2)=x(130)

    do k=1,128
      hh=h(k,2)+h(k,2)
      b=b+0.0078125
      x(k+2)=((h(k,1)+h(k,3)-hh)*b+h(k,3)-h(k,1))*b+hh
    end do

    if (x(1) == -valmanq) then
      x(2)=3.0*(x(3)-x(4))+x(5)
      x(1)=6.0*x(3)-8.0*x(4)+3.0*x(5)
    endif

  end subroutine ma5

  subroutine fft(a)
!
!  Fast Fourier transform
!
    implicit none
    integer                 :: m,nv2,nm1,n,i,j,ip,l,k,le,le1
    real                    :: fij
    complex, dimension(128) :: a
    complex                 :: u,w,t
    m=7
    nv2=64
    nm1=127
    n=128
    j=1
    do i=1,nm1
      if (i < j) then
        t=a(j)
        a(j)=a(i)
        a(i)=t
      end if
      k=nv2
      do while (k < j)
        j=j-k
        k=k/2
      end do
      j=j+k
    end do

    do l=1,m
      le=2**l
      le1=le/2
      u=cmplx(1.0,0.0)
      fij=pi/le1
      w=exp(cmplx(0.0,fij))
      do j=1,le1
        do i=j,n,le
          ip=i+le1
          t=a(ip)*u
          a(ip)=a(i)-t
          a(i)=a(i)+t
        end do
        u=u*w
      end do
    end do
  end subroutine fft

  subroutine ma3(f,v0,q,n,r,g,fnjd,ng,r1,r2,noa)
!
!  Compute the nodal factors (amplitudes and phases)
!
    implicit none
    integer                :: i,ng,n
    integer, dimension(6)  :: nd
    integer, dimension(30) :: noa
    real                   :: fa1,fb1,fa2,fb2,a1,a2
    double precision       :: fnjd
    real, dimension(30)    :: q,v0,f,r,g,r1,r2
    double precision       :: dq
    complex                :: fc

    do i=1,n
      nd(6)=noa(i)
      nd(5)=noa(i)/10
      nd(4)=noa(i)/100
      nd(3)=noa(i)/1000
      nd(2)=noa(i)/10000
      nd(1)=noa(i)/100000
      nd(6)=nd(6)-nd(5)*10
      nd(5)=nd(5)-nd(4)*10
      nd(4)=nd(4)-nd(3)*10
      nd(3)=nd(3)-nd(2)*10
      nd(2)=nd(2)-nd(1)*10
      ng=nd(1)+1
      if (nd(1) == 0) r(i)=r(i)*2.0
      call masfo(nd,dq)
      q(i)=dq
      g(i)=mod(g(i),360.0)
      if (g(i) < 0.0) g(i)=g(i)+360.0
      v0(i)=asfov(nd)+mod(real(dq*fnjd,8),dble(360.))
      fc=1.0
      if (r1(i) /= 0.0 .or. r2(i) /= 0.0) then
        fb2=9.242202e-4
        fb1=-fb2
        fa2=1.760045
        fa1=4.523139
        if (noa(i) == 275555) then
          fb1=2*fb2
          fa1=2*fa2
        endif
        a1=fa1+fb1*fnjd
        a2=fa2+fb2*fnjd
        fc=1.0+r1(i)*cmplx(cos(a1),sin(a1))+                      &
           r2(i)*cmplx(cos(a2),sin(a2))
      end if

      v0(i)=v0(i)+atan2(aimag(fc),real(fc))*57.29578
      f(i)=abs(fc)
      if (noa(i) == 355555 .or. noa(i) == 382555 .or.            &
          noa(i) == 164555) v0(i)=v0(i)-90.0
    end do
  end subroutine ma3

  subroutine masfo(nd,asfo)
!
! Compute the angular speeds of the constituents
!
    implicit none
    integer, dimension(6) :: nd
    double precision      :: asfo

    asfo =(dble(360.)-dble(12.19074939))*(nd(1))                  &
          +dble(13.17639673)*(nd(2)-5)+dble(0.98564734)*(nd(3)-5) &
          +dble(0.11140408)*(nd(4)-5)+dble(0.05295392)*(nd(5)-5)  &
          +dble(0.00004707)*(nd(6)-5)
  end subroutine masfo

  function asfov(nd)
!
! Compute the astronomical arguments
!
    implicit none
    integer, dimension(6) :: nd
    real                  :: fov,asfov

    fov=280.1895*(nd(1)+nd(3)-5)+                                 &
        277.0248*(nd(2)-nd(1)-5)+334.3853*(nd(4)-5)               &
        +100.8432*(nd(5)-5)+281.2209*(nd(6)-5)+(mod(nd(1),2))*90.0
    asfov=(mod(fov,360.0))
  end function asfov

  function julien2(ia,jou,moi)
!
! Compute the number of days elapsed since 1 January 1900
!
    implicit none
    integer :: i,julien2,ia,jou,moi,ibs,jour,mois,iy,m
    integer, dimension(12,2), parameter :: jo=                    &
      reshape((/0,31,59,90,120,151,181,212,243,273,304,334,       &
                0,31,60,91,121,152,182,213,244,274,305,335/),(/12,2/))
    real :: a,b

    jour=jou
    mois=moi
    if (mois == 1) then
      ibs=1
      if (ia >= 1582 .and. (ia /= 1582 .or. jou >= 277)) then
        ibs=mod(ia,4)+2
        if (ibs /= 2) ibs=1
      end if
      do i=1,12
        mois=12-i+1
        jour=jou-jo(mois,ibs)
        if (jour > 0) exit
      end do
    end if

    a=ia+mois/100.0+jour/10000.0
    b=0.0
    iy=ia
    m=mois

    if (mois <= 2) then
      iy=iy-1
      m=m+12
    endif

    if (a >= 1582.1015) then
      b=2-int(iy/100)+int(iy/100)/4
    endif

    julien2=int(365.25*iy)+int(30.6001*(m+1))+jour+1720995+b
  end function julien2
#endif
end module mod_tides_mas

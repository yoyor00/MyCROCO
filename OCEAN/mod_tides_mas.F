!============================================================================================
! Reference : Coastal Tides book by Simon (SHOM)
! https://diffusion.shom.fr/coastal-tides-version-anglaise-de-la-maree-oceanique-cotiere.html
!=============================================================================================
#include "cppdefs.h"
      MODULE MOD_TIDES_MAS
#ifdef TIDES_MAS
#include "param.h"
      INTEGER,PARAMETER :: nmarmax=144
      INTEGER,DIMENSION(nmarmax)                  :: nu
      INTEGER,DIMENSION(30,8)                     :: nomd
      REAL,DIMENSION(30,8)              :: fr1,fr2
      CHARACTER(len=100),DIMENSION(nmarmax) :: nommar
      LOGICAL,DIMENSION(nmarmax)            :: l_presence
      INTEGER,DIMENSION(nmarmax)            :: doodson
      REAL,DIMENSION(nmarmax)     :: omega
      REAL,DIMENSION(nmarmax)     :: frequence,equitide

      INTEGER :: inivm,isa,issa,imsm,imm,imsf,imf,imfm,i2q1,isigma1
      INTEGER :: iq1,irho1,io1
      INTEGER :: ims1,imp1,ia19,im1,ikhi1,ipi1,i2ns2,i3m2s2,ioq2mnk2,i33,imns2
      INTEGER :: ip1,is1,ik1,ipsi1,iphi1,itheta1,ij1,iso1,ioo1,ikq1,i2mn2s2
      INTEGER :: imnus2,i2mk2,i2n22nm2,imu2ms2,i39,in2,inu2,iop2msk2
      INTEGER :: i43,imsk2,im2,imks2,imks22,ilambda2,i2mk3,im3,iso3
      INTEGER :: ims3,imk3,ia87
      INTEGER :: il22mn2,inkm2,it2,is2,ir22,ik2,imsn2,ikj2mkn2,i2sm2
      INTEGER :: iskm2,imq3
      INTEGER :: isp3,is3,isk3,ik3,i2nms4,i2mmus4,i2mns4,i2mnus4,i3mk4,in4
      INTEGER :: i3ms4,imn4,imnu4,i2msk4,im4,i2msn4,i2mkn4,is4,isk4,i3mnk6
      INTEGER :: i2mks4,isn4,i3mn4ml4,ink4,i2smk4,imt4,ims4,imk4,i2snm4
      INTEGER :: i3mns6,i3mnus6,i4mk6,i2nm6,i4ms6,i2mn6,i2mnu6,i3msk6,im6
      INTEGER :: i3mks6,imsn6,i4mn2ml6,imnk6,i2mt6,i3mnu8,i2msk8,im8,i4mks8
      INTEGER :: i2ms6,i2mk6,i2sn6,i3msn6,i3mkn6,i2sm6,imsk6,i2mns8,i2mn8,i3mn8
      INTEGER :: i2msn8,i3ml8,i2mnk8,i3ms8,i3mk8,i2smn8,i4msn8
      INTEGER :: imsnk8,i4mnk8,i2ms8
      INTEGER :: i4mn10,im10,i3msn10,i4ms10,i4mk10,i5msn10,i2msnk10,i3m2s10
      INTEGER :: i2msk82,i5mns10,i3m2n10

      INTEGER,DIMENSION(0:10),PARAMETER :: esp_corres=
     &                         (/1,2,3,4,5,0,6,0,7,0,8/)

!tidepredic
      integer                   :: jourmem=-1
!  Calendar state, refreshed once per step by mas_update_calendar
!  (called under OMP MASTER + BARRIER, before the tile loop): tiles
!  only ever read these, so no per-tile mutation races on jourmem.
      integer                   :: mas_ijour=-1,mas_imois=-1,
     &                              mas_ian=-1
      real                       :: mas_tempis=0.0
      logical                    :: mas_is_new_day=.false.
      integer,dimension(8),save :: num
      integer,dimension(144),PARAMETER :: ide=(/
     &   55555,  56555,  57555,  63655,  65455,  73555,  75555, 125755, 
     &  127555, 135655, 137455, 145555, 146555, 147555, 153655, 155655, 
     &  157455, 162556, 163555, 164555, 165555, 166554, 167555, 173655,
     &  175455, 183555, 185555, 195455, 209655, 217755, 219555, 225655,
     &  225855, 227655, 229455, 235555, 235755, 237555, 238554, 245655,
     &  247455, 253555, 253755, 254556, 255555, 256554, 257555, 263655,
     &  265455, 265655, 272556, 273555, 274554, 275555, 283455, 285455,
     &  291555, 293555, 335655, 345555, 355555, 363555, 364555, 365555,
     &  375555, 381555, 382555, 383555, 385555, 417755, 419555, 427655,
     &  429455, 435555, 435755, 437555, 445655, 447455, 453555, 455555,
     &  457555, 463655, 465455, 465655, 471555, 472556, 473555, 475555,
     &  481655, 483455, 485455, 491555, 493555, 625655, 627655, 629455,
     &  635555, 635755, 637555, 645655, 647455, 653555, 655555, 657555,
     &  663655, 665455, 665655, 672556, 673555, 675555, 681655, 683455,
     &  685455, 691555, 693555, 827655, 835755, 845655, 847455, 853555,
     &  855555, 857555, 863655, 865455, 865655, 873555, 875555, 881655,
     &  883455, 883655, 885455, 891555, 893555,1027655,1035755,1045655,
     &  1055555,1063655,1073555,1075555,1083455,1083655,1091555,85455/)
      integer*4,dimension(144)    ::   rr1=(/
     & 0,0,0,-72,-66,72,0,189,189,189,188,189,0,0,227,0,0,0,0,0,  20,
     & 0,0,0,
     & 0,0,0,0,-112,-75,-112,279,0,-75,-75,279,-37,-37,0,-37,-37,
     & 260,0,0
     & ,-37,0,-37,-37,-37,-37,0,2,0,0,-37,-37,0,0,152, 61,-56,189,
     & 0,0,0,0
     & ,0, 20, 20,-112,-112,-112,-112,186,-75,-112,-75,-75,223,-75,
     & -75,-37,
     & -112,-37,261,-37,-37,-37,-37,-75,-75,0,0,
     & 148,-149,-149,148,-112,-149,-112,-112,186,-112,-112,-75,-149,-75,
     & -75,-75,-75,-37,-112,-112,-37,-37,
     & -112,-149,-149,-149,148,-149,-149,-112,-149,-112,-112,-112,-75,
     & -149,-75,-149,-75,-75,
     & -224,-186,-186,-186,-149,-149,-149,-186,-112,-112,41/)
      integer*4,dimension(144)    :: rr2=(/
     &  0,0,0,-65,-65,-64,414,
     &  0,0,0,0,0,0,218,0,-200,-219,0,0,0,-135,0,0,0,-198,-195,-640,
     &  -640,
     &  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &  -477,298,0,0,260,0,0,0,298,-37,260,-37,260,
     &  0,20,0,0,0,135,439,0,0,-135,-433,
     &  0,0,0,0,0,0,0,0,0,0,0,
     &  298,0,-37,298,0,0,0,298,-37,-37,261,0,298,
     &  0,0,0,0,0,0,0,0,0,0,
     &  297,0,-37,297,0,0,297,0,-37,260,0,297,
     &  0,0,0,0,0,0,297,0,0,297,0,297,0,-37,297,260,0,297,
     &  0,0,0,0,0,0,297,-37,297,0,0/)

      integer*4,dimension(144)    ::  rr0=(/
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,
     &  -1000,-1000,-1000,-1000,1000,1000,1000,
     &  -1000,-1000,-1000,-1000,-1000,-1000,-1000,-1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  -1000,-1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,
     &  -1000,-1000,1000,1000,-1000,-1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,
     &  -1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000,1000,1000,1000,1000,1000,1000,1000,1000,
     &  1000/)

      integer,parameter :: rsh=4
      REAL,PARAMETER    :: pi=3.141592654
      REAL,PARAMETER    :: rad2deg=180.0/pi
      REAL(kind=rsh),PARAMETER    :: valmanq=999.
!
! pi/rad2deg are private: callers already get identically-named
! PARAMETERs from scalars.h, and USE-associating both would conflict.
      PRIVATE :: pi,rad2deg

      CONTAINS
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

        ! Trouve les ondes presente et stocke les infos
        num(:) = 0
        nomd(:,:) = 0
        fr1(:,:) = 0.0
        fr2(:,:) = 0.0

        DO ht = 1,Ntides
          IF (l_presence(nu(ht))) THEN
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
          ENDIF
        ENDDO
        end subroutine mas_init_waves


        subroutine tide_data
!
!  Populate the module-level wave catalog (nommar/doodson/frequence/
!  equitide) and each named species index (inivm, isa, ..., i3m2n10)
!  for the 144 MAS harmonics, then derive frequence/omega for all of
!  them. Call once at initialization, before mas_init_waves.
!
      implicit none

!
      !
      integer              :: i
      integer,dimension(6) :: nd_tmp
      real   :: asfo,vathau,valonmlun,valonmsol,valonperl,valonascl
      real   :: valonpers
!
!************************************************************
!     Executable statements
!************************************************************
!
!      niveau moyen
!------------------
!
      inivm=1
      nommar(inivm) ='niveau_moyen'
      doodson(inivm) =  55555
      frequence(inivm) = 0.00000
      equitide(inivm) = 0.0
!
!      sa
!------------------
!
      isa=2
      nommar(isa) ='sa'
      doodson(isa) =  56555
      frequence(isa) = 0.041068639
      equitide(isa) = 0.0
!
!      ssa
!------------------
!
      issa=3
      nommar(issa) ='ssa'
      doodson(issa) =  57555
      frequence(issa) = 0.082137278
      equitide(issa) = 0.0
!
!      msm
!------------------
!
      imsm=4
      nommar(imsm) ='msm'
      doodson(imsm) =  63655
      frequence(imsm) = 0.471521089
      equitide(imsm) = 0.0
!
!      mm
!------------------
!
      imm=5
      nommar(imm) ='mm'
      doodson(imm) =  65455
      frequence(imm) = 0.544374694
      equitide(imm) = 0.0
!
!      msf
!------------------
!
      imsf=6
      nommar(imsf) ='msf'
      doodson(imsf) =  73555
      frequence(imsf) = 1.015895783
      equitide(imsf) = 0.0
!
!      mf
!------------------
!
      imf=7
      nommar(imf) ='mf'
      doodson(imf) =  75555
      frequence(imf) = 1.098033061
      equitide(imf) = 0.0
!
!      mfm
!------------------
!
      imfm=144
      nommar(imfm) ='mfm'
      doodson(imfm) =  85455
      frequence(imfm) = 1.62240772
      equitide(imfm) = 0.0
!
!      2Q1
!------------------
!
      i2q1=8
      nommar(i2q1) ='2Q1'
      doodson(i2q1) = 125755
      frequence(i2q1) = 12.854286191
      equitide(i2q1) = 0.0
!
!      sigma1
!------------------
!
      isigma1=9
      nommar(isigma1) ='sigma1'
      doodson(isigma1) = 127555
      frequence(isigma1) = 12.927139796
      equitide(isigma1) = 0.0
!
!      Q1
!------------------
!
      iq1=10
      nommar(iq1) ='Q1'
      doodson(iq1) = 135655
      frequence(iq1) = 13.398660885
      equitide(iq1) = 0.019256
!
!   rho1
!------------------
!
      irho1=11
      nommar(irho1) ='rho1'
      doodson(irho1) = 137455
      frequence(irho1) = 13.471514490
      equitide(irho1) = 0.0
!
!   O1
!------------------
!
      io1=12
      nommar(io1) ='O1'
      doodson(io1) = 145555
      frequence(io1) = 13.943035578
      equitide(io1) =0.100574
!
!   MS1
!------------------
!
      ims1=13
      nommar(ims1) ='MS1'
      doodson(ims1) = 146555
      frequence(ims1) = 13.984104217
      equitide(ims1) = 0.0
!
!   MP1
!------------------
!
      imp1=14
      nommar(imp1) ='MP1'
      doodson(imp1) = 147555
      frequence(imp1) = 14.025172857
      equitide(imp1) = 0.0
!
!   a19
!------------------
!
      ia19=15
      nommar(ia19) ='a19'
      doodson(ia19) = 153655
      frequence(ia19) = 14.414556667
      equitide(ia19) = 0.0
!
!   M1
!------------------
!
      im1=16
      nommar(im1) ='M1'
      doodson(im1) = 155655
      frequence(im1) = 14.496693945
      equitide(im1) = 0.0
!
!   khi1
!------------------
!
      ikhi1=17
      nommar(ikhi1) ='khi1'
      doodson(ikhi1) = 157455
      frequence(ikhi1) = 14.569547550
      equitide(ikhi1) = 0.0
!
!   pi1
!------------------
!
      ipi1=18
      nommar(ipi1) ='pi1'
      doodson(ipi1) = 162556
      frequence(ipi1) = 14.917864683
      equitide(ipi1) = 0.0
!
!   p1
!------------------
!
      ip1=19
      nommar(ip1) ='p1'
      doodson(ip1) = 163555
      frequence(ip1) = 14.958931361
      equitide(ip1) =0.046843
!
!   s1
!------------------
!
      is1=20
      nommar(is1) ='s1'
      doodson(is1) = 164555
      frequence(is1) = 15.0
      equitide(is1) = 0.0
!
!   k1
!------------------
!
      ik1=21
      nommar(ik1) ='k1'
      doodson(ik1) = 165555
      frequence(ik1) = 15.041068639
      equitide(ik1) = 0.141565
!
!   psi1
!------------------
!
      ipsi1=22
      nommar(ipsi1) ='psi1'
      doodson(ipsi1) = 166554
      frequence(ipsi1) = 15.082135317
      equitide(ipsi1) = 0.0
!
!   phi1
!------------------
!
      iphi1=23
      nommar(iphi1) ='phi1'
      doodson(iphi1) = 167555
      frequence(iphi1) = 15.123205917
      equitide(iphi1) = 0.0
!
!   tetha1
!------------------
!
      itheta1=24
      nommar(itheta1) ='theta1'
      doodson(itheta1) = 173655
      frequence(itheta1) = 15.512589728
      equitide(itheta1) = 0.0
!
!   J1
!------------------
!
      ij1=25
      nommar(ij1) ='J1'
      doodson(ij1) = 175455
      frequence(ij1) = 15.585443333
      equitide(ij1) = 0.0
!
!  SO1
!------------------
!
      iso1=26
      nommar(iso1) ='SO1'
      doodson(iso1) = 183555
      frequence(iso1) = 16.056964422
      equitide(iso1) = 0.0
!
!  OO1
!------------------
!
      ioo1=27
      nommar(ioo1) ='OO1'
      doodson(ioo1) = 185555
      frequence(ioo1) = 16.1391017
      equitide(ioo1) = 0.0
!
!  KQ1
!------------------
!
      ikq1=28
      nommar(ikq1) ='KQ1'
      doodson(ikq1) = 195455
      frequence(ikq1) = 16.683476394
      equitide(ikq1) = 0.0
!
!  2MN2S2
!------------------
!
      i2mn2s2=29
      nommar(i2mn2s2) ='2MN2S2'
      doodson(i2mn2s2) = 209655
      frequence(i2mn2s2) = 26.407937959
      equitide(i2mn2s2) = 0.0
!
!  2NS2
!------------------
!
      i2ns2=30
      nommar(i2ns2) ='2NS2'
      doodson(i2ns2) = 217755
      frequence(i2ns2) = 26.879459048
      equitide(i2ns2) = 0.0
!
!  3M2S2
!------------------
!
      i3m2s2=31
      nommar(i3m2s2) ='3M2S2'
      doodson(i3m2s2) = 219555
      frequence(i3m2s2) = 26.952312652
      equitide(i3m2s2) = 0.0
!
!  OQ2MNK2
!------------------
!
      ioq2mnk2=32
      nommar(ioq2mnk2) ='OQ2MNK2'
      doodson(ioq2mnk2) = 225655
      frequence(ioq2mnk2) = 27.341696463
      equitide(ioq2mnk2) = 0.0
!
!  33
!------------------
!
      i33=33
      nommar(i33) ='33'
      doodson(i33) = 225855
      frequence(i33) = 27.350980136
      equitide(i33) = 0.0
!
!  MNS2
!------------------
!
      imns2=34
      nommar(imns2) ='MNS2'
      doodson(imns2) = 227655
      frequence(imns2) = 27.423833741
      equitide(imns2) = 0.0
!
!  MNUS2
!------------------
!
      imnus2=35
      nommar(imnus2) ='MNUS2'
      doodson(imnus2) = 229455
      frequence(imnus2) = 27.496687346
      equitide(imnus2) = 0.0
!
!  2MK2
!------------------
!
      i2mk2=36
      nommar(i2mk2) ='2MK2'
      doodson(i2mk2) = 235555
      frequence(i2mk2) = 27.886071157
      equitide(i2mk2) = 0.0
!
!  2N22NM2
!------------------
!
      i2n22nm2=37
      nommar(i2n22nm2) ='2N22NM2'
      doodson(i2n22nm2) = 235755
      frequence(i2n22nm2) = 27.89535483
      equitide(i2n22nm2) = 0.0
!
!  MU22MS2
!------------------
!
      imu2ms2=38
      nommar(imu2ms2) ='MU22MS2'
      doodson(imu2ms2) = 237555
      frequence(imu2ms2) = 27.968208435
      equitide(imu2ms2) = 0.0
!
!  39
!------------------
!
      i39=39
      nommar(i39) ='39'
      doodson(i39) = 238554
      frequence(i39) = 28.009275113
      equitide(i39) = 0.0
!
!  N2
!------------------
!
      in2=40
      nommar(in2) ='N2'
      doodson(in2) = 245655
      frequence(in2) = 28.439729524
      equitide(in2) = 0.046398
!
!  nu2
!------------------
!
      inu2=41
      nommar(inu2) ='nu2'
      doodson(inu2) = 247455
      frequence(inu2) = 28.512583129
      equitide(inu2) = 0.0
!
!  OP2MSK2
!------------------
!
      iop2msk2=42
      nommar(iop2msk2) ='OP2MSK2'
      doodson(iop2msk2) = 253555
      frequence(iop2msk2) = 28.901966939
      equitide(iop2msk2) = 0.0
!
!  43
!------------------
!
      i43=43
      nommar(i43) ='43'
      doodson(i43) = 253755
      frequence(i43) = 28.911250612
      equitide(i43) = 0.0
!
!  M(SK)2
!------------------
!
      imsk2=44
      nommar(imsk2) ='M_SK_2'
      doodson(imsk2) = 254556
      frequence(imsk2) = 28.94303754
      equitide(imsk2) = 0.0
!
!  M2
!------------------
!
      im2=45
      nommar(im2) ='M2'
      doodson(im2) = 255555
      frequence(im2) = 28.984104218
      equitide(im2) = 0.24230
!
!  M(KS)2
!------------------
!
      imks2=46
      nommar(imks2) ='M_KS_2'
      doodson(imks2) = 256554
      frequence(imks2) = 29.025170895
      equitide(imks2) = 0.0
!
!  MKS2
!------------------
!
      imks22=47
      nommar(imks22) ='MKS2'
      doodson(imks22) = 257555
      frequence(imks22) = 29.066241496
      equitide(imks22) = 0.0
!
!  lambda2
!------------------
!
      ilambda2=48
      nommar(ilambda2) ='lambda2'
      doodson(ilambda2) = 263655
      frequence(ilambda2) = 29.455625306
      equitide(ilambda2) = 0.0
!
!  L2 2MN2
!------------------
!
      il22mn2=49
      nommar(il22mn2) ='L2_2MN2'
      doodson(il22mn2) = 265455
      frequence(il22mn2) = 29.528478911
      equitide(il22mn2) = 0.0
!
!  NKM2
!------------------
!
      inkm2=50
      nommar(inkm2) ='NKM2'
      doodson(inkm2) = 265655
      frequence(inkm2) = 29.537762585
      equitide(inkm2) = 0.0
!
!  T2
!------------------
!
      it2=51
      nommar(it2) ='T2'
      doodson(it2) = 272556
      frequence(it2) = 29.958933322
      equitide(it2) = 0.0
!
!  S2
!------------------
!
      is2=52
      nommar(is2) ='S2'
      doodson(is2) = 273555
      frequence(is2) = 30.0
      equitide(is2) = 0.112841
!
!  R2
!------------------
!
      ir22=53
      nommar(ir22) ='R2'
      doodson(ir22) = 274554
      frequence(ir22) = 30.041066678
      equitide(ir22) = 0.0
!
!  K2
!------------------
!
      ik2=54
      nommar(ik2) ='K2'
      doodson(ik2) = 275555
      frequence(ik2) = 30.082137278
      equitide(ik2) =   0.030704
!
!  MSN2
!------------------
!
      imsn2=55
      nommar(imsn2) ='MSN2'
      doodson(imsn2) = 283455
      frequence(imsn2) = 30.544374694
      equitide(imsn2) = 0.0
!
!  KJ2MKN2
!------------------
!
      ikj2mkn2=56
      nommar(ikj2mkn2) ='KJ2MKN2'
      doodson(ikj2mkn2) = 285455
      frequence(ikj2mkn2) = 30.626511972
      equitide(ikj2mkn2) = 0.0
!
!  2SM2
!------------------
!
      i2sm2=57
      nommar(i2sm2) ='2SM2'
      doodson(i2sm2) = 291555
      frequence(i2sm2) = 31.015895782
      equitide(i2sm2) = 0.0
!
!  SKM2
!------------------
!
      iskm2=58
      nommar(iskm2) ='SKM2'
      doodson(iskm2) = 293555
      frequence(iskm2) = 31.098033061
      equitide(iskm2) = 0.0
!
!  MQ3
!------------------
!
      imq3=59
      nommar(imq3) ='MQ3'
      doodson(imq3) = 335655
      frequence(imq3) = 42.382765102
      equitide(imq3) = 0.0
!
!  2MK3
!------------------
!
      i2mk3=60
      nommar(i2mk3) ='2MK3'
      doodson(i2mk3) = 345555
      frequence(i2mk3) = 42.927139796
      equitide(i2mk3) = 0.0
!
!  M3
!------------------
!
      im3=61
      nommar(im3) ='M3'
      doodson(im3) = 355555
      frequence(im3) = 43.476156326
      equitide(im3) = 0.0
!
!  SO3
!------------------
!
      iso3=62
      nommar(iso3) ='SO3'
      doodson(iso3) = 363555
      frequence(iso3) = 43.943035578
      equitide(iso3) = 0.0
!
!  MS3
!------------------
!
      ims3=63
      nommar(ims3) ='MS3'
      doodson(ims3) = 364555
      frequence(ims3) = 43.984104217
      equitide(ims3) = 0.0
!
!  MK3
!------------------
!
      imk3=64
      nommar(imk3) ='MK3'
      doodson(imk3) = 365555
      frequence(imk3) = 44.025172857
      equitide(imk3) = 0.0
!
!  A87
!------------------
!
      ia87=65
      nommar(ia87) ='A87'
      doodson(ia87) = 375555
      frequence(ia87) = 44.574189387
      equitide(ia87) = 0.0
!
!  SP3
!------------------
!
      isp3=66
      nommar(isp3) ='SP3'
      doodson(isp3) = 381555
      frequence(isp3) = 44.958931361
      equitide(isp3) = 0.0
!
!  s3
!------------------
!
      is3=67
      nommar(is3) ='s3'
      doodson(is3) = 382555
      frequence(is3) = 45.0
      equitide(is3) = 0.0
!
!  sk3
!------------------
!
      isk3=68
      nommar(isk3) ='sk3'
      doodson(isk3) = 383555
      frequence(isk3) = 45.041068639
      equitide(isk3) = 0.0
!
!  k3
!------------------
!
      ik3=69
      nommar(ik3) ='k3'
      doodson(ik3) = 385555
      frequence(ik3) = 45.123205918
      equitide(ik3) = 0.0
!
!  2NMS4
!------------------
!
      i2nms4=70
      nommar(i2nms4) ='2NMS4'
      doodson(i2nms4) = 417755
      frequence(i2nms4) = 55.863563265
      equitide(i2nms4) = 0.0
!
!  2MMUS4
!------------------
!
      i2mmus4=71
      nommar(i2mmus4) ='2MMUS4'
      doodson(i2mmus4) = 419555
      frequence(i2mmus4) = 55.936416870
      equitide(i2mmus4) = 0.0
!
!  2MNS4
!------------------
!
      i2mns4=72
      nommar(i2mns4) ='2MNS4'
      doodson(i2mns4) = 427655
      frequence(i2mns4) = 56.407937959
      equitide(i2mns4) = 0.0
!
!  2MNUS4
!------------------
!
      i2mnus4=73
      nommar(i2mnus4) ='2MNUS4'
      doodson(i2mnus4) = 429455
      frequence(i2mnus4) = 56.480791564
      equitide(i2mnus4) = 0.0
!
!  3MK4
!------------------
!
      i3mk4=74
      nommar(i3mk4) ='3MK4'
      doodson(i3mk4) = 435555
      frequence(i3mk4) = 56.870175374
      equitide(i3mk4) = 0.0
!
!  N4
!------------------
!
      in4=75
      nommar(in4) ='N4'
      doodson(in4) = 435755
      frequence(in4) = 56.879459047
      equitide(in4) = 0.0
!
!  3MS4
!------------------
!
      i3ms4=76
      nommar(i3ms4) ='3MS4'
      doodson(i3ms4) = 437555
      frequence(i3ms4) = 56.952312653
      equitide(i3ms4) = 0.0
!
!  MN4
!------------------
!
      imn4=77
      nommar(imn4) ='MN4'
      doodson(imn4) = 445655
      frequence(imn4) = 57.423833741
      equitide(imn4) = 0.0
!
!  MNU4
!------------------
!
      imnu4=78
      nommar(imnu4) ='MNU4'
      doodson(imnu4) = 447455
      frequence(imnu4) = 57.496687346
      equitide(imnu4) = 0.0
!
!  2MSK4
!------------------
!
      i2msk4=79
      nommar(i2msk4) ='2MSK4'
      doodson(i2msk4) = 453555
      frequence(i2msk4) = 57.886071157
      equitide(i2msk4) = 0.0
!
!  M4
!------------------
!
      im4=80
      nommar(im4) ='M4'
      doodson(im4) = 455555
      frequence(im4) = 57.968208435
      equitide(im4) = 0.0
!
!  2MKS4
!------------------
!
      i2mks4=81
      nommar(i2mks4) ='2MKS4'
      doodson(i2mks4) = 457555
      frequence(i2mks4) = 57.968208435
      equitide(i2mks4) = 0.0
!
!  SN4
!------------------
!
      isn4=82
      nommar(isn4) ='SN4'
      doodson(isn4) = 463655
      frequence(isn4) = 57.968208435
      equitide(isn4) = 0.0
!
!  3MN4ML4
!------------------
!
      i3mn4ml4=83
      nommar(i3mn4ml4) ='3MN4ML4'
      doodson(i3mn4ml4) = 465455
      frequence(i3mn4ml4) = 57.968208435
      equitide(i3mn4ml4) = 0.0
!
!  NK4
!------------------
!
      ink4=84
      nommar(ink4) ='NK4'
      doodson(ink4) = 465655
      frequence(ink4) = 57.968208435
      equitide(ink4) = 0.0
!
!  2SMK4
!------------------
!
      i2smk4=85
      nommar(i2smk4) ='2SMK4'
      doodson(i2smk4) = 471555
      frequence(i2smk4) = 57.968208435
      equitide(i2smk4) = 0.0
!
!  MT4
!------------------
!
      imt4=86
      nommar(imt4) ='MT4'
      doodson(imt4) = 472556
      frequence(imt4) = 57.968208435
      equitide(imt4) = 0.0
!
!  MS4
!------------------
!
      ims4=87
      nommar(ims4) ='MS4'
      doodson(ims4) = 473555
      frequence(ims4) = 57.968208435
      equitide(ims4) = 0.0
!
!  MK4
!------------------
!
      imk4=88
      nommar(imk4) ='MK4'
      doodson(imk4) = 475555
      frequence(imk4) = 57.968208435
      equitide(imk4) = 0.0
!
!  2SNM4
!------------------
!
      i2snm4=89
      nommar(i2snm4) ='2SNM4'
      doodson(i2snm4) = 481655
      frequence(i2snm4) = 57.968208435
      equitide(i2snm4) = 0.0
!
!  2MSN4
!------------------
!
      i2msn4=90
      nommar(i2msn4) ='2MSN4'
      doodson(i2msn4) = 483455
      frequence(i2msn4) = 57.968208435
      equitide(i2msn4) = 0.0
!
!  2MKN4
!------------------
!
      i2mkn4=91
      nommar(i2mkn4) ='2MKN4'
      doodson(i2mkn4) = 485455
      frequence(i2mkn4) = 57.968208435
      equitide(i2mkn4) = 0.0
!
!  S4
!------------------
!
      is4=92
      nommar(is4) ='S4'
      doodson(is4) = 491555
      frequence(is4) = 57.968208435
      equitide(is4) = 0.0
!
!  SK4
!------------------
!
      isk4=93
      nommar(isk4) ='SK4'
      doodson(isk4) = 493555
      frequence(isk4) = 57.968208435
      equitide(isk4) = 0.0
!
!  3MNK6
!------------------
!
      i3mnk6=94
      nommar(i3mnk6) ='3MNK6'
      doodson(i3mnk6) = 625655
      frequence(i3mnk6) = 57.968208435
      equitide(i3mnk6) = 0.0
!
!  3MNS6
!------------------
!
      i3mns6=95
      nommar(i3mns6) ='3MNS6'
      doodson(i3mns6) = 627655
      frequence(i3mns6) = 57.968208435
      equitide(i3mns6) = 0.0
!
!  3MNUS6
!------------------
!
      i3mnus6=96
      nommar(i3mnus6) ='3MNUS6'
      doodson(i3mnus6) = 629455
      frequence(i3mnus6) = 57.968208435
      equitide(i3mnus6) = 0.0
!
!  4MK6
!------------------
!
      i4mk6=97
      nommar(i4mk6) ='4MK6'
      doodson(i4mk6) = 635555
      frequence(i4mk6) = 57.968208435
      equitide(i4mk6) = 0.0
!
!  2NM6
!------------------
!
      i2nm6=98
      nommar(i2nm6) ='2NM6'
      doodson(i2nm6) = 635755
      frequence(i2nm6) = 57.968208435
      equitide(i2nm6) = 0.0
!
!  4MS6
!------------------
!
      i4ms6=99
      nommar(i4ms6) ='4MS6'
      doodson(i4ms6) = 637555
      frequence(i4ms6) = 57.968208435
      equitide(i4ms6) = 0.0
!
!  2MN6
!------------------
!
      i2mn6=100
      nommar(i2mn6) ='2MN6'
      doodson(i2mn6) = 645655
      frequence(i2mn6) = 57.968208435
      equitide(i2mn6) = 0.0
!
!  2MNU6
!------------------
!
      i2mnu6=101
      nommar(i2mnu6) ='2MNU6'
      doodson(i2mnu6) = 647455
      frequence(i2mnu6) = 57.968208435
      equitide(i2mnu6) = 0.0
!
!  3MSK6
!------------------
!
      i3msk6=102
      nommar(i3msk6) ='3MSK6'
      doodson(i3msk6) = 653555
      frequence(i3msk6) = 57.968208435
      equitide(i3msk6) = 0.0
!
!  M6
!------------------
!
      im6=103
      nommar(im6) ='M6'
      doodson(im6) = 655555
      frequence(im6) = 57.968208435
      equitide(im6) = 0.0
!
!  3MKS6
!------------------
!
      i3mks6=104
      nommar(i3mks6) ='3MKS6'
      doodson(i3mks6) = 657555
      frequence(i3mks6) = 57.968208435
      equitide(i3mks6) = 0.0
!
!  MSN6
!------------------
!
      imsn6=105
      nommar(imsn6) ='MSN6'
      doodson(imsn6) = 663655
      frequence(imsn6) = 57.968208435
      equitide(imsn6) = 0.0
!
!  4MN.2ML6
!------------------
!
      i4mn2ml6=106
      nommar(i4mn2ml6) ='4MN.2ML6'
      doodson(i4mn2ml6) = 665455
      frequence(i4mn2ml6) = 57.968208435
      equitide(i4mn2ml6) = 0.0
!
!  MNK6
!------------------
!
      imnk6=107
      nommar(imnk6) ='MNK6'
      doodson(imnk6) = 665655
      frequence(imnk6) = 57.968208435
      equitide(imnk6) = 0.0
!
!  2MT6
!------------------
!
      i2mt6=108
      nommar(i2mt6) ='2MT6'
      doodson(i2mt6) = 672556
      frequence(i2mt6) = 57.968208435
      equitide(i2mt6) = 0.0
!
!  2MS6
!------------------
!
      i2ms6=109
      nommar(i2ms6) ='2MS6'
      doodson(i2ms6) = 673555
      frequence(i2ms6) = 57.968208435
      equitide(i2ms6) = 0.0
!
!  2MK6
!------------------
!
      i2mk6=110
      nommar(i2mk6) ='2MK6'
      doodson(i2mk6) = 675555
      frequence(i2mk6) = 57.968208435
      equitide(i2mk6) = 0.0
!
!  2SN6
!------------------
!
      i2sn6=111
      nommar(i2sn6) ='2SN6'
      doodson(i2sn6) = 681655
      frequence(i2sn6) = 57.968208435
      equitide(i2sn6) = 0.0
!
!  3MSN6
!------------------
!
      i3msn6=112
      nommar(i3msn6) ='3MSN6'
      doodson(i3msn6) = 683455
      frequence(i3msn6) = 57.968208435
      equitide(i3msn6) = 0.0
!
!  3MKN6
!------------------
!
      i3mkn6=113
      nommar(i3mkn6) ='3MKN6'
      doodson(i3mkn6) = 685455
      frequence(i3mkn6) = 57.968208435
      equitide(i3mkn6) = 0.0
!
!  2SM6
!------------------
!
      i2sm6=114
      nommar(i2sm6) ='2SM6'
      doodson(i2sm6) = 691555
      frequence(i2sm6) = 57.968208435
      equitide(i2sm6) = 0.0
!
!  MSK6
!------------------
!
      imsk6=115
      nommar(imsk6) ='MSK6'
      doodson(imsk6) = 693555
      frequence(imsk6) = 57.968208435
      equitide(imsk6) = 0.0
!
!  2MNS8
!------------------
!
      i2mns8=116
      nommar(i2mns8) ='2MNS8'
      doodson(i2mns8) = 827655
      frequence(i2mns8) = 57.968208435
      equitide(i2mns8) = 0.0
!
!  2(MN)8
!------------------
!
      i2mn8=117
      nommar(i2mn8) ='2_MN_8'
      doodson(i2mn8) = 835755
      frequence(i2mn8) = 57.968208435
      equitide(i2mn8) = 0.0
!
!  3MN8
!------------------
!
      i3mn8=118
      nommar(i3mn8) ='3MN8'
      doodson(i3mn8) = 845655
      frequence(i3mn8) = 57.968208435
      equitide(i3mn8) = 0.0
!
!  3MNU8
!------------------
!
      i3mnu8=119
      nommar(i3mnu8) ='3MNU8'
      doodson(i3mnu8) = 847455
      frequence(i3mnu8) = 57.968208435
      equitide(i3mnu8) = 0.0
!
!  2MSK8
!------------------
!
      i2msk8=120
      nommar(i2msk8) ='2MSK8'
      doodson(i2msk8) = 853555
      frequence(i2msk8) = 57.968208435
      equitide(i2msk8) = 0.0
!
!  M8
!------------------
!
      im8=121
      nommar(im8) ='M8'
      doodson(im8) = 855555
      frequence(im8) = 57.968208435
      equitide(im8) = 0.0
!
!  4MKS8
!------------------
!
      i4mks8=122
      nommar(i4mks8) ='4MKS8'
      doodson(i4mks8) = 857555
      frequence(i4mks8) = 57.968208435
      equitide(i4mks8) = 0.0
!
!  2MSN8
!------------------
!
      i2msn8=123
      nommar(i2msn8) ='2MSN8'
      doodson(i2msn8) = 863655
      frequence(i2msn8) = 57.968208435
      equitide(i2msn8) = 0.0
!
!  3ML8
!------------------
!
      i3ml8=124
      nommar(i3ml8) ='3ML8'
      doodson(i3ml8) = 865455
      frequence(i3ml8) = 57.968208435
      equitide(i3ml8) = 0.0
!
!  2MNK8
!------------------
!
      i2mnk8=125
      nommar(i2mnk8) ='2MNK8'
      doodson(i2mnk8) = 865655
      frequence(i2mnk8) = 57.968208435
      equitide(i2mnk8) = 0.0
!
!  3MS8
!------------------
!
      i3ms8=126
      nommar(i3ms8) ='3MS8'
      doodson(i3ms8) = 873555
      frequence(i3ms8) = 57.968208435
      equitide(i3ms8) = 0.0
!
!  3MK8
!------------------
!
      i3mk8=127
      nommar(i3mk8) ='3MK8'
      doodson(i3mk8) = 875555
      frequence(i3mk8) = 57.968208435
      equitide(i3mk8) = 0.0
!
!  2SMN8
!------------------
!
      i2smn8=128
      nommar(i2smn8) ='2SMN8'
      doodson(i2smn8) = 881655
      frequence(i2smn8) = 57.968208435
      equitide(i2smn8) = 0.0
!
!  4MSN8
!------------------
!
      i4msn8=129
      nommar(i4msn8) ='4MSN8'
      doodson(i4msn8) = 883455
      frequence(i4msn8) = 57.968208435
      equitide(i4msn8) = 0.0
!
!  MSNK8
!------------------
!
      imsnk8=130
      nommar(imsnk8) ='MSNK8'
      doodson(imsnk8) = 883655
      frequence(imsnk8) = 57.968208435
      equitide(imsnk8) = 0.0
!
!  4MNK8
!------------------
!
      i4mnk8=131
      nommar(i4mnk8) ='4MNK8'
      doodson(i4mnk8) = 885455
      frequence(i4mnk8) = 57.968208435
      equitide(i4mnk8) = 0.0
!
!  2(MS)8
!------------------
!
      i2ms8=132
      nommar(i2ms8) ='2_MS_8'
      doodson(i2ms8) = 891555
      frequence(i2ms8) = 57.968208435
      equitide(i2ms8) = 0.0
!
!  2MSK8
!------------------
!
      i2msk82=133
      nommar(i2msk82) ='2MSK8'
      doodson(i2msk82) = 893555
      frequence(i2msk82) = 57.968208435
      equitide(i2msk82) = 0.0
!
!  5MNS10
!------------------
!
      i5mns10=134
      nommar(i5mns10) ='5MNS10'
      doodson(i5mns10) =1027655
      frequence(i5mns10) = 57.968208435
      equitide(i5mns10) = 0.0
!
!  3M2N10
!------------------
!
      i3m2n10=135
      nommar(i3m2n10) ='3M2N10'
      doodson(i3m2n10) =1035755
      frequence(i3m2n10) = 57.968208435
      equitide(i3m2n10) = 0.0
!
!  4MN10
!------------------
!
      i4mn10=136
      nommar(i4mn10) ='4MN10'
      doodson(i4mn10) =1045655
      frequence(i4mn10) = 57.968208435
      equitide(i4mn10) = 0.0
!
!  M10
!------------------
!
      im10=137
      nommar(im10) ='M10'
      doodson(im10) =1055555
      frequence(im10) = 57.968208435
      equitide(im10) = 0.0
!
!  3MSN10
!------------------
!
      i3msn10=138
      nommar(i3msn10) ='3MSN10'
      doodson(i3msn10) =1063655
      frequence(i3msn10) = 57.968208435
      equitide(i3msn10) = 0.0
!
!  4MS10
!------------------
!
      i4ms10=139
      nommar(i4ms10) ='4MS10'
      doodson(i4ms10) =1073555
      frequence(i4ms10) = 57.968208435
      equitide(i4ms10) = 0.0
!
!  4MK10
!------------------
!
      i4mk10=140
      nommar(i4mk10) ='4MK10'
      doodson(i4mk10) =1075555
      frequence(i4mk10) = 57.968208435
      equitide(i4mk10) = 0.0
!
!  5MSN10
!------------------
!
      i5msn10=141
      nommar(i5msn10) ='5MSN10'
      doodson(i5msn10) =1083455
      frequence(i5msn10) = 57.968208435
      equitide(i5msn10) = 0.0
!
!  2MSNK10
!------------------
!
      i2msnk10=142
      nommar(i2msnk10) ='2MSNK10'
      doodson(i2msnk10) =1083655
      frequence(i2msnk10) = 57.968208435
      equitide(i2msnk10) = 0.0
!
!  3M2S10
!------------------
!
      i3m2s10=143
      nommar(i3m2s10) ='3M2S10'
      doodson(i3m2s10) =1091555
      frequence(i3m2s10) = 57.968208435
      equitide(i3m2s10) = 0.0
!
!-------------------------------------------------
!
      vathau=347.80925080
!     vathau=347.8092506160
!     vathau=360./1.03505
      valonmlun=13.17639673
!     valonmlun=13.1763967440
!     valonmlun=360/27.321582
      valonmsol=0.98564734
!     valonmsol=0.985647359999986
!     valonmsol=0.98564736
!     valonmsol=360./365.242199
      valonperl=0.11140408
!     valonperl=360./8.847309/365.242199
      valonascl=0.05295392
!     valonascl=360./18.612904/365.242199
      valonpers=0.00004707
!     valonpers=360./209.4021/36524.2199
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
        asfo=vathau*nd_tmp(1)+valonmlun*(nd_tmp(2)-5)+
     &       valonmsol*(nd_tmp(3)-5)+valonperl*(nd_tmp(4)-5)
     &       +valonascl*(nd_tmp(5)-5)+valonpers*(nd_tmp(6)-5)
        frequence(i)=max(mod(asfo/24.0,360.0),1e-10)
        omega(i)=2*pi*frequence(i)*(pi/180.0/3600.0)
       end do
        end subroutine tide_data


        subroutine mas_init_ssh(Istr,Iend,Jstr,Jend,IstrR,IendR,
     &                          JstrR,JendR,Ntides,k_of,nn_of,
     &                          SSH_Tamp,SSH_Tphase,mrssh,mgssh,mhssh)
!
!  Convert the standard SSH_Tamp/SSH_Tphase harmonic amplitude/phase
!  into SHOM reduced amplitude/phase tables, at every grid point.
!
        implicit none
        integer, intent(in) :: Istr,Iend,Jstr,Jend,IstrR,IendR,
     &                          JstrR,JendR,Ntides
        integer, intent(in) :: k_of(:),nn_of(:)
        real, intent(in)    :: SSH_Tamp(GLOBAL_2D_ARRAY,Ntides)
        real, intent(in)    :: SSH_Tphase(GLOBAL_2D_ARRAY,Ntides)
        real, intent(inout)   :: mrssh(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mgssh(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mhssh(130,GLOBAL_2D_ARRAY)
        integer :: i,j,ht,k,nn,ir0

        mhssh(129,:,:) = -999.0
        DO ht=1,Ntides
          IF (l_presence(nu(ht))) THEN
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
          ENDIF
        ENDDO
        end subroutine mas_init_ssh


        subroutine mas_init_pot(Istr,Iend,Jstr,Jend,IstrR,IendR,
     &                          JstrR,JendR,Ntides,k_of,nn_of,
     &                          POT_Tamp,POT_Tphase,mrpot,mgpot,mhpot)
!
        implicit none
        integer, intent(in) :: Istr,Iend,Jstr,Jend,IstrR,IendR,
     &                          JstrR,JendR,Ntides
        integer, intent(in) :: k_of(:),nn_of(:)
        real, intent(in)    :: POT_Tamp(GLOBAL_2D_ARRAY,Ntides)
        real, intent(in)    :: POT_Tphase(GLOBAL_2D_ARRAY,Ntides)
        real, intent(inout)   :: mrpot(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mgpot(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mhpot(130,GLOBAL_2D_ARRAY)

        call mas_init_ssh(Istr,Iend,Jstr,Jend,IstrR,IendR,
     &                     JstrR,JendR,Ntides,k_of,nn_of,
     &                     POT_Tamp,POT_Tphase,mrpot,mgpot,mhpot)
        end subroutine mas_init_pot


# if defined UV_TIDES
        subroutine mas_init_uv(Istr,Iend,Jstr,Jend,IstrR,IendR,
     &                         JstrR,JendR,Ntides,k_of,nn_of,
     &                         UV_Tmajor,UV_Tminor,UV_Tangle,UV_Tphase,
     &                         angler,mru,mgu,mrv,mgv,
     &                         mhu,mhv)
!
!  Decompose the U/V tidal ellipses into complex U/V harmonics
!  (ep2ap), then convert their rho-to-u/v averages into SHOM
!  reduced amplitude/phase tables, at every u- and v-point.
!
        implicit none
        integer, intent(in) :: Istr,Iend,Jstr,Jend,IstrR,IendR,
     &                          JstrR,JendR,Ntides
        integer, intent(in) :: k_of(:),nn_of(:)
        real, intent(in)    :: UV_Tmajor(GLOBAL_2D_ARRAY,Ntides)
        real, intent(in)    :: UV_Tminor(GLOBAL_2D_ARRAY,Ntides)
        real, intent(in)    :: UV_Tangle(GLOBAL_2D_ARRAY,Ntides)
        real, intent(in)    :: UV_Tphase(GLOBAL_2D_ARRAY,Ntides)
        real, intent(in)    :: angler(GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mru(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mgu(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mrv(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mgv(30,8,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mhu(130,GLOBAL_2D_ARRAY)
        real, intent(inout)   :: mhv(130,GLOBAL_2D_ARRAY)
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
 
        DO ht=1,Ntides
          IF (l_presence(nu(ht))) THEN
            ir0=rr0(nu(ht))
            k =k_of(ht)
            nn=nn_of(ht)

            do j=Jstr-1,JendR
              do i=Istr-1,IendR
                angle=UV_Tangle(i,j,ht)
#  ifdef CURVGRID
     &                -angler(i,j)*rad2deg
#  endif
                call ep2ap(UV_Tmajor(i,j,ht),UV_Tminor(i,j,ht),
     &                     angle,UV_Tphase(i,j,ht),
     &                     u_cmplx(i,j),v_cmplx(i,j))
              enddo
            enddo
            do j=JstrR,JendR
              do i=Istr,IendR
                U_interp = 0.5*(u_cmplx(i,j)+u_cmplx(i-1,j))
                mru(nn,k,i,j) = abs(U_interp)*100.0
                mgu(nn,k,i,j) = atan2(aimag(U_interp),real(U_interp))
     &                         *rad2deg
                if (mgu(nn,k,i,j) < 0.0) mgu(nn,k,i,j)=mgu(nn,k,i,j)
     &                                                +360.0
                if (ir0 < 0) mgu(nn,k,i,j)=mgu(nn,k,i,j)+180.0
              enddo
            enddo
            do j=Jstr,JendR
              do i=IstrR,IendR
                V_interp = 0.5*(v_cmplx(i,j)+v_cmplx(i,j-1))
                mrv(nn,k,i,j) = abs(V_interp)*100.0
                mgv(nn,k,i,j) = atan2(aimag(V_interp),real(V_interp))
     &                         *rad2deg
                if (mgv(nn,k,i,j) < 0.0) mgv(nn,k,i,j)=mgv(nn,k,i,j)
     &                                                +360.0
                if (ir0 < 0) mgv(nn,k,i,j)=mgv(nn,k,i,j)+180.0
              enddo
            enddo
          ENDIF
        ENDDO
        end subroutine mas_init_uv
# endif


        subroutine mas_update_calendar(time)
!
!  Decomposes `time` into a calendar date and flags whether a new
!  day has started since the last call, publishing the result via
!  the module's mas_ijour/mas_imois/mas_ian/mas_tempis/mas_is_new_day.
!  Must be called exactly once per time step, before any tile reads
!  those values (e.g. under OMP MASTER + BARRIER, ahead of the tile
!  loop that calls clm_tides): this is the only place jourmem is
!  touched, so calling it per-tile instead would race across threads.
!
        implicit none
        real, intent(in) :: time
        CHARACTER(len=19) :: cdate,tool_sectodat
        integer :: iheure,iminu,isec

        cdate = tool_sectodat(time)
        CALL tool_decompdate(cdate,mas_ijour,mas_imois,mas_ian,
     &                        iheure,iminu,isec)
        mas_tempis = iheure*3600.0+iminu*60.0+isec
        mas_is_new_day = (mas_ijour /= jourmem)
        if (mas_is_new_day) jourmem = mas_ijour
        end subroutine mas_update_calendar


        subroutine mas_daily_update(imin,imax,jmin,jmax,
     &                               new_day,ijour,imois,ian,r,g,h)
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
              call ma1(ijour,imois,ian,
     &                 h(:,i,j),r(:,:,i,j),g(:,:,i,j))
            enddo
          enddo
        endif
        end subroutine mas_daily_update


        subroutine mas_interpolate(imin,imax,jmin,jmax,
     &                              Istr,Iend,Jstr,Jend,tempis,h,X)
!
!  Interpolate the current 128-heights/day table to the current time
!  of day (ma7), at every grid point, for a single h/X field pair
!  (e.g. mhssh/Etide, or mhu/Utide). Call once per active field.
!
!  X is a per-tile PRIVATE_2D_SCRATCH_ARRAY (Istr-2:Iend+2,Jstr-2:
!  Jend+2), unlike h which is a persistent GLOBAL_2D_ARRAY: its
!  bounds are a runtime offset of the caller's own tile, so Istr/
!  Iend/Jstr/Jend (not just imin/imax/jmin/jmax) must be passed in
!  to declare it with matching bounds.
!
        implicit none
        integer :: i,j
        integer, intent(in) :: imin,imax,jmin,jmax
        integer, intent(in) :: Istr,Iend,Jstr,Jend
        real, intent(in)    :: tempis
        real, intent(in)    :: h(130,GLOBAL_2D_ARRAY)
        real, intent(out) :: X(PRIVATE_2D_SCRATCH_ARRAY)

         do j=jmin,jmax
          do i=imin,imax
            call ma7(tempis,h(:,i,j),X(i,j))
          enddo
        enddo
        end subroutine mas_interpolate


# if defined UV_TIDES
        subroutine ep2ap(semi_major, semi_minor, inclination, phase,
     &                    u_complex, v_complex)
     !
     ! Compute U/V amplitude and phase from ellipse parameters
     !
        implicit none
        real(8), intent(in)  :: semi_major, semi_minor
        real(8), intent(in)  :: inclination, phase
        complex, intent(out) :: u_complex, v_complex 

        complex(8) :: wp, wm
        real(8) :: incl_rad, phase_rad

        real pi, deg2rad, rad2deg
        parameter (pi=3.14159265358979323846, deg2rad=pi/180.,
     &                                      rad2deg=180./pi)

        ! Convert to radians
        incl_rad = inclination * deg2rad
        phase_rad = phase * deg2rad

        ! Decompose ellipse into counter-rotating circles
        wp = (semi_major + semi_minor) / 2.0d0 *
     &       exp(cmplx(0.0d0,-1.D0)*(incl_rad - phase_rad))
        wm = (semi_major - semi_minor) / 2.0d0 *
     &       exp(cmplx(0.0d0,-1.D0)*(incl_rad + phase_rad))

        ! Reconstruct U and V components
        u_complex = wp + conjg(wm)
        v_complex = cmplx(0.0d0, 1.0d0) * (wp - conjg(wm))

        end subroutine ep2ap
# endif
!******************************************************
!
! Interpolations entres les 128 hauteurs precalculees
!
        subroutine ma7(seci,hn,haut)
        implicit none
        integer                      :: idk,l,k
        real               :: d1,haut,fk,d,hh,d2
        real,dimension(130):: hn
        real               :: seci

        fk=seci/60.0/11.25+2.0
        k=fk+0.5
        d=fk-k
        if(k >= 130)then
         idk=k-129
         k=k-idk
         d=d+idk
        endif
        l=k
        if(l == 130)l=129
        hh=hn(l)+hn(l)
        d1=hn(l+1)-hn(l-1)
        d2=hn(l+1)+hn(l-1)-hh

        haut=(hh+d*(d*d2+d1))/800.0
        end subroutine ma7



        subroutine ma1(jour,mois,ia,hn,rr,gg)

!*****  calcul de hn :128 hauteurs reduite de maree par jour
        implicit none
        integer :: l,nj,i,nd0,n,m,jour,mois,ia,kg,ng,j,k
        integer :: ll,kk
        integer,dimension(30)             :: noa
        real                    :: fq,co,so
        double precision        :: fnjd
        real,dimension(30)      :: f,q,q0,r1,r2,v0
        real,dimension(31)      :: g,r
        real,dimension(130)     :: hn
        real,dimension(362)     :: c
        real,dimension(91)      :: cs=(/
     &
     &  1.000000,0.999848,0.999391,0.998630,0.997564,0.996195,0.994522,
     &  0.992546,0.990268,0.987688,0.984808,0.981627,0.978148,0.974370,
     &  0.970296,0.965926,0.961262,0.956305,0.951057,0.945519,0.939693,
     &  0.933580,0.927184,0.920505,0.913545,0.906308,0.898794,0.891007,
     &  0.882948,0.874620,0.866025,0.857167,0.848048,0.838671,0.829038,
     &  0.819152,0.809017,0.798636,0.788011,0.777146,0.766044,0.754710,
     &  0.743145,0.731354,0.719340,0.707107,0.694658,0.681998,0.669131,
     &  0.656059,0.642788,0.629320,0.615661,0.601815,0.587785,0.573576,
     &  0.559193,0.544639,0.529919,0.515038,0.500000,0.484810,0.469472,
     &  0.453991,0.438371,0.422618,0.406737,0.390731,0.374607,0.358368,
     &  0.342020,0.325568,0.309017,0.292372,0.275637,0.258819,0.241922,
     &  0.224951,0.207912,0.190809,0.173648,0.156434,0.139173,0.121869,
     &  0.104528,0.087156,0.069756,0.052336,0.034899,0.017452,0.000000/)
        real,dimension(30,8)    :: rr,gg
        real,dimension(11,3)    :: x,y
        equivalence (c(1),cs(1))

        l=1

        x(:,:)=0.0
        y(:,:)=0.0
   
        c(271)=0.0
        do i=1,90
         c(182-i)=-c(i)
         c(362-i)=c(i)
         c(180+i)=-c(i)
        end do
        nj=julien2(ia,jour,mois)-2415021
        fnjd=real(nj+l,8)-DBLE(1.4921875)
        nj=mod(nj,7)+1
        if(nj <= 0)nj=nj+7
        nd0=-1
        n=0
        m=1

        do kg=1,8
         n=num(kg)
         if(n /= 0) then
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
           if(q0(i) < 0.0)q0(i)=q0(i)+360.0
          end do
          do i=l,3
           do j=1,n
            fq=q0(j)
            k=fq+1
            fq=fq-k+1
            co=((c(k+1)-c(k))*fq+c(k))*r(j)
            k=k-90
            if(k <= 0) k=k+360
            so=((c(k+1)-c(k))*fq+c(k))*r(j)
            x(ng,i)=x(ng,i)+co
            y(ng,i)=y(ng,i)+so
            q0(j)=q0(j)+q(j)
            if(q0(j) > 360.0)q0(j)=q0(j)-360.0
           end do
          end do
         end if
        end do
        call ma4(hn,l,x,y)
        end subroutine ma1


        subroutine ma4(hn,l,hm1,hm2)
!  Calcul des hauteurs reduites h pour les jours jour-1,jour et jour+1
!  et appel de la subroutine de calcul des 128 hauteurs par jour (ma5)
        implicit none
        integer                        :: k,l,i
        real                 :: fnm
        real,dimension(130)  :: hn
        real,dimension(11,3) :: hm1,hm2
        real,dimension(128,3):: h
        complex,dimension(128)         :: x

        do  k=l,3
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
         h(:,k)=x(:)
        end do
        call ma5(h,hn)
        end subroutine ma4


        subroutine ma5(h,x)
!
!  Calcul des 128 hauteurs par interpolation entre les hauteurs reduites
!
        implicit none
        integer                        :: k
        real                 :: b,hh
        real,dimension(128,3):: h
        real,dimension(130)  :: x

        b=-0.5078125

        x(1)=x(129)
        x(2)=x(130)

        do k=1,128
         hh=h(k,2)+h(k,2)
         b=b+0.0078125
         x(k+2)=((h(k,1)+h(k,3)-hh)*b+h(k,3)-h(k,1))*b+hh
        end do

        if(x(1) == -valmanq)then
         x(2)=3.0*(x(3)-x(4))+x(5)
         x(1)=6.0*x(3)-8.0*x(4)+3.0*x(5)
        endif

        end subroutine ma5

        subroutine fft(a)
!
!  Transformee de Fourier rapide
!
        implicit none
        integer                :: m,nv2,nm1,n,i,j,ip,l,k,le,le1
        real         :: fij
        complex,dimension(128) :: a
        complex                :: u,w,t
        m=7
        nv2=64
        nm1=127
        n=128
        j=1
        do  i=1,nm1
         if(i < j) then
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
         w=cexp(cmplx(0.0,fij))
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
!  Calcul des facteurs nodaux (amplitudes et phases)
!
        implicit none
        integer                     :: i,ng,n
        integer,dimension(6)        :: nd
        integer,dimension(30)       :: noa
        real              :: fa1,fb1,fa2,fb2,a1,a2
        double precision            :: fnjd
        real,dimension(30):: q,v0,f,r,g,r1,r2
        double precision            :: dq
        complex                     :: fc
   
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
         if(nd(1) == 0) r(i)=r(i)*2.0
         call masfo(nd,dq)
         q(i)=dq
         g(i)=mod(g(i),360.0)
         if(g(i) < 0.0) g(i)=g(i)+360.0
         v0(i)=asfov(nd)+mod(REAL(dq*fnjd,8),DBLE(360.))
         fc=1.0
         if( r1(i) /= 0.0 .or. r2(i) /= 0.0 ) then
          fb2=9.242202e-4
          fb1=-fb2
          fa2=1.760045
          fa1=4.523139
          if(noa(i) == 275555)then
           fb1=2*fb2
           fa1=2*fa2
          endif
          a1=fa1+fb1*fnjd
          a2=fa2+fb2*fnjd
          fc=1.0+r1(i)*cmplx(cos(a1),sin(a1))+
     &    r2(i)*cmplx(cos(a2),sin(a2))
        end if
  
        v0(i)=v0(i)+atan2(aimag(fc),real(fc))*57.29578
        f(i)=cabs(fc)
        if(noa(i) == 355555.or.noa(i) == 382555 .or.noa(i) == 164555)
     &          v0(i)=v0(i)-90.0
        end do
        end subroutine ma3
 
        subroutine masfo(nd,asfo)
!
! calcul des vitesses angulaires des composantes
!
        implicit none
        integer,dimension(6) :: nd
        double precision       :: asfo
   
        asfo =(DBLE(360.)-DBLE(12.19074939))*(nd(1))
     &     +DBLE(13.17639673)*(nd(2)-5)+DBLE(0.98564734)*(nd(3)-5)
     &     +DBLE(0.11140408)*(nd(4)-5)+DBLE(0.05295392)*(nd(5)-5)+
     &      DBLE(0.00004707)*(nd(6)-5)
        end subroutine masfo

        function asfov(nd)
!
! Calcul des arguments astronomiques
!
        implicit none
        integer,dimension(6) :: nd
        real       :: fov,asfov
   
        fov=280.1895*(nd(1)+nd(3)-5)+
     &  277.0248*(nd(2)-nd(1)-5)+334.3853*(nd(4)-5)
     &  +100.8432*(nd(5)-5)+281.2209*(nd(6)-5)+(mod(nd(1),2))*90.0
        asfov=(mod(fov,360.0))
        end function asfov

        function julien2(ia,jou,moi)
!
! Calcul du nombre de jours ecoules depuis le 1er Janvier 1900
!
        implicit none
        integer                 :: i,julien2,ia,jou,moi,ibs,jour,mois,iy,m
        integer,dimension(12,2),PARAMETER :: jo=
     &     reshape((/0,31,59,90,120,151,181,212,243,273,304,334,
     &     0,31,60,91,121,152,182,213,244,274,305,335/),(/12,2/))
        real          :: a,b

        jour=jou
        mois=moi
        if(mois==1) then
        ibs=1
        if( ia >= 1582 .and.(ia /= 1582 .or. jou >= 277) ) then
         ibs=mod(ia,4)+2
         if(ibs /= 2)ibs=1
        end if
        do i=1,12
         mois=12-i+1
         jour=jou-jo(mois,ibs)
         if(jour > 0) exit
        end do
        end if

        a=ia+mois/100.0+jour/10000.0
        b=0.0
        iy=ia
        m=mois

        if(mois <= 2)then
         iy=iy-1
         m=m+12
        endif

        if(a >= 1582.1015)then
         b=2-int(iy/100)+int(iy/100)/4
        endif

        julien2=int(365.25*iy)+int(30.6001*(m+1))+jour+1720995+b
        end function julien2
#endif   
      END MODULE MOD_TIDES_MAS

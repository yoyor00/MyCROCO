!======================================================================
! CROCO is derived from the ROMS-AGRIF branch of ROMS.
! ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
! from the UCLA branch (Shchepetkin et al.) and the Rutgers
! University branch (Arango et al.), both under MIT/X style license.
! Copyright (C) 2005-2026 CROCO Development Team
! License: CeCILL-2.1 - see LICENSE.txt
!
! CROCO website : https://www.croco-ocean.org
!======================================================================
!
! Parameters not modified by USER
!
!======================================================================

      integer  LLm,Lm,MMm,Mm
#if defined AGRIF
      integer LLmm2, MMmm2
#endif
#ifdef AGRIF
      common /scrum_physical_grid/ LLm,Lm,LLmm2,MMm,Mm,MMmm2
#else
      parameter (LLm=LLm0,  MMm=MMm0)
#endif

!
!----------------------------------------------------------------------
! MPI related variables
!----------------------------------------------------------------------
      integer Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
#ifdef MPI
#ifdef OPENACC
      integer my_acc_device
      logical compute_on_device
      common/comm_my_device/my_acc_device,compute_on_device
#endif
#elif defined OPENMP
#else
#ifdef OPENACC
      integer, parameter :: my_acc_device = 0
      logical, parameter :: compute_on_device = .true.
#endif
#endif


!
!----------------------------------------------------------------------
! Number maximum of weights for the barotropic mode
!----------------------------------------------------------------------
!
      integer, parameter :: NWEIGHT = 1000

!
!----------------------------------------------------------------------
! Derived dimension parameters
!----------------------------------------------------------------------
!
      integer stdout, padd_X,padd_E
#ifdef AGRIF
      common/scrum_deriv_param/padd_X,padd_E
#endif
#ifdef LOGFILE
      common /stdout/stdout
#else
      parameter (stdout=6)
#endif
      integer, parameter :: Np = N+1
      integer, parameter :: NpHz = (N+1+N_sl)
#ifndef AGRIF
# ifdef MPI
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
# else
      parameter (Lm=LLm, Mm=MMm)
# endif
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
#endif

#if defined AGRIF
      integer N2d,N3d,N3dHz,N1dXI,N1dETA
# if !defined NBQ
      integer, parameter :: NSA = 28
# else
      integer, parameter :: NSA = 35
# endif
      common /scrum_private_param/ N2d,N3d,N3dHz,N1dXI,N1dETA
#else
# if !defined NBQ
      integer, parameter :: NSA = 28
# else
      integer, parameter :: NSA = 35
# endif
# ifdef ALLOW_SINGLE_BLOCK_MODE
      integer, parameter :: size_XI = 6+Lm
      integer, parameter :: size_ETA = 6+Mm
# else
      integer, parameter :: size_XI = 7+(Lm+NSUB_X-1)/NSUB_X
      integer, parameter :: size_ETA = 7+(Mm+NSUB_E-1)/NSUB_E
# endif
      integer, parameter :: sse = size_ETA/Np
      integer, parameter :: ssz = Np/size_ETA
      integer, parameter :: se = sse/(sse+ssz)
      integer, parameter :: sz = 1-se
      integer, parameter :: N2d = size_XI*(se*size_ETA+sz*Np)
      integer, parameter :: N3d = size_XI*size_ETA*Np
      integer, parameter :: N3dHz = size_XI*size_ETA*NpHz
# ifdef ABL1D
      integer, parameter :: sse_abl = size_ETA/(N_abl+1)
      integer, parameter :: ssz_abl = (N_abl+1)/size_ETA
      integer, parameter :: se_abl = sse_abl/(sse_abl+ssz_abl)
      integer, parameter :: sz_abl = 1-se_abl
      integer, parameter :: N2dabl = size_XI*(se_abl*size_ETA+sz_abl*(N_abl+1))
      integer, parameter :: N3dabl = size_XI*size_ETA*(N_abl+1)
# endif
#endif

!
!----------------------------------------------------------------------
! I/O : flag for type sigma vertical transformation
!----------------------------------------------------------------------
!
#ifdef NEW_S_COORD
      real, parameter :: Vtransform = 2
#else
      real, parameter :: Vtransform = 1
#endif


#ifdef STATIONS
      ! Maximum number of stations
      integer, parameter :: Msta = 1000
#endif

!
!----------------------------------------------------------------------
! Number of tracers
!----------------------------------------------------------------------
!
#ifdef SOLVE3D
!
# ifdef TEMPERATURE
      integer, parameter :: itemp = 1
      integer, parameter :: ntrc_temp = 1
# else
      integer, parameter :: itemp = 0
      integer, parameter :: ntrc_temp = 0
# endif
# ifdef SALINITY
      integer, parameter :: ntrc_salt = 1
# else
      integer, parameter :: ntrc_salt = 0
# endif
# ifdef STOGEN
      integer, parameter :: ntrc_stogen = 3
# else
      integer, parameter :: ntrc_stogen = 0
#endif

# ifdef BIOLOGY
#  ifdef PISCES
#   ifdef key_pisces_npzd
         integer, parameter :: ntrc_bio = 9
#   elif defined key_pisces_quota
#    ifdef key_ligand
         integer, parameter :: ntrc_bio = 40
#    else
         integer, parameter :: ntrc_bio = 39
#    endif
#   else
#    ifdef key_ligand
         integer, parameter :: ntrc_bio = 25
#    else
         integer, parameter :: ntrc_bio = 24
#    endif
#   endif
#  elif defined BIO_NChlPZD
#   ifdef OXYGEN
      integer, parameter :: ntrc_bio = 6
#   else
      integer, parameter :: ntrc_bio = 5
#   endif
#  elif defined BIO_N2ChlPZD2
      integer, parameter :: ntrc_bio = 7
#  elif defined BIO_BioEBUS
#   ifdef NITROUS_OXIDE
      integer, parameter :: ntrc_bio = 12
#   else
      integer, parameter :: ntrc_bio = 11
#   endif
#  endif
# else
      integer, parameter :: ntrc_bio = 0
# endif /* BIOLOGY */

# ifdef GLS_MIXING
      integer, parameter :: NGLS = 2
      integer, parameter :: itke = 1
      integer, parameter :: igls = 2
# endif

#endif /* SOLVE3D */



!
!----------------------------------------------------------------------
! Tracer identification indices
!----------------------------------------------------------------------
!
#if defined SOLVE3D
      integer   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer   ntrc_diaeddy, ntrc_surf
# ifdef BIOLOGY
     &          , itrc_bio
# endif
# ifdef SEDIMENT
     &          , itrc_sed, itrc_sand, itrc_mud, itrc_grav
# endif
# ifdef SALINITY
     &          , isalt
# endif
# ifdef PASSIVE_TRACER
     &          , itpas
# endif
!
# ifdef BIOLOGY
#  ifdef PISCES
     &          , iDIC_, iTAL_, iOXY_, iCAL_, iPO4_
     &          , iPOC_, iSIL_, iPHY_, iZOO_, iDOC_
     &          , iDIA_, iMES_, iDSI_, iFER_
     &          , iBFE_, iGOC_, iSFE_, iDFE_, iGSI_
     &          , iNFE_, iNCH_, iDCH_, iNO3_, iNH4_
     &          , iLGW_, iDON_, iDOP_, iPON_, iPOP_
     &          , iNPH_, iPPH_, iNDI_, iPDI_, iPIC_
     &          , iNPI_, iPPI_, iPFE_, iPCH_, iGON_
     &          , iGOP_
#   ifdef DIAGNOSTICS_BIO
#    ifdef key_trc_diaadd
     &          , Nhi,Nco3,Naksp,Netot,Nprorca
     &          , Nprorcad,Npronew,Npronewd
     &          , Nprobsi,Nprofed,Nprofen
     &          , Ngrapoc,Ngrapoc2
     &          , Nmico2,Nmeso2
     &          , Nnitrifo2,Nfixo2,Nremino2
     &          , Npronewo2,Nprorego2
     &          , Nfld,Nflu16,Nkgco2,Natcco2,Nsinking
     &          , Nsinkfer,Nsinksil,Nironsed
     &          , Nsinkcal,Nheup,Nnitrpot
     &          , Nirondep,Nsildep,Npo4dep
     &          , Nno3dep,Nnh4dep
#    endif
#   endif
     &          , NumFluxTerms,NumVSinkTerms,NumGasExcTerms
#  elif defined BIO_NChlPZD
     &          , iNO3_, iChla, iPhy1, iZoo1
     &          , iDet1
#   ifdef OXYGEN
     &          , iO2
#   endif
     &          , NFlux_NewProd, NFlux_Grazing, NFlux_SlopFeed
     &          , NFlux_Pmort, NFlux_Zmetab, NFlux_Zmort, NFlux_ReminD
     &          , NumFluxTermsN
#   ifdef OXYGEN
     &          , OGain_NewProd, OLoss_Zmetab
     &          , OLoss_ReminD, NumFluxTermsO, OFlux_GasExc
     &          , NumGasExcTerms, NumFluxTerms
#   else
     &          , NumGasExcTerms
     &          , NumFluxTerms
#   endif
     &          , NFlux_VSinkP1, NFlux_VSinkD1
     &          , NumVSinkTerms
#  elif defined BIO_N2ChlPZD2
     &          , iNO3_, iNH4_, iChla, iPhy1, iZoo1
     &          , iDet1, iDet2
     &          , NFlux_NewProd
     &          , NFlux_RegProd
     &          , NFlux_Nitrific
     &          , NFlux_Grazing
     &          , NFlux_SlopFeed
     &          , NFlux_Pmort
     &          , NFlux_Zmetab
     &          , NFlux_Zmort
     &          , NFlux_ReminD1, NFlux_ReminD2
     &          , NFlux_CoagPhy, NFlux_CoagSDet
     &          , NumFluxTermsN, NumFluxTerms, NumGasExcTerms
     &          , NFlux_VSinkP1
     &          , NFlux_VSinkD1, NFlux_VSinkD2
     &          , NumVSinkTerms

#  elif defined BIO_BioEBUS
     &          , iNO3_, iNO2_, iNH4_, iPhy1, iPhy2, iZoo1, iZoo2
     &          , iDet1, iDet2, iDON, iO2
#   ifdef NITROUS_OXIDE
     &          , iN2O
#   endif
     &          , NFlux_lightlimitP1, NFlux_lightlimitP2
     &          , NFlux_templimitP1, NFlux_templimitP2
     &          , NFlux_NO3limitP1, NFlux_NO2limitP1
     &          , NFlux_NH4limitP1, NFlux_NO3limitP2
     &          , NFlux_NO2limitP2, NFlux_NH4limitP2
     &          , NFlux_ProdNO3P1, NFlux_ProdNO3P2
     &          , NFlux_ProdNO2P1, NFlux_ProdNO2P2
     &          , NFlux_Nitrif1, NFlux_Nitrif2, NFlux_ProdNH4P1
     &          , NFlux_ProdNH4P2, NFlux_P1Z1Grazing
     &          , NFlux_P2Z1Grazing, NFlux_P1mort, NFlux_P2mort
     &          , NFlux_P1Z2Grazing, NFlux_P2Z2Grazing
     &          , NFlux_Z1Z2Grazing, NFlux_Z1metab, NFlux_Z1mort
     &          , NFlux_Z2metab, NFlux_Z2mort, NFlux_HydrolD1
     &          , NFlux_ReminOxyD1, NFlux_Denitr1D1
     &          , NFlux_Denitr2D1
     &          , NFlux_HydrolD2, NFlux_ReminOxyD2
     &          , NFlux_Denitr1D2, NFlux_Denitr2D2
     &          , NFlux_ReminOxyDON
     &          , NFlux_Denitr1DON, NFlux_Denitr2DON
     &          , NFlux_NO2anammox
     &          , NFlux_NH4anammox, O2Flux_GasExc, NumFluxTermsN
#   ifdef NITROUS_OXIDE
     &          , NFlux_paramN2O, N2OFlux_GasExc
#   endif
     &          , NumFluxTerms, NumGasExcTerms
     &          , NFlux_VSinkP2, NFlux_VSinkD1
     &          , NFlux_VSinkD2, NumVSinkTerms
#  endif
# endif   /* BIOLOGY */

# ifdef SEDIMENT
     &          ,isand1,imud1,isand2,imud2,igrav1,igrav2
# endif

!
! ================  Parameters  =====================
!

# ifdef SALINITY
      parameter (isalt=itemp+1)
# endif
# ifdef PASSIVE_TRACER
      parameter (itpas=itemp+ntrc_salt+1)
# endif

!
! ===  BIOLOGY  ===
!
# ifdef BIOLOGY
#  ifdef PISCES
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+1)
      parameter (iDIC_=itrc_bio, iTAL_=iDIC_+1, iOXY_=iDIC_+2)
#   ifdef key_pisces_npzd
      parameter ( iPOC_=iDIC_+3,  iPHY_=iDIC_+4, iZOO_=iDIC_+5,
     &            iDOC_=iDIC_+6,  iNO3_=iDIC_+7, iFER_=iDIC_+8)
#   endif
#   if ! defined key_pisces_npzd
      parameter ( iCAL_=iDIC_+3,  iPO4_=iDIC_+4,
     &            iPOC_=iDIC_+5,  iSIL_=iDIC_+6,  iPHY_=iDIC_+7,
     &            iZOO_=iDIC_+8,  iDOC_=iDIC_+9,  iDIA_=iDIC_+10,
     &            iMES_=iDIC_+11, iDSI_=iDIC_+12, iFER_=iDIC_+13,
     &            iBFE_=iDIC_+14, iGOC_=iDIC_+15, iSFE_=iDIC_+16,
     &            iDFE_=iDIC_+17, iGSI_=iDIC_+18, iNFE_=iDIC_+19,
     &            iNCH_=iDIC_+20, iDCH_=iDIC_+21, iNO3_=iDIC_+22,
     &            iNH4_=iDIC_+23)
#    ifdef key_ligand
      parameter (iLGW_=iDIC_+24)
#    endif
#   endif
#   ifdef key_pisces_quota
#    ifdef key_ligand
      parameter (iDON_=iDIC_+25, iDOP_=iDIC_+26, iPON_=iDIC_+27,
     &	         iPOP_=iDIC_+28, iNPH_=iDIC_+29, iPPH_=iDIC_+30,
     &	         iNDI_=iDIC_+31, iPDI_=iDIC_+32, iPIC_=iDIC_+33,
     &	         iNPI_=iDIC_+34, iPPI_=iDIC_+35, iPFE_=iDIC_+36,
     &	         iPCH_=iDIC_+37, iGON_=iDIC_+38, iGOP_=iDIC_+39)
#    else
      parameter (iDON_=iDIC_+24, iDOP_=iDIC_+25, iPON_=iDIC_+26,
     &           iPOP_=iDIC_+27, iNPH_=iDIC_+28, iPPH_=iDIC_+29,
     &           iNDI_=iDIC_+30, iPDI_=iDIC_+31, iPIC_=iDIC_+32,
     &           iNPI_=iDIC_+33, iPPI_=iDIC_+34, iPFE_=iDIC_+35,
     &           iPCH_=iDIC_+36, iGON_=iDIC_+37, iGOP_=iDIC_+38)
#    endif
#   endif
#   ifdef key_trc_diaadd
      parameter (Nhi        = 1,
     &            Nco3      = 2,
     &            Naksp     = 3,
     &            Netot     = 4,
     &            Nprorca   = 5,
     &            Ngrapoc   = 6,
     &            Nmico2    = 7,
     &            Nremino2  = 8,
     &            Nfixo2    = 9,
     &            Nirondep  = 10,
     &            Nironsed  = 11,
     &            Npronewo2 = 12,
     &            Npronew   = 13,
#    if defined key_pisces_npzd
     &            NumFluxTerms = Npronew)
#    else
     &            Npronewd  = 14,
     &            Nprorcad  = 15,
     &            Nprobsi   = 16,
     &            Nprofen   = 17,
     &            Nprofed   = 18,
     &            Nprorego2 = 19,
     &            Ngrapoc2  = 20,
     &            Nmeso2    = 21,
     &            Nnitrifo2 = 22,
     &            NumFluxTerms = Nnitrifo2)
#    endif

       parameter (Nfld      = 1,
     &            Nflu16    = 2,
     &            Nkgco2    = 3,
     &            Natcco2   = 4,
     &            Nsinking  = 5,
     &            Nheup     = 6,
     &            Nno3dep   = 7,
     &            Nnitrpot  = 8,
#    if defined key_pisces_npzd
     &            NumGasExcTerms = 0,
     &            NumVSinkTerms = Nnitrpot)
#    else
     &            Nsinkfer  = 9,
     &            Nsinksil  = 10,
     &            Nsinkcal  = 11,
     &            Nsildep   = 12,
     &            Npo4dep   = 13,
     &            Nnh4dep   = 14,
     &            NumGasExcTerms = 0,
     &            NumVSinkTerms = Nnh4dep)
#    endif
#   else
       parameter (NumFluxTerms = 0)
       parameter (NumGasExcTerms = 0, NumVSinkTerms = 0)
#   endif

#  elif defined BIO_NChlPZD
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+1)
      parameter (iNO3_=itrc_bio, iChla=iNO3_+1,
     &           iPhy1=iNO3_+2,
     &           iZoo1=iNO3_+3,
     &           iDet1=iNO3_+4)
#   ifdef OXYGEN
      parameter (iO2=iNO3_+5)
#   endif
      parameter (NFlux_NewProd  = 1,
     &           NFlux_Grazing  = 2,
     &           NFlux_SlopFeed = 3,
     &           NFlux_Pmort    = 4,
     &           NFlux_Zmetab   = 5,
     &           NFlux_Zmort    = 6,
     &           NFlux_ReminD   = 7,
     &           NumFluxTermsN  = 7,
#   ifdef OXYGEN
     &           OGain_NewProd  = 8,
     &           OLoss_Zmetab   = 9,
     &           OLoss_ReminD   = 10,
     &           NumFluxTermsO  = 3,
     &           OFlux_GasExc   = 1,
     &           NumGasExcTerms = 1,
     &           NumFluxTerms   = 10,
#   else
     &           NumGasExcTerms = 0,
     &           NumFluxTerms   = 7,
#   endif
     &           NFlux_VSinkP1  = 1,
     &           NFlux_VSinkD1  = 2,
     &           NumVSinkTerms  = 2)

#  elif defined BIO_N2ChlPZD2
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+1)
      parameter (iNO3_=itrc_bio, iNH4_=iNO3_+1, iChla=iNO3_+2,
     &           iPhy1=iNO3_+3,
     &           iZoo1=iNO3_+4,
     &           iDet1=iNO3_+5, iDet2=iNO3_+6)
      parameter (NFlux_NewProd  = 1,
     &           NFlux_RegProd  = 2,
     &           NFlux_Nitrific = 3,
     &           NFlux_Grazing  = 4,
     &           NFlux_SlopFeed = 5,
     &           NFlux_Pmort    = 6,
     &           NFlux_Zmetab   = 7,
     &           NFlux_Zmort    = 8,
     &           NFlux_CoagPhy  = 9,
     &           NFlux_CoagSDet = 10,
     &           NFlux_ReminD1  = 11,
     &           NFlux_ReminD2  = 12,
     &           NumFluxTermsN  = 12,
     &           NumFluxTerms   = 12,
     &           NumGasExcTerms = 0,
     &           NFlux_VSinkP1  = 1,
     &           NFlux_VSinkD1  = 2,
     &           NFlux_VSinkD2  = 3,
     &           NumVSinkTerms  = 3)

#  elif defined BIO_BioEBUS
#   ifdef NITROUS_OXIDE
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+1)
#   else
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+1)
#   endif

      parameter (iNO3_=itrc_bio, iNO2_=iNO3_+1, iNH4_=iNO3_+2,
     &           iPhy1=iNO3_+3,  iPhy2=iNO3_+4,
     &           iZoo1=iNO3_+5,  iZoo2=iNO3_+6,
     &           iDet1=iNO3_+7,  iDet2=iNO3_+8,
     &		 iDON=iNO3_+9,   iO2=iNO3_+10)

#   ifdef NITROUS_OXIDE
      parameter (iN2O=iNO3_+11)
#   endif

      parameter(  NFlux_lightlimitP1=1
     &          , NFlux_lightlimitP2=2
     &          , NFlux_templimitP1=3
     &          , NFlux_templimitP2=4
     &          , NFlux_NO3limitP1=5
     &          , NFlux_NO2limitP1=6
     &          , NFlux_NH4limitP1=7
     &          , NFlux_NO3limitP2=8
     &          , NFlux_NO2limitP2=9
     &          , NFlux_NH4limitP2=10
     &          , NFlux_ProdNO3P1=11
     &          , NFlux_ProdNO3P2=12
     &          , NFlux_ProdNO2P1=13
     &          , NFlux_ProdNO2P2=14
     &          , NFlux_Nitrif1=15
     &          , NFlux_Nitrif2=16
     &          , NFlux_ProdNH4P1=17
     &          , NFlux_ProdNH4P2=18
     &          , NFlux_P1Z1Grazing=19
     &          , NFlux_P2Z1Grazing=20
     &          , NFlux_P1mort=21
     &          , NFlux_P2mort=22
     &          , NFlux_P1Z2Grazing=23
     &          , NFlux_P2Z2Grazing=24
     &          , NFlux_Z1Z2Grazing=25
     &          , NFlux_Z1metab=26
     &          , NFlux_Z1mort=27
     &          , NFlux_Z2metab=28
     &          , NFlux_Z2mort=29
     &          , NFlux_HydrolD1=30
     &          , NFlux_ReminOxyD1=31
     &          , NFlux_Denitr1D1=32
     &          , NFlux_Denitr2D1=33
     &          , NFlux_HydrolD2=34
     &          , NFlux_ReminOxyD2=35
     &          , NFlux_Denitr1D2=36
     &          , NFlux_Denitr2D2=37
     &          , NFlux_ReminOxyDON=38
     &          , NFlux_Denitr1DON=39
     &          , NFlux_Denitr2DON=40
     &          , NFlux_NO2anammox=41
     &          , NFlux_NH4anammox=42
#   ifdef NITROUS_OXIDE
     &          , NFlux_paramN2O=43
     &          , NumFluxTermsN=NFlux_paramN2O
#   else
     &          , NumFluxTermsN=NFlux_NH4anammox
#   endif
     &          , NumFluxTerms=NumFluxTermsN
     &          , O2Flux_GasExc=1
#   ifdef NITROUS_OXIDE
     &          , N2OFlux_GasExc=2
     &          , NumGasExcTerms=2
#   else
     &          , NumGasExcTerms=1
#   endif
     &          , NFlux_VSinkP2=1
     &          , NFlux_VSinkD1=2
     &          , NFlux_VSinkD2=3
     &          , NumVSinkTerms=3)
#  endif  /* PISCES ... */

!
! ===  BIOLOGY DIAGS ===
!

#  if defined BIO_NChlPZD || defined BIO_N2ChlPZD2 || defined PISCES \
                          || defined BIO_BioEBUS
#   ifdef DIAGNOSTICS_BIO
      parameter (ntrc_diabio=NumFluxTerms+
     &                       NumGasExcTerms+NumVSinkTerms)
#   else
      parameter (ntrc_diabio=0)
#   endif
#  else
      parameter (ntrc_diabio=0)
#  endif
# else
      parameter (ntrc_diabio=0)
# endif /* BIOLOGY */

!
! Total number of active tracers
!
      integer, parameter :: NTA = itemp+ntrc_salt
!
! Total number of tracers
!
# ifdef SUBSTANCE
      integer, parameter :: NT = itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_subs
      integer, parameter :: NTot = NT+ntfix
      integer, parameter :: itsubs1 = itemp+ntrc_salt+ntrc_pas+ntrc_bio+1
      integer, parameter :: itsubs2 = itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_subs
# else
      integer, parameter :: NT = itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed
      integer, parameter :: NTot = NT
# endif /* SUBSTANCE */


# ifdef SEDIMENT
      parameter (itrc_sed=itemp+ntrc_salt+ntrc_pas+ntrc_bio+1)
      parameter (itrc_sand=itrc_sed,itrc_mud=itrc_sand+NSAND)
      parameter (itrc_grav=itrc_mud+NGRAV)
      parameter (isand1=1,imud1=isand1+NSAND)
      parameter (igrav1=imud1+NMUD)
      parameter (isand2=isand1+NSAND-1,imud2=imud1+NMUD-1)
      parameter (igrav2=igrav1+NGRAV-1)
# endif

# ifdef MUSTANG
! Vertical dimension (ksdmin:ksdmax) of variables in sediment
      integer, parameter :: ksdmin = 1
# endif /* MUSTANG */
!
! ===  u,v and tracer equations Diagnostics  ===
# ifdef DIAGNOSTICS_TS
#  ifdef DIAGNOSTICS_TS_MLD
#   ifdef DIAGNOSTICS_TS_MLD_DENS
      parameter (ntrc_diats=18*NT)
#   else
      parameter (ntrc_diats=17*NT)
#   endif
#  else
      parameter (ntrc_diats=8*NT)
#  endif
# else
      parameter (ntrc_diats=0)
# endif
# ifdef DIAGNOSTICS_UV
      parameter (ntrc_diauv=26)
# else
      parameter (ntrc_diauv=0)
# endif
# ifdef DIAGNOSTICS_VRT
      parameter (ntrc_diavrt=16)
# else
      parameter (ntrc_diavrt=0)
# endif
# ifdef DIAGNOSTICS_KE
      parameter (ntrc_diaek=18)
# else
      parameter (ntrc_diaek=0)
# endif
# ifdef DIAGNOSTICS_PV
      parameter (ntrc_diapv=12)
# else
      parameter (ntrc_diapv=0)
# endif
# if defined DIAGNOSTICS_EDDY && ! defined XIOS
      parameter (ntrc_diaeddy=15)
# else
      parameter (ntrc_diaeddy=0)
# endif
# if defined OUTPUTS_SURFACE && ! defined XIOS
      parameter (ntrc_surf=5)
# else
      parameter (ntrc_surf=0)
# endif
#endif /*SOLVE3D */


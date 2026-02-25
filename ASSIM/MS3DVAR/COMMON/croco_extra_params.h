!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA,
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
! croco_extra_params.h - Additional CROCO parameters for MS3DVAR
! ==============================================================
!
! This file contains the CROCO parameters from OCEAN/param.h that are
! NOT included in param_ms3dvar.h (which covers grid/tiling only).
!
! Specifically, this provides:
! - Vtransform (sigma coordinate type)
! - Tracer counts (NT, ntrc_temp, ntrc_salt, ntrc_bio, etc.)
! - Tracer identification indices (isalt, itrc_bio, etc.)
! - Diagnostic tracer counts (ntrc_diats, ntrc_diauv, etc.)
! - BSTRESS_FAST parameter
!
! This file is included at the end of param_ms3dvar.h so all MS3DVAR
! variants automatically get the full CROCO parameter set.
!
! SOURCE: Extracted from OCEAN/param.h lines 427-1089
!         (everything after the grid/tiling section)
!
!======================================================================

!
!----------------------------------------------------------------------
! I/O : flag for type sigma vertical transformation
!----------------------------------------------------------------------
!
#ifdef NEW_S_COORD
      real Vtransform
      parameter (Vtransform=2)
#else
      real Vtransform
      parameter (Vtransform=1)
#endif

!
!----------------------------------------------------------------------
! Number of tracers
!----------------------------------------------------------------------
!
#ifdef SOLVE3D
      integer   NT, NTA, itemp, NTot
      integer   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer   ntrc_subs, ntrc_substot
      integer   ntrc_mld
!
# ifdef TEMPERATURE
      parameter (itemp=1)
      parameter (ntrc_temp=1)
# else
      parameter (itemp=0)
      parameter (ntrc_temp=0)
# endif
# ifdef SALINITY
      parameter (ntrc_salt=1)
# else
      parameter (ntrc_salt=0)
# endif
#if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT
      parameter (ntrc_mld=3)
# else
      parameter (ntrc_mld=0)
# endif
# ifdef PASSIVE_TRACER
#  ifdef KH_INST
      parameter (ntrc_pas=2)
#  else
      parameter (ntrc_pas=1)
#  endif
# else
      parameter (ntrc_pas=0)
# endif
# ifdef BIOLOGY
#  ifdef PISCES
#   ifdef key_pisces_npzd
         parameter (ntrc_bio=9)
#   elif defined key_pisces_quota
#    ifdef key_ligand
         parameter (ntrc_bio=40)
#    else
         parameter (ntrc_bio=39)
#    endif
#   else
#    ifdef key_ligand
         parameter (ntrc_bio=25)
#    else
         parameter (ntrc_bio=24)
#    endif
#   endif
#  elif defined BIO_NChlPZD
#   ifdef OXYGEN
      parameter (ntrc_bio=6)
#   else
      parameter (ntrc_bio=5)
#   endif
#  elif defined BIO_N2ChlPZD2
      parameter (ntrc_bio=7)
#  elif defined BIO_BioEBUS
#   ifdef NITROUS_OXIDE
      parameter (ntrc_bio=12)
#   else
      parameter (ntrc_bio=11)
#   endif
#  endif
# else
      parameter (ntrc_bio=0)
# endif /* BIOLOGY */

/*! === SUBSTANCE ===*/
!
# if defined SUBSTANCE
! ntrc_subs : number of advected substances (not fixed, neither benthic)
      integer  itsubs1, itsubs2, ntfix
#  ifdef SED_TOY
#   if defined SED_TOY_FLOC_0D || defined SED_TOY_FLOC_1D
      parameter (ntrc_subs=15 , ntfix=0, ntrc_substot=ntrc_subs+ntfix )
#   else
      parameter (ntrc_subs=6 , ntfix=0, ntrc_substot=ntrc_subs+ntfix )
#   endif
#  elif defined TIDAL_FLAT
      parameter (ntrc_subs=3 , ntfix=0, ntrc_substot=ntrc_subs+ntfix )
#  elif defined ESTUARY
      parameter (ntrc_subs=2 , ntfix=0, ntrc_substot=ntrc_subs+ntfix )
#  elif defined VILAINE
      parameter (ntrc_subs=3 , ntfix=0, ntrc_substot=ntrc_subs+ntfix )
#  else
      parameter (ntrc_subs=2 , ntfix=0, ntrc_substot=ntrc_subs+ntfix )
#  endif
      parameter (itsubs1= itemp+ntrc_salt+ntrc_pas+ntrc_mld+ntrc_bio+1 )
      parameter (itsubs2= itemp+ntrc_salt+ntrc_pas+ntrc_mld+ntrc_bio+ntrc_subs )
# else
      parameter (ntrc_subs=0, ntrc_substot=0)
# endif /* SUBSTANCE */

!
# ifdef SEDIMENT
! NSAND          Number of sand classes
! NMUD           Number of mud classes
! NGRAV          Number of gravel classes (not implemented...)
! NST            Number of sediment (tracer) size classes
! NLAY           Number of layers in sediment bed
!
      integer NSAND, NMUD, NGRAV, NST, NLAY
#  ifdef DUNE
#   ifdef ANA_DUNE
      parameter (NSAND=1, NMUD=0, NGRAV=0)
      parameter (NLAY=11)
#   else
      parameter (NSAND=2, NMUD=0, NGRAV=0)
      parameter (NLAY=10)
#   endif
#  elif defined SED_TOY
#   if defined SED_TOY_RESUSP || defined SED_TOY_CONSOLID
      parameter (NSAND=2, NMUD=2, NGRAV=0)
      parameter (NLAY=41)
#   elif defined SED_TOY_FLOC_0D || defined SED_TOY_FLOC_1D
      parameter (NSAND=4, NMUD=15, NGRAV=0)
      parameter (NLAY=20)
#   elif defined SED_TOY_ROUSE
      parameter (NSAND=0, NMUD=6, NGRAV=0)
      parameter (NLAY=1)
#   endif
#  else
      parameter (NSAND=2, NMUD=0, NGRAV=0)
      parameter (NLAY=1)
#  endif
      parameter (ntrc_sed=NSAND+NMUD+NGRAV)
      parameter (NST=ntrc_sed)
# else
      parameter (ntrc_sed=0)
# endif /* SEDIMENT */
!
! Total number of active tracers
!
      parameter (NTA=itemp+ntrc_salt)

!
! Total number of tracers
!
# ifdef SUBSTANCE
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_subs+ntrc_mld)
      parameter (NTot=NT+ntfix)
# else
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_mld)
      parameter (NTot=NT)
# endif /* SUBSTANCE */

# ifdef MUSTANG
   ! vertical dimension (ksdmin:ksdmax) of variables in sediment
   ! (ksdmax=max number of layers)
      integer ksdmin,ksdmax
#  if defined ANA_DUNE || defined key_ANA_bedload
      parameter (ksdmin=1,ksdmax=11)
#  elif defined TIDAL_FLAT || defined ESTUARY
      parameter (ksdmin=1,ksdmax=3)
#  else
      parameter (ksdmin=1,ksdmax=10)
#  endif
# endif /* MUSTANG */

# if defined SEDIMENT && defined AGRIF
      integer Agrif_lev_sedim
      parameter (Agrif_lev_sedim=0)
# endif

# ifdef GLS_MIXING
      integer NGLS
      parameter(NGLS=2)
      integer itke
      parameter(itke=1)
      integer igls
      parameter(igls=2)
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
# if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT
     &          , iCRT2, iCRT3, iCRT4
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
# if defined DIAGNOSTICS_TS_MLD && defined DIAGNOSTICS_TS_MLD_CRIT
      parameter (iCRT2=itemp+ntrc_salt+1)
      parameter (iCRT3=itemp+ntrc_salt+2)
      parameter (iCRT4=itemp+ntrc_salt+3)
# endif
# ifdef PASSIVE_TRACER
      parameter (itpas=itemp+ntrc_salt+ntrc_mld+1)
# endif

!
! ===  BIOLOGY  ===
!
# ifdef BIOLOGY
#  ifdef PISCES
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+ntrc_mld+1)
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
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+ntrc_mld+1)
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
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+ntrc_mld+1)
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
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+ntrc_mld+1)
#   else
      parameter (itrc_bio=itemp+ntrc_salt+ntrc_pas+ntrc_mld+1)
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
! === SEDIMENTS ===
!

# ifdef SEDIMENT
      parameter (itrc_sed=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_mld+1)
      parameter (itrc_sand=itrc_sed,itrc_mud=itrc_sand+NSAND)
      parameter (itrc_grav=itrc_mud+NGRAV)
      parameter (isand1=1,imud1=isand1+NSAND)
      parameter (igrav1=imud1+NMUD)
      parameter (isand2=isand1+NSAND-1,imud2=imud1+NMUD-1)
      parameter (igrav2=igrav1+NGRAV-1)
# endif

!
! ===  u,v and tracer equations Diagnostics  ===
!
# ifdef DIAGNOSTICS_TS
#  ifdef DIAGNOSTICS_TS_MLD
      parameter (ntrc_diats=19*NT)
#  else
      parameter (ntrc_diats=8*NT)
#  endif
# else
      parameter (ntrc_diats=0)
# endif
# ifdef DIAGNOSTICS_UV
      parameter (ntrc_diauv=24)
# else
      parameter (ntrc_diauv=0)
# endif
# ifdef DIAGNOSTICS_VRT
      parameter (ntrc_diavrt=16)
# else
      parameter (ntrc_diavrt=0)
# endif
# ifdef DIAGNOSTICS_EK
#  ifdef DIAGNOSTICS_EK_MLD
      parameter (ntrc_diaek=28)
#  else
      parameter (ntrc_diaek=16)
#  endif
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

!
!----------------------------------------------------------------------
! Max time increment for computing bottom stress at the 3D fast time
! steps
!----------------------------------------------------------------------
!
#ifdef BSTRESS_FAST
      integer inc_faststep_max
      parameter(inc_faststep_max = 10)
#endif

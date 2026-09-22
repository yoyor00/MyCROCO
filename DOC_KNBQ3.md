# Documentation de la fonctionnalité KNBQ3

## 1. Vue d'ensemble

La fonctionnalité **KNBQ3** (Keyed NBQ-3) est un **solveur de vitesse non-hydrostatique compressible** pour le modèle de circulation océanique CROCO (CROCO is a branch of ROMS). Il implémente une décomposition en modes rapide/lent (fast/slow mode splitting) pour résoudre efficacement les équations de Navier-Stokes compressibles en océanographie.

**Activer le solveur** : via la clé de préprocesseur `K3FAST` dans `cppdefs.h`.

## 2. Principes physiques

### 2.1 Décomposition en modes

Le solveur sépare la dynamique en :
- **Mode externe (2D)** : mouvement de surface (zeta), courant moyen barotrope (ubar, vbar) — se propage rapidement.
- **Mode interne (3D)** : composantes rapides des flux de quantité de mouvement (`qdmu_nbq`, `qdmv_nbq`, `qdmw_nbq`) et de la densité (`rho_nbq`).

### 2.2 Variables principales (définies dans `nbq.h`)

| Variable | Description |
|---|---|
| `qdmu_nbq`, `qdmv_nbq`, `qdmw_nbq` | Flux de quantité de mouvement rapides (U, V, W) |
| `rho_nbq` | Anomalie de densité rapide (compressibilité) |
| `thetadiv_nbq` | Terme "pression-viscosité" θ ou divergence du mouvement |
| `ru_nbq`, `rv_nbq`, `rw_nbq` | Intégrales verticales des forçages RHS |
| `ru_int_nbq`, `rv_int_nbq`, `rw_int_nbq` | RHS couplés lent→rapide |
| `soundspeed_nbq`, `soundspeed2_nbq` | Vitesse du son (compressibilité) |
| `visc2_nbq` | Coefficient de viscosité de second ordre |
| `Hzw_nbq`, `zw_nbq` | Épaisseur des couches et positions w |
| `dzeta_nbq`, `wsurf_nbq`, `usurf_nbq`, `vsurf_nbq` | Cinématique de surface |
| `CFL_nbq`, `dtnbq` | Rapport CFL et pas de temps du mode rapide |

## 3. Architecture du solveur

Le cœur du solveur est le fichier `step3d_fast.F` qui décompose le pas de temps rapide en **5 parties** :

### PART B : Traitement du mode externe 2D
1. **B.0** — Préparation de l'intégration 2D (`k3fast_2Dmode_prep.h`)
2. **B.1** — Mise à jour AGRIF pour le nesting (`k3fast_AGRIF0.h`)
3. **B.2** — Calcul de `zeta(m+1)` par la relation cinématique de surface (`k3fast_zeta_update.h`)
4. **B.3** — Mise à jour de la grille verticale (`set_depth_tile`, `grid_nbq_tile`)
5. **B.4** — Pas arrière AM4 pour les RHS 2D : gradient de pression de surface, advection, Coriolis, frottement (`k3fast_2Dmode_rhs.h`)
6. **B.5** — Pré-étape 2D avec couplage (`k3fast_pre_step2d.h`)

### PART C : Résolution des équations compressibles 3D
1. **C.0** — Initialisations et frottement de fond rapide (`k3fast_fbf.h`, `k3fast_init.h`)
2. **C.1** — Force de Coriolis non-traditionnelle (`k3fast_ntcoriolis.h`)
3. **C.2** — Calcul de θ = pression-compressibilité + viscosité de second ordre (`k3fast_qdmuv_update.h`)
4. **C.3** — Mise à jour des moments U et V (`k3fast_qdmuv_update.h`)
5. **C.4** — Mise à jour du moment W :
   - Schéma explicite direct, ou
   - Schéma implicite avec élimination de Gauss tridiagonale (`k3fast_qdmw_implicit.h`)
6. **C.5** — Opérateur de divergence horizontal + vertical (`k3fast_divh.h`, `k3fast_divv.h`)
7. **C.6** — Conservation de la masse : `rho_nbq(m+1) = rho_nbq(m) - dtfast * DIV(m+1)` (`k3fast_mass_update.h`)

### PART D : Post-traitement du mode rapide
1. **D.1** — Moyenne temporelle et ajustement des variables barotropes (`k3fast_post.h`)
2. **D.2** — Ajustement final de `zeta` par conservation de masse intégrée (`k3fast_zeta_correct.h`)
3. **D.3** — Conservation AGRIF pour le nesting (`k3fast_AGRIF1.h`)
4. **D.4** — Correction de `Hz` par inversion de la continuité interne (`k3fast_hz_correct.h`)

### PART E : Test CFL
Vérification que `VMAX < 100 m/s` et détection de NaN pour arrêter proprement la simulation.

## 4. Schémas de intégration temporelle

Le solveur utilise l'algorithme **Généralisé Forward-Backward AB3-AM4** (Shchepetkin & McWilliams, 2009) :
- **Forward (explicite)** : zeta, qdmu_nbq, qdmw_nbq
- **Backward (implicite)** : rho_nbq (et qdmw_nbq si `NBQ_IMP`)

Les paramètres AB3-AM4 sont :
```
myalpha   = 0.1 (0.01 pour TANK/AgAc)
myepsilon = 0.00976186 - 0.13451357*myalpha
mygamma   = 0.08344500 - 0.51358400*myalpha
mybeta    = 0.281105
```

Des options supplémentaires `K3FAST_AM4`, `K3FAST_AM4b`, `K3FAST_AM4d` activent des variantes du schéma.

## 5. Intégration W-momentum (mode lent K3SLOW_W)

Le fichier `step3d_w.F` gère la mise à jour de `wz` (vitesse verticale) avec :
- **Diffusion verticale turbulente implicite** (tridiagonale) quand `NBQ_WDIF` est actif
- Conditions aux limites :
  - Fond : `wz=0` (no-slip) ou `NBQ_FREESLIP` (glissement libre via relation cinématique)
  - Surface : relation cinématique
- Le RHS est calculé dans `rhs3d_w_nh.F` avec :
  - Advection horizontale (C2, TVD, WENO5, C6, UP5)
  - Advection verticale (C2, TVD, SPLINES, WENO5, C6, UP5)

## 6. Options (clés de préprocesseur)

### 6.1 Compressibilité et vitesse du son
| Clé | Description |
|---|---|
| `K3FAST_RHO` | Active la densité rapide compressible |
| `K3FAST_CSVISC2K` | Visco-sound speed 2D → 3D |
| `NBQ_RCSOUND` | Vitesse du son recomputée |
| `K3FAST_NOBPG` | Pas de gradient de pression de fond |

### 6.2 Couplage lent/rapide
| Clé | Description |
|---|---|
| `K3FAST_C3D_UVSF` | Couplage 3D complet U/V |
| `K3FAST_C3D_WSF` | Couplage 3D complet W |
| `K3FAST_COUPLING_SCH0/1/2` | Schéma de couplage temporal |
| `K3FAST_COUPLING2D`, `K3FAST_COUPLING3D` | Couplage 2D/3D |
| `K3FAST_COUPLINGW_SCH1/2` | Couplage W |
| `K3FAST_DISSCOUPLING` | Dissipation de couplage (`alpha_uv`, `alpha_w`) |

### 6.3 Schémas numériques
| Clé | Description |
|---|---|
| `NBQ_IMP` | Intégration implicite de W (tridiagonale) |
| `NBQ_THETAIMP` | θ implicite pour viscosité (θimp_nbq = 0.6 ou 1.0) |
| `K3FAST_AB3` | Schéma Adams-Bashforth 3 |
| `K3FAST_UV`, `K3FAST_W` | Activation U/V ou W rapide |
| `K3FAST_PG2` | Schéma de gradient de pression NBQ |

### 6.4 Grille et topographie
| Clé | Description |
|---|---|
| `NBQ_GRID_SLOW` | Mise à jour grille à fréquence réduite |
| `NBQ_HZCORRECT` | Correction de Hz par continuité |
| `NBQ_HZCORR_DEBUG` | Debug correction Hz |
| `NBQ_HZ_PROGNOSTIC` | Hz prognostique (couches sédimentaires) |
| `K3FAST_SEDLAYERS` | Couches sédimentaires |
| `CUVE_BATHY`, `BATHY_SLOPE` | Bathymétrie spécifique |

### 6.5 Conditions aux limites et autres
| Clé | Description |
|---|---|
| `NBQ_FREESLIP` | Glissement libre aux fond/surface |
| `OBC_NBQ` | Conditions aux limites NBQ |
| `NBQ_NUDGING`, `NBQ_NUDGING_W` | Nudging |
| `KNHINT_CORR` | Coriolis/nuising dans la couche de fond |
| `KNHINT_3M` | Intégration conditionnelle 3D mode |
| `WET_DRY` | Mouillage/assèchement |
| `MVB` | Bathymétrie mobile |
| `NBQ_SPONGE` | Couche amortissante |
| `NBQ_MASS` | Conservation de masse |
| `NBQ_GRAV` | Gravité |
| `K3FAST_SACOUS` | Acoustique sous-marine |
| `K3FAST_DIAGACOUS` | Diagnostics acoustiques |
| `CENTRIFUGE` | Force centrifuge |
| `PSOURCE` | Sources ponctuelles |

## 7. Intégration avec AGRIF (nesting multi-résolution)

Les fichiers `update2D.F` et `update3D.F` gèrent la mise à jour bidirectionnelle entre grille fine (KNBQ3) et grille parente (hydrostatique classique) :

- **update2D** : `Updateunbq`, `Updatevnbq`, `Updatewnbq`, `Updaterhonbq`, `Updateubar`, `Updatevbar`, `Updatezeta`
- **update3D** : `Agrif_update_np1` (traceurs), `Agrif_update_uv_np1` (U, V, W)
- Conservation du volume : `AGRIF_CONSERV_VOL`

## 8. Structure des fichiers

```
KNBQ3/
├── Makefile                         # Compilation
├── nbq.h                            # Variables communes du solveur
├── step3d_fast.F                    # Cœur du solveur rapide 3D
├── step3d_w.F                       # Mise à jour W (mode lent)
├── rhs3d_w_nh.F                     # RHS vertical velocity
├── pre_step3d_KNBQ3.h               # Pré-étape : prédicteurs U/V/W
├── update2D.F                       # Nesting AGRIF 2D
├── update3D.F                       # Nesting AGRIF 3D
├── grid_nbq.F                       # Grille NBQ
├── initial_nbq.F                    # Initialisation
├── nbq_bry_store.F                  # Stockage BC NBQ
├── zoom.F, zoom.h                   # Support AGRIF zoom
├── zoombc_3D.F, zoombc_3Dfast.F     # BC zoom 3D
├── w3dbc.F                          # BC verticales
├── unbq_bc.F, vnbq_bc.F, wnbq_bc.F, rnbq_bc.F  # BC NBQ
├── rnbq_bc.F                        # BC densité
├── prepro_KXX.py                    # Pré-traitement KXX AGRIF
├── common2device.exe                # Utilitaire OpenACC
│
├── k3fast_*.h                       # Modules du solveur rapide
│   ├── k3fast_2Dmode_prep.h         # Prep mode 2D
│   ├── k3fast_2Dmode_rhs.h          # RHS 2D (AM4 backward)
│   ├── k3fast_zeta_update.h         # Mise à jour zeta
│   ├── k3fast_zeta_correct.h        # Correction zeta
│   ├── k3fast_qdmuv_update.h        # Mise à jour U/V rapide
│   ├── k3fast_qdmw_update.h         # Mise à jour W rapide
│   ├── k3fast_qdmw_implicit.h       # W implicite (tridiag)
│   ├── k3fast_mass_update.h         # Conservation masse
│   ├── k3fast_divh.h                # Divergence horizontale
│   ├── k3fast_divv.h                # Divergence verticale
│   ├── k3fast_ntcoriolis.h          # Coriolis non-traditionnel
│   ├── k3fast_fbf.h                 # Frottement fond rapide
│   ├── k3fast_init.h                # Initialisations
│   ├── k3fast_post.h                # Post-traitement
│   ├── k3fast_hz_correct.h          # Correction Hz
│   ├── k3fast_AGRIF0.h, k3fast_AGRIF1.h  # Support AGRIF
│   ├── k3fast_sacous.h              # Acoustique
│   └── k3fast_diagacousOA.h         # Diag acoustique OpenAnalysis
│
├── MY_XIOS/                         # Config XIOS personnalisée
│   ├── iodef.xml
│   ├── context_myxios.xml
│   ├── domain_def_myxios.xml
│   ├── grid_def_myxios.xml
│   ├── field_def_myxios.xml
│   ├── file_def_myxios.xml
│   └── axis_def_myxios.xml
│
└── XIOS_XMLFILES/                   # Config XIOS standard sphérique
    ├── iodef.xml
    ├── context_spherical.xml
    ├── domain_def_spherical.xml
    ├── field_def.xml
    └── file_def_spherical.xml
```

## 9. Workflow d'un pas de temps rapide KNBQ3

```
┌─────────────────────────────────────────────────────┐
│  step3d_fast (cœur du solveur)                       │
├─────────────────────────────────────────────────────┤
│  PART B : Mode externe 2D                            │
│  ├── k3fast_2Dmode_prep.h                           │
│  ├── [AGRIF0] Load RHS AGRIF                         │
│  ├── [HZCORRECT] Restore Hz                          │
│  ├── grid_nbq (variables dérivées)                   │
│  ├── k3fast_zeta_update.h  → zeta(m+1)              │
│  ├── set_depth + grid_nbq (grille)                   │
│  ├── k3fast_2Dmode_rhs.h  → rubar, rvbar (AM4)      │
│  └── k3fast_pre_step2d.h  → couplage 2D              │
├─────────────────────────────────────────────────────┤
│  PART C : Mode interne 3D compressible               │
│  ├── [FBF] Frottement fond rapide                    │
│  ├── k3fast_init.h                                   │
│  ├── [OBC_NBQ] nbq_bry_store                         │
│  ├── [KNHINT_3M] dtfast *= nsdtnbq                   │
│  ├── [UV_COR_NT] k3fast_ntcoriolis.h                 │
│  ├── k3fast_qdmuv_update.h → θ + qdmu, qdmv         │
│  ├── [NBQ_IMP] k3fast_divh.h → divergence h          │
│  ├── k3fast_qdmw_update.h → qdmw                     │
│  ├── [NBQ_IMP] k3fast_qdmw_implicit.h → tridiag W    │
│  ├── [NUDGING/CORR] Nudging W                        │
│  ├── [OBC_NBQ] wnbq_bc                               │
│  ├── [!NBQ_IMP] k3fast_divh.h                        │
│  ├── k3fast_divv.h → divergence v                     │
│  ├── [MPI] exchange divergence                        │
│  └── k3fast_mass_update.h → rho_nbq(m+1)             │
│  └── [KNHINT_3M] dtfast /= nsdtnbq                   │
├─────────────────────────────────────────────────────┤
│  PART D : Post-traitement                            │
│  ├── k3fast_post.h  → moyenne + ajustement baro      │
│  ├── k3fast_zeta_correct.h → zeta final              │
│  ├── [AGRIF] k3fast_AGRIF1.h                         │
│  └── [HZCORRECT] k3fast_hz_correct.h → Hz            │
├─────────────────────────────────────────────────────┤
│  PART E : Test CFL + détection NaN/blow-up           │
└─────────────────────────────────────────────────────┘
```

## 10. Support OpenACC (GPU)

Le code est annoté avec des directives OpenACC (`$acc kernels`, `$acc loop`, `$acc declare create`) pour l'exécution sur GPU. Les tableaux critiques sont créés sur dispositif :
- `Hzw_nbq_inv`, `Hzr_nbq_inv`
- `FX`, `FY`, `FC`, `CF`
- `dZdxq_u`, `dZdy_v`
- `ntcoru/v/w`
- `grad`

## 11. Références

- **Algorithme FB AB3-AM4** : Shchepetkin, A.F., and J.C. McWilliams, 2009 : "Computational kernel algorithms for fine-scale, multiprocess, longtime oceanic simulations." In *Handbook of Numerical Analysis: Computational Meteorology for the Atmosphere and Oceans*, R.M. Teman and J.J. Tribbia, eds, Elsevier Science.
- **CROCO** : http://www.croco-ocean.org — Branche de ROMS développée à l'IRD, INRIA, Ifremer, CNRS et Université Toulouse III.

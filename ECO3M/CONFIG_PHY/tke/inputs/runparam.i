version
263
   directory sauvegarde - adresse absolue:  cnmlog
../OUTPUTS/
   directory local de travail: cdir
../CONFIG_PHY/tke/
   titre du run: ctit
RUN
   sous titre du run: csbtit
Chla
   pas de temps en MINUTES: dt
4.
   temps initial en JOURS (par rapport a 1er Janvier a 0h = 0): tdebut0
0.    
   duree totale du run en JOURS: tfin0 1095 = 3y 1460 = 4y, 1825 = 5 y, 3650 = 10y
10.
   temps sauvegarde variables - flux en MINUTES: timesave0
360.   
flag ecriture: 0 ==> minimum ; 1 ==> ecritures evolution (nwrite)
1
   Nbr d intervalles de sauvegarde pour ecriture-check: nwrite Unites=timesave0
200
   nbr de traceurs biogeochimiques pronostiques: nbrprono
6
   nom des variables biogeochimiques pronostiques: character*3 CNMTRA
nut
phy
zoo
det
dic
alk
   nbr de traceurs biogeochimiques diagnostiques: nbrdiag 
5
   nom des variables biogeochimiques diagnostiques: character*3 CNMTRA
chl
cch
par
pur
ccq
   flag stockage des variables biogeochimiques: MSTOCK
1 1 1 1 1 1 1 1 1 1 1 
   flag stockage des tendaces diffusion prognostiques: MBIODIF
1 1 1 1 1 1
   flag stockage des tendances dues a sedimentation: mbiosed
1
   nbr de flux biogeochimiques (sauvegarde potentielle): nbrflux
16
   nom des flux biogeochimiques: character*6 CNMFLX
detnut
detzoo
nutphy
phydet
phynut
phyzoo
rapalk
rapdic
rapnut
zoodet
zoopel
zoonut
carphy
carres
alkphy
alkcal
   flag stockage des flux biogeochimiques: MFLUX
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
   Flag sauvegarde TKE: mtke
1
   Flag sauvegarde mixing length, Brunt Vaissaila, csaillement: mlength
1
   Flag sauvegarde diffision verticale traceurs: mdifft
1
   Flag sauvegarde temperature: mt
1
   Flag sauvegarde salinite et densite 0: ms
1
   Flag sauvegarde tendances equation TKE: mtket
1
   Flag sauvegarde tendances equation temperature: mtt
1
   Flag sauvegarde tendances equation salinite: mst
1
   nbr de niveaux verticaux: nzt
42
   nom fichier definition de la grille verticale (93): cgrille
inputs/dyfamed_grid.dat
C
C
C **********************************************************************
C **********************************************************************
       subroutine BIOMASSE
C
C Appele par biologie.f a chaque passage: TOUTES LES SUBROUTINES APPELEES
C    PAR biomasse.f SE TROUVENT DANS biomodel.f
C
C Calculs des termes non conservatifs: boucle principale en z.
C
C Programme en parties:
C   1. Echanges des variables de formes vectorielles a forme scalaire
C         au niveau nk = BIOVARIABLE
C   2. Calculs des differents termes non conservatifs = BIOMODEL 
C             subroutines dans biomodel.f
C   3. Calculs des rappels dynamiques = BIORAP
C   4. Integration avant et verification de la positivite des champs
C   5. Inverse de 1: retour de l'information scalaire sous forme
C         vectorielle = BIOELBAIRAV
C
C            de traceurs updates = BIOAVANT
C
C ======================================================================
C
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl

C Debut boucle en z
C 	barrier.n
C	No reason why the BIOAVANT function should be called for each z-level
C	In BIOAVANT, it is already done in 3D; i.e including the z-dimension
C       do jk=1,nzt
C          write(*,*) "+++++++++++++++++++++++++++++++jk"
C         nk=jk
C         call BIOAVANT
C       enddo

      return
      end

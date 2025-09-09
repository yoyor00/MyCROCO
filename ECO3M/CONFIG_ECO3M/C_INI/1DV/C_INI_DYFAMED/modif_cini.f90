 Implicit none
 Character(30)::fic
integer::i,j,k,l
real::conc

i=1
j=1

 write(*,*) 'fichier ?'
 read(*,*) fic
 open (15,file=fic)
 open (16,file=trim(fic)//"_OK")
 Do l=1,42
  read(15,*)k,conc
  write(16,'(I2,2X,I2,2X,I2,2X,E13.5)')i,j,k,conc
Enddo
End


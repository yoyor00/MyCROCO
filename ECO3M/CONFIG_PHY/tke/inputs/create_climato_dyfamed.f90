Program climato
IMPLICIT None
integer,dimension(1:96360)::yy,mm,dd,hh
real,dimension(1:96360)::a,c,g,h,b,d,e,f
integer,dimension(1:2920)::moy_yy,moy_mm,moy_dd,moy_hh
real,dimension(1:2920)::moy_a,moy_c,moy_g,moy_h,moy_b,moy_d,moy_e,moy_f
INTEGER::err,i,j,cpt

open(15,file="dyfamed_forcings_1d.dat")

err = 0

DO i= 1,96360
  read(15,*,IOSTAT=err) yy(i),mm(i),dd(i),hh(i),a(i),b(i),c(i),d(i),e(i),f(i),g(i),h(i)
ENDDO

Do j=1,365*8

   moy_yy(j) =0.0 
   moy_mm(j) =0.0 
   moy_dd(j) =0.0 
   moy_hh(j) =0.0 
   moy_a(j) = 0.0
   moy_b(j) = 0.0
   moy_c(j) = 0.0
   moy_d(j) = 0.0
   moy_e(j) = 0.0
   moy_f(j) = 0.0
   moy_g(j) = 0.0
   moy_h(j) = 0.0
   cpt = 0 
Do i=1,96360
     if (mm(i) == mm(j) .and. dd(i) == dd(j) .and. hh(i) ==hh(j)) then
             moy_yy(j) = moy_yy(j) + yy(i)
             moy_mm(j) = moy_mm(j) + mm(i)
             moy_dd(j) = moy_dd(j) + dd(i)
             moy_hh(j) = moy_hh(j) + hh(i)
             moy_a(j)=moy_a(j)+a(i)
             moy_b(j)=moy_b(j)+b(i)
             moy_c(j)=moy_c(j)+c(i)
             moy_d(j)=moy_d(j)+d(i)
             moy_e(j)=moy_e(j)+e(i)
             moy_f(j)=moy_f(j)+f(i)
             moy_g(j)=moy_g(j)+g(i)
             moy_h(j)=moy_h(j)+h(i)
            cpt = cpt + 1
     endif

enddo
     write(*,*) 'cpt =',cpt
enddo

Do j = 1,365*8
   moy_yy(j) = moy_yy(j)/(2012-1980+1)
   moy_mm(j) = moy_mm(j)/(2012-1980+1)
   moy_dd(j) = moy_dd(j)/(2012-1980+1)
   moy_hh(j) = moy_hh(j)/(2012-1980+1)
   moy_a(j) = moy_a(j)/(2012-1980+1) 
   moy_b(j) = moy_b(j)/(2012-1980+1)
   moy_c(j) = moy_c(j)/(2012-1980+1)
   moy_d(j) = moy_d(j)/(2012-1980+1)
   moy_e(j) = moy_e(j)/(2012-1980+1)
   moy_f(j) = moy_f(j)/(2012-1980+1)
   moy_g(j) = moy_g(j)/(2012-1980+1)
   moy_h(j) = moy_h(j)/(2012-1980+1)
Enddo

Open(16,file="dyfamed_climato_forcings_1d.dat")
Do j=1,365*8
  write(16,*) moy_yy(j),moy_mm(j),moy_dd(j),moy_hh(j),moy_a(j),moy_b(j),moy_c(j),&
 moy_d(j),moy_e(j),moy_f(j),moy_g(j),moy_h(j)
Enddo
close (15)
close(16)
 
END Program climato

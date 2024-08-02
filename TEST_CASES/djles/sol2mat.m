%this script save few interested varaible in a struct object .mat
s.u=u;
s.w=w;
s.rho=density; 
s.eta=eta; 
s.celerity = c; 
s.wavelenght=wavelength; 
s.amplitude=wave_ampl; 
s.ke=kewave; 
s.d_ke=kepert;
s.ape=apedens; 
save(sprintf("%s",name),"s")

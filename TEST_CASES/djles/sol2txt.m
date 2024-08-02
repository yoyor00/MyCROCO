%This script convert solution to .txt file croco readable % 

nx = length(x);
nz = length(z);
n = nx*nz;

%physical variables 
u_r = reshape(u, [1,n]);
w_r = reshape(w, [1,n]); 
d_r = reshape(density, [1,n]); 

%grid variables 
x_r = repmat(x, [nz,1]); 
x_grid = reshape(x_r, [1,n]); 
zj = linspace(1,length(z),nz); 
zj_r = repmat(zj,[1,nx]); 

%creating matrix of solution 
Tab = [x_grid; zj_r; u_r; w_r; d_r*1028]; 

%saving matrix 
writematrix(transpose(Tab), 'sol_test.txt', 'Delimiter', ' ')



clear all

%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = 6*10e3;%6.6*10e3; % APE for wave (m^4/s^2)
L  = 45000;%180000;%87570; %50000.0;  % domain width (m)
H  = 600;  % domain depth (m)
NX = 100;   % grid
NZ = 50;   % grid

% The unitless density profile (normalized by a reference density rho0)
a_d=1.5/1028.0; z0_d=100.0; d_d=30.0;%30.
%a_d=1.5; z0_d=200.0; d_d=30.0; yyp
rho =@(z) 1-a_d*tanh((z+z0_d)/d_d);
%% go 
rhoz =@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;
%rho0 =1028;

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;
% Find the solution of the DJL equation
djles_refine_solution

% Increase the resolution, and iterate to convergence
epsilon=1e-6;
NX=200; NZ=50;
djles_refine_solution
%sol2txt

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics
djles_plot
sol2txt
%name=sprintf('/Users/roustan/Documents/These/FRONT_SOLITON/OUTPUT/djl_H_%i_APE_%i.mat',z0_d,A);
%sol2mat



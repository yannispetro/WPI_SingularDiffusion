clc
clear all
close all

m1 = 1;
m2 = 1;
c1 = 0.1;
c2 = 0.1;
k1 = 1;
k2 = 1;

ec1 = 0.0;
ec2 = 0.0;
ek1 = 0.5;
ek2 = 0.5;

S0 = 0.1;

par = [m1,m2,c1,c2,k1,k2,ec1,ec2,ek1,ek2,S0];

ndof = 2;
ord  = 2;
h    = 7;

% Time parameters
ti = 0; % Initial time
tf = 1.0; % Final time
Fs = 100;

% Create the original grid
points = [11 11 11 11];

domain = [-1.3 1.3 -0.2 0.2 -2.1 2.1 -0.6 0.6]; % 1.0 sec   x1 x2 x1' x2'
% domain = [-2.0 2.0 -2.0 2.0 -3.0 3.0 -2.0 2.0]; % 3.0 sec   x1 x2 x1' x2'
% domain = [-2.8 2.8 -2.8 2.8 -4.0 4.0 -3.5 3.5]; % 10.0 sec   x1 x2 x1' x2'


[exponent,constraint_val,exit_flags,C_MPP,LM,elt] = a4_functional_minimization_Ritz(...
                        domain, points, ndof, ord, h, ...
                        par, ti, tf, Fs);

save(['files/WPI_t' num2str(tf) '_h' num2str(h) '_fs' num2str(Fs) ...
'_p' num2str(points(1)) '.mat'], ...
'domain','points','m1','m2','c1','c2','k1','k2','ec1','ec2','ek1','ek2','S0',...
'exponent','constraint_val','exit_flags','C_MPP','LM','elt');

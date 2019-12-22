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

M = [m1 0; 0 m2];
C = [c1+c2 -c2; -c2 c2];
K = [k1+k2 -k2; -k2 k2];

par = [m1,m2,c1,c2,k1,k2,ec1,ec2,ek1,ek2,S0];

Fs   = 100;  % Sampling frequency      
Tot  = 5;  % Total time  (sec)
% Time vector (Dimensional)
[tt, dt, nt] = time( Fs,Tot); 

% ----- Modulating function ----
At1 = ones([1,nt]); % sin(2.5*pi*tt/tt(end)).^2 + 0.2;
At2 = zeros([1,nt]); % sin(2.5*pi*tt/tt(end)).^2 + 0.2;

% nw = nt;
% cut_off = pi/dt;
% dw = cut_off/nw;
% omega = linspace(-cut_off,cut_off,nw);

tol = 1e-6;
err = 10^5;
tspan = [tt(1) tt(end)];
z0 = zeros([14,1]);
Ce = zeros([2,2,nt]);
Ke = zeros([2,2,nt]);
Ce_new = zeros([2,2,nt]);
Ke_new = zeros([2,2,nt]);
errCe = zeros([nt,1]);
errKe = zeros([nt,1]);
while err > tol
    
    sol = ode45(@(t,z) ODE_t_dep_14(t,z,par,Ce,Ke,tt,At1,At2),tspan,z0);
    z_new = deval(sol,tt);
    
    mu1  = z_new(1,:);
    mu2  = z_new(2,:);
    mu3  = z_new(3,:);
    mu4  = z_new(4,:);
    v11 = z_new(5,:);
    v22 = z_new(9,:);
    v33 = z_new(12,:);
    v44 = z_new(14,:);
    
    Ce_new(1,1,:) = 9*c1*ec1*mu3(:).^2 + 3*c1*ec1*v33(:);
    Ce_new(2,2,:) = 9*c2*ec2*mu4(:).^2 + 3*c2*ec2*v44(:);
    
    Ke_new(1,1,:) = 9*k1*ek1*mu1(:).^2 + 3*k1*ek1*v11(:);
    Ke_new(2,2,:) = 9*k2*ek2*mu2(:).^2 + 3*k2*ek2*v22(:);
    
    for i = 1:nt
        errCe(i) = ( norm(Ce(:,:,i) - Ce_new(:,:,i),2) )/norm(C);
        errKe(i) = ( norm(Ke(:,:,i) - Ke_new(:,:,i),2) )/norm(K);
    end
    errAvCe = mean(errCe)/nt;
    errAvKe = mean(errKe)/nt;
    
    err = ( errAvCe + errAvKe )/2;
    
    disp([errAvCe,errAvKe, err])
    
    Ce = Ce_new;
    Ke = Ke_new;
    
end

mean_Z = z_new(1:4,:);
Cov_Z  = zeros([4,4,nt]);
for i = 1:nt
    Cov_Z(:,:,i) = Vec10toMat4(z_new(5:end,i));
end
Cov_Z(:,:,end)

save(['files/EqLin_tot' num2str(Tot) ...
      '_e' num2str(ec1) '_' num2str(ec2) '_' num2str(ek1) '_' num2str(ek2) ...
      '_S' num2str(S0) '_fs' num2str(Fs)  '.mat'] ...
      ,'tt','mean_Z','Cov_Z','Ke','Ce');

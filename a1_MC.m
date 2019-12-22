clc
clear all
close all
tic

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

g_c = [ec1*c1 0; 0 ec2*c2];
g_k = [ek1*k1 0; 0 ek2*k2];

G1 = [zeros(2,2), eye(2); -inv(M)*K, -inv(M)*C];
G2 = [zeros(2,2), zeros(2,2); -inv(M)*g_k, -inv(M)*g_c];

% Time vector (Dimensional)
Fs   = 100;  % Sampling frequency      
Tot  = 5;  % Total time  (sec)
[tt, dt, nt] = time( Fs,Tot);   

% ----- Modulating function ----
At1 = ones([1,nt]); % sin(2.5*pi*tt/tt(end)).^2 + 0.2;
At2 = zeros([1,nt]); % sin(2.5*pi*tt/tt(end)).^2 + 0.2;

% Number of realizations
NR = 10000;
fsr1 = zeros(nt,NR);
fsr2 = zeros(nt,NR);
tspan = [tt(1) tt(end)];
x0 = [0 0 0 0];
z1 = zeros(nt,NR);
z2 = zeros(nt,NR);
z3 = zeros(nt,NR);
z4 = zeros(nt,NR);
for jj = 1:NR  
    jj
    
    fsr1(:,jj) = wgn(nt,1,db(sqrt(2*pi*S0/(dt))));
    fsr2(:,jj) = wgn(nt,1,db(sqrt(2*pi*S0/(dt))));
    
    f1 = (-1)*(At1.').*fsr1(:,jj);
    f2 = (-1)*(At2.').*fsr2(:,jj);

    [t1,response_temp] = ode45(@(t,x) ODE_MC(t,x,G1,G2,tt,f1,f2,M),tspan,x0);
    z1(:,jj) = interp1(t1,response_temp(:,1),tt,'linear');
    z2(:,jj) = interp1(t1,response_temp(:,2),tt,'linear');
    z3(:,jj) = interp1(t1,response_temp(:,3),tt,'linear');
    z4(:,jj) = interp1(t1,response_temp(:,4),tt,'linear');
end

toc

save(['files/MC45_' num2str(NR) '_tot' num2str(Tot) ...
      '_e' num2str(ec1) '_' num2str(ec2) '_' num2str(ek1) '_' num2str(ek2) ...
      '_S' num2str(S0) '_fs' num2str(Fs) '.mat'],...
      'z1','z2','z3','z4','S0','ec1','ec2','ek1','ek2','Fs','Tot');
  
plot_mean_var( z1,z2,z3,z4,dt )
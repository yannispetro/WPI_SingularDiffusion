clc 
clear all
close all

ec1 = 0.0;
ec2 = 0.0;
ek1 = 0.5;
ek2 = 0.5;

S0 = 0.1;

NR = 1000;
Tot = 20;
Fs = 100.;

tf = 1.;
Fs_WPI = 100.;
points = 11;
h = 7;

MC_file = ['files/MC45_' num2str(NR) '_tot' num2str(Tot) ...
          '_e' num2str(ec1) '_' num2str(ec2) '_' num2str(ek1) '_' num2str(ek2) ...
          '_S' num2str(S0) '_fs' num2str(Fs) '.mat'];

WPI_Ritz_file = ['files/WPI_t' num2str(tf) '_h' num2str(h) ...
            '_fs' num2str(Fs_WPI) '_p' num2str(points(1)) '.mat'];
        
WPI_EL_file = ['files/WPI_EL_tf' num2str(tf) '_S' num2str(S0) ...
               '_fs' num2str(Fs_WPI) '_p' num2str(points) '.mat'];

LinT_file = ['files/EqLin_tot' num2str(Tot) ...
             '_e' num2str(ec1) '_' num2str(ec2) '_' num2str(ek1) '_' num2str(ek2) ...
             '_S' num2str(S0) '_fs' num2str(Fs) '.mat'];

         
 % --------------------- MC ---------------------------------------
if exist(MC_file, 'file') == 2
    plot_MC( MC_file,LinT_file, tf )
    
end

if exist(WPI_Ritz_file, 'file') == 2 
    load(WPI_Ritz_file)
    x1 = linspace(domain(1), domain(2), points(1));
    x2 = linspace(domain(3), domain(4), points(2));
    x3 = linspace(domain(5), domain(6), points(3));
    x4 = linspace(domain(7), domain(8), points(4));
    
    PDF = reshape(exp(-exponent), points(1),points(2),points(3),points(4));

    plot_iter( x1, x2, x3, x4, PDF );
        
end
clc 
clear all
close all

ec1 = 0.0;
ec2 = 0.0;
ek1 = 0.5;
ek2 = 0.5;

S0 = 0.1;

NR = 10000;
Tot = 5;
Fs = 100.;

tf = 1.;
Fs_WPI = 100.;
points = 15;
h = 7;

global FIG
FIG = figure('pos',[100 300 1200 400]);

tfs = [1, 3];
for it = 1:length(tfs)
    
    tf = tfs(it);
    
    MC_file = ['files/MC45_' num2str(NR) '_tot' num2str(Tot) ...
              '_e' num2str(ec1) '_' num2str(ec2) '_' num2str(ek1) '_' num2str(ek2) ...
              '_S' num2str(S0) '_fs' num2str(Fs) '.mat'];

    WPI_Ritz_file = ['files/WPI_t' num2str(tf) '_h' num2str(h) '_S' num2str(S0) ...
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
        x1 = linspace(bounds(1), bounds(2), points(1));
        x2 = linspace(bounds(3), bounds(4), points(2));
        x3 = linspace(bounds(5), bounds(6), points(3));
        x4 = linspace(bounds(7), bounds(8), points(4));
        fprintf('==================================================================================== \n')
        fprintf('iter  mu        mean f(x)    mean c(x)     max c(x)    Nef0    Nef1     Nef2    Nef5 \n')
        fprintf('==================================================================================== \n')
        for ik = 1:7 % length(Mu)
            PDF = reshape(exp(-exponent(ik,:)), points(1),points(2),points(3),points(4));

            plot_iter_2( x1, x2, x3, x4, PDF,it,ik );

            fprintf(' %d    %4.4f      %4.2f      %5.4f        %5.4f   %4.0f       %4.0f    %4.0f    %4.0f \n',...
            ik, Mu(ik),mean(exponent(ik,:)),  mean(constraint_val(ik,:)), max(constraint_val(ik,:)), ...
            length(find(exit_flags(ik,:)==0)),length(find(exit_flags(ik,:)==1)), ...
            length(find(exit_flags(ik,:)==2)),length(find(exit_flags(ik,:)==5)))

            pause()
        end
    end
end
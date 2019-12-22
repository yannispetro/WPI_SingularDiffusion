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

tf = 3.;
Fs_WPI = 100.;
points = 15;
h = 7;
             
type = 'x2x4';
ik = 6;

colim = 'yes';
% colax = [0,17];   % x1x2 @ 1 sec
% colax = [0,0.45];   % x1x2 @ 3 sec
% colax = [0,1.05]; % x1x3 @ 1 sec
% colax = [0,0.40];   % x1x3 @ 3 sec
% colax = [0,70];   % x2x4 @ 1 sec
colax = [0,0.6];   % x2x4 @ 3 sec
% colax = [0,1.85]; % x3x4 @ 1 sec
% colax = [0,0.55]; % x3x4 @ 3 sec

MC_file = ['files/MC45_' num2str(NR) '_tot' num2str(Tot) ...
          '_e' num2str(ec1) '_' num2str(ec2) '_' num2str(ek1) '_' num2str(ek2) ...
          '_S' num2str(S0) '_fs' num2str(Fs) '.mat'];

WPI_Ritz_file = ['files/WPI_t' num2str(tf) '_h' num2str(h) '_S' num2str(S0) ...
            '_fs' num2str(Fs_WPI) '_p' num2str(points(1)) '.mat']; 
load(WPI_Ritz_file);
boundz = bounds;
      
if exist(MC_file, 'file') == 2
    nkde2 = 2^9;
%     boundz = [-5.0 5.0 -1.5 1.5 -2.7 2.7];% 15.0 sec

    load(MC_file)
    n_tf = round(tf*Fs);
    [tt, dt, N] = time( Fs,Tot );
    x1 = z1(n_tf,:);
    x2 = z2(n_tf,:);
    x3 = z3(n_tf,:);
    x4 = z4(n_tf,:);

    if strcmp(type,'x1x2')
        [~,density,Xmesh,Ymesh]=kde2d([x1.',x2.'],nkde2,[bounds(1), boundz(3)], [boundz(2), boundz(4)]);
    elseif strcmp(type,'x1x3')
        [~,density,Xmesh,Ymesh]=kde2d([x1.',x3.'],nkde2,[bounds(1), boundz(5)], [boundz(2), boundz(6)]);
    elseif strcmp(type,'x1x4')
        [~,density,Xmesh,Ymesh]=kde2d([x1.',x4.'],nkde2,[bounds(1), boundz(7)], [boundz(2), boundz(8)]);
    elseif strcmp(type,'x2x3')
        [~,density,Xmesh,Ymesh]=kde2d([x2.',x3.'],nkde2,[bounds(3), boundz(5)], [boundz(4), boundz(6)]);
    elseif strcmp(type,'x2x4')
        [~,density,Xmesh,Ymesh]=kde2d([x2.',x4.'],nkde2,[bounds(3), boundz(7)], [boundz(4), boundz(8)]);
    elseif strcmp(type,'x3x4')
        [~,density,Xmesh,Ymesh]=kde2d([x3.',x4.'],nkde2,[bounds(5), boundz(7)], [boundz(6), boundz(8)]);
    end

end

figure();
surf(Xmesh,Ymesh,density,'FaceColor','interp','LineStyle','none')
if strcmp(type,'x1x2')
    xlabel('$$x_1$$','Interpreter','latex','FontSize',40)
    ylabel('$$x_2$$','Interpreter','latex','FontSize',40)
    name = 'x1x2';
elseif strcmp(type,'x1x3')
    xlabel('$$x_1$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_1}$$','Interpreter','latex','FontSize',40)
    name = 'x1x3';
elseif strcmp(type,'x1x4')
    xlabel('$$x_1$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_2}$$','Interpreter','latex','FontSize',40)
    name = 'x1x4';
elseif strcmp(type,'x2x3')
    xlabel('$$x_2$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_1}$$','Interpreter','latex','FontSize',40)
    name = 'x2x3';
elseif strcmp(type,'x2x4')
    xlabel('$$x_2$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_2}$$','Interpreter','latex','FontSize',40)
    name = 'x2x4';
elseif strcmp(type,'x3x4')
    xlabel('$$\dot{x_1}$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_2}$$','Interpreter','latex','FontSize',40)
    name = 'x3x4';
end
zlabel('$PDF$','Interpreter','latex')
xlim([Xmesh(1),Xmesh(end)]);
ylim([Ymesh(1),Ymesh(end)]);
set(gcf,'renderer','zbuffer') 
set(gca,'fontsize',30)
% pbaspect([1 1 1])
view(2)
colorbar()
if strcmp(colim,'yes')
    caxis(colax)
end
saveas(gcf,[type,'_MCS_',num2str(tf),'sec'],'epsc')


% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3) - 4*ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

% stop


WPI_Ritz_file = ['files/WPI_t' num2str(tf) '_h' num2str(h) '_S' num2str(S0) ...
            '_fs' num2str(Fs_WPI) '_p' num2str(points(1)) '.mat'];


load(WPI_Ritz_file);
XX1 = linspace(bounds(1), bounds(2), points(1));
XX2 = linspace(bounds(3), bounds(4), points(2));
XX3 = linspace(bounds(5), bounds(6), points(3));
XX4 = linspace(bounds(7), bounds(8), points(4));
PDF = reshape(exp(-exponent(ik,:)), points(1),points(2),points(3),points(4));
s = size(PDF);

[x1,x2,f] = joint2D_PDF(XX1,XX2,XX3,XX4,PDF,s,type);

x1q = linspace(x1(1),x1(end),1001);
x2q = linspace(x2(1),x2(end),1001);
fq = interp2(x1,x2,f,x1q,x2q.','spline');

% zz = trapz(x2,f,2);
% plot(x1,zz)

figure();
surf(x1q,x2q,fq.','FaceColor','interp','LineStyle','none')
if strcmp(type,'x1x2')
    xlabel('$$x_1$$','Interpreter','latex','FontSize',40)
    ylabel('$$x_2$$','Interpreter','latex','FontSize',40)
    name = 'x1x2';
elseif strcmp(type,'x1x3')
    xlabel('$$x_1$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_1}$$','Interpreter','latex','FontSize',40)
    name = 'x1x3';
elseif strcmp(type,'x1x4')
    xlabel('$$x_1$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_2}$$','Interpreter','latex','FontSize',40)
    name = 'x1x4';
elseif strcmp(type,'x2x3')
    xlabel('$$x_2$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_1}$$','Interpreter','latex','FontSize',40)
    name = 'x2x3';
elseif strcmp(type,'x2x4')
    xlabel('$$x_2$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_2}$$','Interpreter','latex','FontSize',40)
    name = 'x2x4';
elseif strcmp(type,'x3x4')
    xlabel('$$\dot{x_1}$$','Interpreter','latex','FontSize',40)
    ylabel('$$\dot{x_2}$$','Interpreter','latex','FontSize',40)
    name = 'x3x4';
end
zlabel('$PDF$','Interpreter','latex')
xlim([x1(1),x1(end)]);
ylim([x2(1),x2(end)]);
set(gcf,'renderer','zbuffer') 
set(gca,'fontsize',30)
% pbaspect([1 1 1])
view(2)
colorbar()
if strcmp(colim,'yes')
    caxis(colax)
end
saveas(gcf,[type,'_WPI_',num2str(tf),'sec'],'epsc')
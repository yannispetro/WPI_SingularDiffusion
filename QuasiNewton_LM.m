function [X,lam,exitflag] = QuasiNewton_LM(objfun,confun,C0,Tols,bb,eta,rr)
n = length(C0);

h = 1e-6;
x1 = C0; % C0 + h*ones(n,1);

Df1 = Grad_FD(objfun,x1,h);
Dc1 = Grad_FD(confun,x1,h);

B1 = Hess_FD(objfun,x1,sqrt(h));
% B1 = diag(diag(Hess_FD(objfun,x1,sqrt(h))));
H1 = inv(B1);
% H1 = diag(rand(n,1));

% CON = zeros(1,1000);
% LLL = zeros(1,1000);
% ALP = zeros(1,1000);

Max_iter = 1000;
k = 1;
alpha = 1;
dx = 100*ones(n,1);
mu = -100000;
dFtol = 10000;
% adx = [];
% while ( dFtol > Tols(1) && alpha*norm(dx) > Tols(2) ) && k < Max_iter  
while ( dFtol > Tols(1) && norm(dx) > Tols(2) ) && k < Max_iter   
%     adx = [adx alpha*norm(dx)];

    CC1 = confun(x1);
    
    if (Dc1*H1*Dc1.') == 0
        lam = 0;
    else
        lam = (-CC1 + Dc1*H1*Df1.')/(Dc1*H1*Dc1.');
    end
    
    dx = H1*( lam*Dc1.' - Df1.' );
    
%     if sum(eig(B1)) > 0
%         sigma = 1;
%     else
%         sigma = 0;
%     end
    sigma = 1;
    mu_new = ( Df1*dx + 0.5*sigma*dx.'*B1*dx )/( (1-rr)*abs(CC1) );
    if mu < mu_new
        mu = mu_new + 0.001;
    end
    
    alpha = 1;
    phi11 = objfun(x1 + alpha*dx) + mu*abs(confun(x1 + alpha*dx));
    phi12 = objfun(x1) + mu*abs(CC1);
    Dphi1 = Df1*dx - mu*abs(CC1);
    while phi11 > phi12 + eta*alpha*Dphi1 % && alpha > 1e-4
        alpha = bb*alpha;

        phi11 = objfun(x1 + alpha*dx) + mu*abs(confun(x1 + alpha*dx));
    end
    
    x2 = x1 + alpha*dx;
    
    dFtol = abs(objfun(x2) - objfun(x1))/abs(objfun(x1));
    
    Df2 = Grad_FD(objfun,x2,h); 
    Dc2 = Grad_FD(confun,x2,h); 
    
    s = x2 - x1;
    y = Df2.' - Df1.' - lam*( Dc2.' - Dc1.' );
    
%     disp(s.'*y)
    
    rho = (y.'*s)^(-1);
    H2 = (eye(n) - rho*s*y.')*H1*(eye(n) - rho*y*s.') + rho*(s*s.');
%     B2 = B1 - ( B1*(s*s.')*B1 )/(s.'*B1*s) + y*y.'/(y.'*s);
    
%     CC2 = confun(x2);
%     LL2 = objfun(x2) - lam*CC2;
    
    x1 = x2;
    Df1 = Df2;
    Dc1 = Dc2;
    
%     disp(Df1 - lam*Dc1)

%     H1 = H2;
%     B1 = B2;
    H1 = H2;
    B1 = inv(H2);
    
%     LLL(k) = LL2;
%     CON(k) = CC2;
%     ALP(k) = objfun(x2);
%     fprintf('%d  fun =  %d, con =   %d, lambda = %d, L = %d, alpha = %d \n',k,objfun(x2),CC2,lam,LL2,alpha)
        
    k = k + 1;
end
% figure(1)
% plot(adx)
% 
% figure(2)
% plot(CON(1:k))
% 
% figure(3)
% plot(ALP(1:k))

if k == 1000
    exitflag = 0;
elseif dFtol <= Tols(1)
    exitflag = 1;
elseif alpha*norm(dx) <= Tols(2)
    exitflag = 2;
end

X = x2;


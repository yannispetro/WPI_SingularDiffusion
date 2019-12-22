function [exponent,constraint_val,exit_flags,C_MPP,LM,elt] = a4_functional_minimization_Ritz(...
          domain, points, ndof, ord, h, ...
          par, ti, tf, Fs)

%     m1 = par(1);
%     m2 = par(2);
%     c1 = par(3);
%     c2 = par(4);
%     k1 = par(5);
%     k2 = par(6);
% 
%     ec1 = par(7);
%     ec2 = par(8);
%     ek1 = par(9);
%     ek2 = par(10);
    
    S0 = par(11);
      
    % Create the original grid 
    x1 = linspace(domain(1), domain(2), points(1));
    x2 = linspace(domain(3), domain(4), points(2));
    x3 = linspace(domain(5), domain(6), points(3));
    x4 = linspace(domain(7), domain(8), points(4));

    [X1, X2, X3, X4] = ndgrid(x1, x2, x3, x4);
    %
    num = points(1)*points(2)*points(3)*points(4);
    X1 = reshape(X1, num, 1);
    X2 = reshape(X2, num, 1);
    X3 = reshape(X3, num, 1);
    X4 = reshape(X4, num, 1);

    grid = [X1 X2 X3 X4];

    nmax_needed = num;
      
    % - Time discretization
    dt = 1/Fs;                          % Sampling period (Original - sec)
    nt  = floor((tf - ti)/dt) + 1;      % Length of signal
    tt = (ti + (0:nt-1)*dt).';          % Time vector
    
    % - Construction shifted Legendre polynomials -
    syms t
    P = sym('P', [1 h]);
    P(1) = 1;
    P(2) = (2 * t - ti - tf) / (tf - ti);

    for kk = 3:h
        p = kk - 2;
        P(kk) = expand((2 * p + 1) / (p + 1) * (2 * t - ti - tf)...
            / (tf - ti) * P(kk - 1) - p / (p + 1) * P(kk - 2));
    end
    clear kk p

    % - g is the function that each of its elements multiplies one c_i 
    g = ((t - ti) ^ ord) * ((t - tf) ^ ord) * P;
    dgdt = diff(g, t);
    d2gdt2 = diff(g, t, 2);

    % - Construction Hermite interpolating polynomial -
    Ht = sym('Ht', [1 2 * ord]);
    Ht(1) = 1;
    for kk = 2:2*ord
        Ht(kk) = Ht(kk - 1) * t;
    end
    dnHdtn = sym('dnHdtn', [ord 2*ord]);
    for kk = 1:ord
        dnHdtn(kk, :) = diff(Ht, t, (kk - 1));
    end
    B = sym('B',[2*ord 1]);
    A = [subs(dnHdtn, t, ti)
        subs(dnHdtn, t, tf)];

    aa = A \ B;
    H = Ht * aa;
    dHdt = diff(H, t);
    d2Hdt2 = diff(dHdt, t);
    clear Ht kk dnHdtn A
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    g0 = double(subs(g, tt));
    g1 = double(subs(dgdt, tt));
    g2 = double(subs(d2gdt2, tt));

    H0 = subs(H, tt);
    H1 = subs(dHdt, tt);
    H2 = subs(d2Hdt2, tt);
    
    MTRX_H = zeros([3, nt, 2*ord]);
    MTRX_H(1,:,:) = double(equationsToMatrix(H0, B));
    MTRX_H(2,:,:) = double(equationsToMatrix(H1, B));
    MTRX_H(3,:,:) = double(equationsToMatrix(H2, B));
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Nvar = h;
%     figure(); plot(tt,g0)
%     figure(); plot(tt,g1)
%     figure(); plot(tt,g2)
%     pause()
% - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Boundary conditions
    % 4 because we have 2 initial constraints and 2 final
    BC = zeros(ndof, 2*ord, nmax_needed); 
    BC(1,3,:) = grid(:, 1);
    BC(2,3,:) = grid(:, 2);
    BC(1,4,:) = grid(:, 3);
    BC(2,4,:) = grid(:, 4);
    
    % --- MAIN LOOP ----
    Constraint_integral = @(C, Herm) trapz(tt,(Constraint_function_C( ...
                             Nvar,C,g0,g1,g2,Herm,par)).^2 );

    func = @(C, Herm) stochastic_action(Nvar, C, g0, g1, g2, Herm, ...
                                             par, tt);

    
    exponent = zeros([1, nmax_needed]);
    exit_flags = zeros([1, nmax_needed]);
    constraint_val = zeros([1, nmax_needed]);
    C_MPP = zeros([nmax_needed, ndof*Nvar]);
    C0 = zeros(ndof*Nvar,1);
%     C0 = 0.1*[ones([Nvar 1]); -ones([Nvar 1])];
    LM = zeros([1,nmax_needed]);
    
    tic;
    parfor_progress(nmax_needed);
    parfor jj = 1:nmax_needed
        warning('off')
        % Hermitian polynomial that takes into account the BCs
        Herm = zeros([nt, 3*ndof]);
        for i = 1:ndof
            for j = 1:3
                Herm(:,3*(i-1)+j) = squeeze(MTRX_H(j,:,:)) * BC(i,:,jj).';
            end
        end 

        objfun  = @(q) func(q, Herm);
        confun  = @(q) Constraint_integral(q,Herm);

        bb  = 0.3; % (0,1) 
        eta = 0.1; % (0,0.5)
        rr  = 0.5; % (0,1)
        TolsQN = [1e-15,1e-1];
        
        [X,lam,exitflag] = QuasiNewton_LM(objfun,confun,C0,TolsQN,bb,eta,rr);

        C_MPP(jj,:) = X;
        LM(jj) = lam;


        exponent(jj) = stochastic_action(Nvar, X, g0, g1, g2, Herm, ...
                                         par, tt);

        exit_flags(jj) = exitflag;
        constraint_val(jj) = confun(X);

        parfor_progress;
    end
    parfor_progress(0);
    
    elt=toc;
    fprintf('Total elapsed time = %f seconds \n', elt)
end
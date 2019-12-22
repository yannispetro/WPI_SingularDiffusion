function J = stochastic_action(h, C, g0, g1, g2, Herm, ...
                               par, tt)
    
    m1 = par(1);
%     m2 = par(2);
    c1 = par(3);
    c2 = par(4);
    k1 = par(5);
    k2 = par(6);

    ec1 = par(7);
%     ec2 = par(8);
    ek1 = par(9);
%     ek2 = par(10);
    
    S0 = par(11);
                           
    x10 = g0 * C(1:h) + Herm(:,1);
    x11 = g1 * C(1:h) + Herm(:,2);
    x12 = g2 * C(1:h) + Herm(:,3);
    
    x20 = g0 * C(h+1:2*h) + Herm(:,4);
    x21 = g1 * C(h+1:2*h) + Herm(:,5);
%     x22 = g2 * C(h+1:2*h) + Herm(:,6);

    L = (1/4)*pi^(-1)*S0^(-1)*((k1+k2)*x10+ek1*k1*x10.^3+(-1)*k2*x20+(c1+c2)*x11+c1*ec1*x11.^3+(-1)*c2*x21+m1*x12).^2;
    
    J = trapz(tt,L);
end

function x_dot = ODE_MC( t,x,G1,G2,tt,f1,f2,My )

    ft1 = interp1(tt,f1,t,'linear');
    ft2 = interp1(tt,f2,t,'linear');
    M = inv(My);
    
    x_dot = G1*x + G2*x.^3 + [ 0; 0; M(1,1)*ft1 + M(1,2)*ft2; M(2,1)*ft1 + M(2,2)*ft2];

end


function Constr = Constraint_function_X(x10,x20,x11,x21,x22,par)

%     m1 = par(1);
    m2 = par(2);
%     c1 = par(3);
    c2 = par(4);
%     k1 = par(5);
    k2 = par(6);

%     ec1 = par(7);
    ec2 = par(8);
%     ek1 = par(9);
    ek2 = par(10);
    
%     S0 = par(11);

    Constr = (-1)*k2*x10+k2*x20+ek2*k2*x20.^3+(-1)*c2*x11+c2*x21+c2*ec2*x21.^3+m2*x22;

end



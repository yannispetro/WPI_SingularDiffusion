function Constr = Constraint_function_C(h, C, g0, g1, g2, Herm, par)
                           
    x10 = g0 * C(1:h) + Herm(:,1);
    x11 = g1 * C(1:h) + Herm(:,2);
%     x12 = g2 * C(1:h) + Herm(:,3);
    
    x20 = g0 * C(h+1:2*h) + Herm(:,4);
    x21 = g1 * C(h+1:2*h) + Herm(:,5);
    x22 = g2 * C(h+1:2*h) + Herm(:,6);

    Constr = Constraint_function_X(x10,x20,x11,x21,x22,par);

end


function Df = Grad_FD(fun,x0,h)
    n = length(x0);
    Df = zeros(1,n);
    for i = 1:n
        xi = zeros(n,1);
        xi(i) = 1;
        x1 = x0 - h*xi;
        x2 = x0 + h*xi;
        Df(i) = ( fun(x2) - fun(x1) )/( 2*h );
    end
end

function Hf = Hess_FD(f,x0,h)
    n = length(x0);
    Hf = zeros(n);
    for i = 1:n
        for j = i:n
            xi = zeros(n,1);
            xj = zeros(n,1);
            xi(i) = 1;
            xj(j) = 1;
            if j == i
                Hf(i,i) = ( f(x0 + h*xi) - f(x0) + f(x0 - h*xi) )/( h^2 );
            else
                Hf(i,j) = ( f(x0 + h*xi + h*xj) - f(x0 - h*xi + h*xj) - f(x0 + h*xi - h*xj) + f(x0 - h*xi - h*xj) )/( 4*h^2 );
                Hf(j,i) = Hf(i,j); 
            end
        end
    end
end

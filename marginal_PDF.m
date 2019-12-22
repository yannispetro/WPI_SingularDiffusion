function [ f1,f2,f3,f4 ] = marginal_PDF( x1,x2,x3,x4,PDF )
s = size(PDF);

f1 = zeros(1, s(1));
f2 = zeros(1, s(2));
f3 = zeros(1, s(3));
f4 = zeros(1, s(4));
dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);
dx3 = x3(2)-x3(1);
dx4 = x4(2)-x4(1);
for i=1:s(1)
    for j=1:s(2)
        for k=1:s(3)
            for l=1:s(4)
                f1(i) = f1(i) + PDF(i,j,k,l)*dx2*dx3*dx4;
                f2(j) = f2(j) + PDF(i,j,k,l)*dx1*dx3*dx4;
                f3(k) = f3(k) + PDF(i,j,k,l)*dx1*dx2*dx4;
                f4(l) = f4(l) + PDF(i,j,k,l)*dx1*dx2*dx3;
            end
        end
    end
end
total1 = trapz(x1, f1);
total2 = trapz(x2, f2);
total3 = trapz(x3, f3);
total4 = trapz(x4, f4);
f1 = f1/total1;
f2 = f2/total2;
f3 = f3/total3;
f4 = f4/total4;
end


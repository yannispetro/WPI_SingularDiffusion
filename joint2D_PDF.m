function [xmesh1,xmesh2,f] = joint2D_PDF(x1,x2,x3,x4,PDF,s,type)
dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);
dx3 = x3(2)-x3(1);
dx4 = x4(2)-x4(1);
if strcmp(type,'x1x2')
    f = zeros([s(1),s(2)]);
    for i=1:s(1)
        for j=1:s(2)
            for k=1:s(3)
                for l=1:s(4)
                    f(i,j) = f(i,j) + PDF(i,j,k,l)*dx3*dx4;
                end
            end
        end
    end
    total = trapz(x2,trapz(x1,f,1),2);
    f = f/total;
    xmesh1 = x1;
    xmesh2 = x2;
elseif strcmp(type,'x1x3')
    f = zeros([s(1),s(3)]);
    for i=1:s(1)
        for j=1:s(2)
            for k=1:s(3)
                for l=1:s(4)
                    f(i,k) = f(i,k) + PDF(i,j,k,l)*dx2*dx4;
                end
            end
        end
    end
    total = trapz(x3,trapz(x1,f,1),2);
    f = f/total;
    xmesh1 = x1;
    xmesh2 = x3;
elseif strcmp(type,'x1x4')
    f = zeros([s(1),s(4)]);
    for i=1:s(1)
        for j=1:s(2)
            for k=1:s(3)
                for l=1:s(4)
                    f(i,l) = f(i,l) + PDF(i,j,k,l)*dx2*dx3;
                end
            end
        end
    end
    total = trapz(x4,trapz(x1,f,1),2);
    f = f/total;
    xmesh1 = x1;
    xmesh2 = x4;
elseif strcmp(type,'x2x3')
    f = zeros([s(2),s(3)]);
    for i=1:s(1)
        for j=1:s(2)
            for k=1:s(3)
                for l=1:s(4)
                    f(j,k) = f(j,k) + PDF(i,j,k,l)*dx1*dx4;
                end
            end
        end
    end
    total = trapz(x3,trapz(x2,f,1),2);
    f = f/total;
    xmesh1 = x2;
    xmesh2 = x3;
elseif strcmp(type,'x2x4')
    f = zeros([s(2),s(4)]);
    for i=1:s(1)
        for j=1:s(2)
            for k=1:s(3)
                for l=1:s(4)
                    f(j,l) = f(j,l) + PDF(i,j,k,l)*dx1*dx3;
                end
            end
        end
    end
    total = trapz(x4,trapz(x2,f,1),2);
    f = f/total;
    xmesh1 = x2;
    xmesh2 = x4;
elseif strcmp(type,'x3x4')
    f = zeros([s(3),s(4)]);
    for i=1:s(1)
        for j=1:s(2)
            for k=1:s(3)
                for l=1:s(4)
                    f(k,l) = f(k,l) + PDF(i,j,k,l)*dx1*dx2;
                end
            end
        end
    end
    total = trapz(x4,trapz(x3,f,1),2);
    f = f/total;
    xmesh1 = x3;
    xmesh2 = x4;
end


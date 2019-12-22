function [Xmesh, Ymesh, density] = kde2d_modified(t, dataX, dataY, Rep, dt, limitX, limitY, points)

    [~, density, Xmesh, Ymesh] = kde2d([dataX(1 : Rep, round(t / dt + 1)) dataY(1 : Rep, round(t / dt + 1))], points, [-limitX, -limitY], [limitX, limitY]);

end
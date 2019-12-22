function [] = plot_iter_forPaper( x1, x2, x3, x4, PDF, mrk )
    global H1 H2 H3 H4 FIG1 FIG2 FIG3 FIG4
    
    [f1,f2,f3,f4] = marginal_PDF(x1,x2,x3,x4,PDF);
%     figure(FIG);
    
    if exist('H1', 'var')
%         delete(H1)
        figure(1)
        H1 = plot(x1,f1, mrk,'LineStyle',':','Linewidth',2,'MarkerSize',10); hold on
    end
    if exist('H2', 'var')
%         delete(H2)
       figure(2)
        H2 = plot(x2,f2, mrk,'LineStyle',':','Linewidth',2,'MarkerSize',10); hold on
    end
    if exist('H3', 'var')
%         delete(H3)
        figure(3)
        H3 = plot(x3,f3, mrk,'LineStyle',':','Linewidth',2,'MarkerSize',10); hold on
    end
    if exist('H4', 'var')
%         delete(H4)
        figure(4)
        H4 = plot(x4,f4, mrk,'LineStyle',':','Linewidth',2,'MarkerSize',10); hold on
    end
    drawnow
end


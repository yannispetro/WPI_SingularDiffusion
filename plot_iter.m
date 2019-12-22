function [] = plot_iter( x1, x2, x3, x4, PDF )
    global H1 H2 H3 H4 FIG
    
    [f1,f2,f3,f4] = marginal_PDF(x1,x2,x3,x4,PDF);
%     figure(FIG);
    
    if exist('H1', 'var')
        delete(H1)
%         figure(1);
        subplot(1,4,1);
        H1 = plot(x1,f1, 'x','LineStyle',':','Color',0.6*[1 1 1],'Linewidth',2,'MarkerSize',10); hold on
        xlim([x1(1),x1(end)])
        xlabel('$$x$$','Interpreter','latex','FontSize',20)
    end
    if exist('H2', 'var')
        delete(H2)
%         figure(2);
        subplot(1,4,2);
        H2 = plot(x2,f2, 'x','LineStyle',':','Color',0.6*[1 1 1],'Linewidth',2,'MarkerSize',10); hold on
        xlim([x2(1),x2(end)])
        xlabel('$$y$$','Interpreter','latex','FontSize',20)
    end
    if exist('H3', 'var')
        delete(H3)
%         figure(3);
        subplot(1,4,3);
        H3 = plot(x3,f3, 'x','LineStyle',':','Color',0.6*[1 1 1],'Linewidth',2,'MarkerSize',10); hold on
        xlim([x3(1),x3(end)])
        xlabel('$$\dot{x}$$','Interpreter','latex','FontSize',20)
    end
    if exist('H4', 'var')
        delete(H4)
        subplot(1,4,4);
        H4 = plot(x4,f4, 'x','LineStyle',':','Color',0.6*[1 1 1],'Linewidth',2,'MarkerSize',10); hold on
        xlim([x4(1),x4(end)])
        xlabel('$$\dot{y}$$','Interpreter','latex','FontSize',20)
    end
    drawnow
end


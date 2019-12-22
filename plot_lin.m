function [] = plot_lin( x1, x2, x3, x4, f1, f2 ,f3, f4 )
    global FIG
    
%     if exist('FIG', 'var')
%         disp('')
%     else
%         figure(FIG);
%     end
    
    
    
    subplot(1,4,1);
    plot(x1,f1, 'x','LineStyle','-','Color',0.2*[1 1 1],'Linewidth',1,'MarkerSize',1); hold on

    subplot(1,4,2);
    plot(x2,f2, 'x','LineStyle','-','Color',0.2*[1 1 1],'Linewidth',1,'MarkerSize',1); hold on

    subplot(1,4,3);
    plot(x3,f3, 'x','LineStyle','-','Color',0.2*[1 1 1],'Linewidth',1,'MarkerSize',1); hold on

    subplot(1,4,4);
    plot(x4,f4, 'x','LineStyle','-','Color',0.2*[1 1 1],'Linewidth',1,'MarkerSize',1); hold on

    drawnow
end


function [] = plot_MC( MC_file,LinT_file, tf )

    if exist(MC_file, 'file') == 2
        load(MC_file)
        n_tf = round(tf*Fs);
        x     = z1(1:n_tf,:);
        y     = z2(1:n_tf,:);
        x_dot = z3(1:n_tf,:);
        y_dot = z4(1:n_tf,:);

        [tt, dt, N] = time( Fs,Tot );

        [bandwidth,pdf_x_MC,xmesh1,cdf] = kde(x(end,:));
        [bandwidth,pdf_y_MC,xmesh2,cdf] = kde(y(end,:));
        [bandwidth,pdf_x_dot_MC,xmesh3,cdf] = kde(x_dot(end,:));
        [bandwidth,pdf_y_dot_MC,xmesh4,cdf] = kde(y_dot(end,:));
        
        global FIG
%         if exist('FIG', 'var')
%             disp('')
%         else
%             FIG = figure('pos',[100 300 1200 400]);
%         end
%         FIG = figure('pos',[100 300 600 600]);
        set(gcf, 'Position', [100 300 1200 400])
        
        subplot(1,4,1); plot(xmesh1,pdf_x_MC    ,'LineStyle','-','Color',0.8*[1 1 1],'Linewidth',5); hold on
        subplot(1,4,2); plot(xmesh2,pdf_y_MC    ,'LineStyle','-','Color',0.8*[1 1 1],'Linewidth',5); hold on
        subplot(1,4,3); plot(xmesh3,pdf_x_dot_MC,'LineStyle','-','Color',0.8*[1 1 1],'Linewidth',5); hold on
        subplot(1,4,4); plot(xmesh4,pdf_y_dot_MC,'LineStyle','-','Color',0.8*[1 1 1],'Linewidth',5); hold on

%         subplot(2,2,1); plot(xmesh1,pdf_x_MC,'LineStyle','-','Linewidth',3); hold on
%         subplot(2,2,2); plot(xmesh2,pdf_y_MC,'LineStyle','-','Linewidth',3); hold on
%         subplot(2,2,3); plot(xmesh3,pdf_x_dot_MC,'LineStyle','-','Linewidth',3); hold on
%         subplot(2,2,4); plot(xmesh4,pdf_y_dot_MC,'LineStyle','-','Linewidth',3); hold on

        clear x y x_dot y_dot
        
%         if exist(LinT_file, 'file') == 2
%             load(LinT_file)
% 
%             m_x1 = mean_Z(1,n_tf);
%             m_x2 = mean_Z(2,n_tf);
%             m_x3 = mean_Z(3,n_tf);
%             m_x4 = mean_Z(4,n_tf);
%             var_x1 = Cov_Z(1,1,n_tf);
%             var_x2 = Cov_Z(2,2,n_tf);
%             var_x3 = Cov_Z(3,3,n_tf);
%             var_x4 = Cov_Z(4,4,n_tf);
% 
%             pdf_x1_Lin     = 1/(sqrt(var_x1*2*pi))    *exp(-(xmesh1-m_x1).^2/(2*var_x1));
%             pdf_x2_Lin     = 1/(sqrt(var_x2*2*pi))    *exp(-(xmesh2-m_x2).^2/(2*var_x2));
%             pdf_x3_Lin     = 1/(sqrt(var_x3*2*pi))    *exp(-(xmesh3-m_x3).^2/(2*var_x3));
%             pdf_x4_Lin     = 1/(sqrt(var_x4*2*pi))    *exp(-(xmesh4-m_x4).^2/(2*var_x4));
% 
%             plot_lin(xmesh1,xmesh2,xmesh3,xmesh4,pdf_x1_Lin,pdf_x2_Lin,pdf_x3_Lin,pdf_x4_Lin)
% 
%         end
        
        
    end

end


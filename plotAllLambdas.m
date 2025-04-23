function plotAllLambdas()
    clc; clear; close all;
    T = readtable('results_table.csv');
    lambdas    = unique(T.Lambda);
    leftPhotos = [1,3,4];
    kernels    = [1,3,5];
    
    cmap       = lines(2);
    grayColor  = [0.3 0.3 0.3];
    barWidth   = 0.4;
    colorShift = 0.5;

    for li = 1:numel(lambdas)
        lam = lambdas(li);
        sub = T(T.Lambda==lam, :);

        % build matrices
        IterG = nan(3,3); IterC = nan(3,3);
        Egray = nan(3,3); Ecol  = nan(3,3);
        for i=1:3
          for j=1:3
            r = sub(sub.LeftPhoto==leftPhotos(j) & sub.Kernel==kernels(i), :);
            IterG(i,j) = r.Iter_gray;
            IterC(i,j) = r.Iter_color;
            Egray(i,j) = r.E_gray_end;
            Ecol(i,j)  = r.E_color_end;
          end
        end

        % wide figure
        fig = figure('Units','normalized','Position',[0.1 0.1 0.8 0.6],'Name',sprintf('\\lambda = %.2f',lam));

        % Left subplot: Iterations
        ax1 = subplot('Position',[0.05 0.1 0.425 0.8]);
        hold(ax1,'on')
        hG = bar3(ax1,IterG,barWidth);
        for k=1:numel(hG), set(hG(k),'FaceColor',cmap(1,:),'EdgeColor','k'); end
        hC = bar3(ax1,IterC,barWidth);
        for k=1:numel(hC)
          yd = get(hC(k),'YData'); set(hC(k),'YData',yd+colorShift);
          set(hC(k),'FaceColor',cmap(2,:),'EdgeColor','k');
        end
        for i=1:3, for j=1:3
          text(ax1,j,i,IterG(i,j),sprintf('%d',IterG(i,j)), ...
               'HorizontalAlignment','center','VerticalAlignment','bottom','Color',grayColor);
          text(ax1,j,i+colorShift,IterC(i,j),sprintf('%d',IterC(i,j)), ...
               'HorizontalAlignment','center','VerticalAlignment','bottom','Color',grayColor);
        end,end
        view(ax1,45,25)
        xlabel(ax1,'Distance between the scenes'); ylabel(ax1,'Kernel'); zlabel(ax1,'Iterations')
        title(ax1,sprintf('\\lambda=%.2f — Iterations',lam))
        lg1=legend(ax1,[hG(1),hC(1)],{'Gray','Color'},'Location','northeastoutside');
        lg1.TextColor = grayColor;
        xticks(ax1,1:3); xticklabels(ax1,string(5-leftPhotos));
        yticks(ax1,1:3); yticklabels(ax1,string(kernels));
        axis(ax1,'square')
        hold(ax1,'off')

        % Right subplot: Energy
        ax2 = subplot('Position',[0.525 0.1 0.425 0.8]);
        hold(ax2,'on')
        hG2 = bar3(ax2,Egray,barWidth);
        for k=1:numel(hG2), set(hG2(k),'FaceColor',cmap(1,:),'EdgeColor','k'); end
        hC2 = bar3(ax2,Ecol,barWidth);
        for k=1:numel(hC2)
          yd = get(hC2(k),'YData'); set(hC2(k),'YData',yd+colorShift);
          set(hC2(k),'FaceColor',cmap(2,:),'EdgeColor','k');
        end
        for i=1:3, for j=1:3
          text(ax2,j,i,Egray(i,j),sprintf('%.1f',Egray(i,j)), ...
               'HorizontalAlignment','center','VerticalAlignment','bottom','Color',grayColor);
          text(ax2,j,i+colorShift,Ecol(i,j),sprintf('%.1f',Ecol(i,j)), ...
               'HorizontalAlignment','center','VerticalAlignment','bottom','Color',grayColor);
        end,end
        view(ax2,45,25)
        xlabel(ax2,'Distance between the scenes'); ylabel(ax2,'Kernel'); zlabel(ax2,'Energy')
        title(ax2,sprintf('\\lambda=%.2f — Final Energy',lam))
        lg2=legend(ax2,[hG2(1),hC2(1)],{'Gray','Color'},'Location','northeastoutside');
        lg2.TextColor = grayColor;
        xticks(ax2,1:3); xticklabels(ax2,string(5-leftPhotos));
        yticks(ax2,1:3); yticklabels(ax2,string(kernels));
        axis(ax2,'square')
        hold(ax2,'off')

        % save
        saveas(fig, sprintf('lambda_%.2f_summary.png',lam))
        close(fig)
    end
end

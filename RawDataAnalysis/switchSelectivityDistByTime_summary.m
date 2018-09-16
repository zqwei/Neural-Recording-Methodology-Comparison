frac_ = [0.1319    0.6028    0.2653; 
         0.2934    0.6557    0.0509; 
         0.4592    0.4980    0.0427;
         0.4292    0.5194    0.0514;
         0.4625    0.4917    0.0458];
     
sum_size = [720; 1493; 2293; 720; 720];
ylim_    = [0.7, 0.7, 0.3];

figure
for nPlot = 1:3
    subplot(1, 3, nPlot)
    hold on
    bar(1:5, frac_(:, nPlot), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    errorbar(1:5, frac_(:, nPlot), ...
        sqrt(frac_(:, nPlot).*(1-frac_(:, nPlot))./sum_size), 'k')
    set(gca, 'TickDir', 'out')
    xlim([0.5, 5.5])
    ylim([0 ylim_(nPlot)])
    set(gca, 'XTick', 1:5)
    set(gca, 'YTick', 0:0.1:ylim_(nPlot))
end
setPrint(4*3, 3, 'cell_type', 'pdf')

frac_ = [0.1689    0.6267    0.2044; 
         0.5210    0.4528    0.0262; 
         0.5867    0.3689    0.0444];
     
sum_size = [225; 2672; 225];

figure
for nPlot = 1:3
    subplot(1, 3, nPlot)
    hold on
    bar(1:3, frac_(:, nPlot), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    errorbar(1:3, frac_(:, nPlot), ...
        sqrt(frac_(:, nPlot).*(1-frac_(:, nPlot))./sum_size), 'k')
    set(gca, 'TickDir', 'out')
    xlim([0.5, 3.5])
    ylim([0 ylim_(nPlot)])
    set(gca, 'XTick', 1:3)
    set(gca, 'YTick', 0:0.1:ylim_(nPlot))
end
setPrint(4*3, 3, 'cell_type_fast', 'pdf')

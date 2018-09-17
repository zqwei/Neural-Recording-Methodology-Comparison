addpath('../Func');
setDir;
load('refMat.mat')
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params        = DataSetList(3).params;

% dataSetNames  = {DataSetList(4).name, DataSetList(3).name, DataSetList(10).name};
% dataSetNames_   = {'Modeled_6s_AAV_nRSParaSet_', 'Modeled_GP43_nRSParaSet_', 'Modeled_GP517_nRSParaSet_'};

dataSetNames  = {DataSetList(10).name, DataSetList(3).name, DataSetList(4).name};
dataSetNames_   = {'Modeled_GP517_nRSParaSet_','Modeled_GP43_nRSParaSet_', 'Modeled_6s_AAV_nRSParaSet_'};


line_color = [     0    0.4470    0.7410;
              0.3010    0.7450    0.9330;
              0.4940    0.1840    0.5560];

figure;
subplot(2, 3, 1)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_s = tmp(2,:)'/num_';
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(2,:)'/num_';
    
    xi = 0.3:0.001:0.7; 
    fs = ksdensity(tmp_s, xi);
    fc = ksdensity(tmp_c, xi);
    f_ = max([fs, fc]);
    
    fill([xi, xi(1)], [fs, fs(1)]/f_*0.97+nData, nData, 'facecolor', line_color(nData, :), 'facealpha', 0.5, 'edgecolor', 'none');
    plot(xi, fc/f_*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(quantile(tmp_s, [.25 .75]), [nData, nData], 'linewid', 2, 'color', [0.5 0.5 0.5])
    plot(median(tmp_s), nData, 'o', 'linewid', 2, 'markerfacecolor', line_color(nData, :), 'color', [0.5 0.5 0.5])
    plot(quantile(tmp_c, [.25 .75]), [nData, nData], 'linewid', 2, 'color', line_color(nData, :))
    plot(median(tmp_c), nData, 'o', 'linewid', 2, 'markerfacecolor', 'w', 'color', line_color(nData, :))
end
ylim([0.9, 4.1])
xlim([0.299 0.701])
set(gca,'Ytick', 1.5:3.5)
plot([performanceMat.mono], 3:-1:1, 'k+', 'linewid', 1)
plot([refEphys.mono, refEphys.mono, refEphysFast.mono], 3:-1:1, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 2)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_s = tmp(3,:)'/num_';
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(3,:)'/num_';
    
    xi = 0:0.001:0.32; 
    fs = ksdensity(tmp_s, xi);
    fc = ksdensity(tmp_c, xi);
    f_ = max([fs, fc]);
    
    fill([xi, xi(1)], [fs, fs(1)]/f_*0.97+nData, nData, 'facecolor', line_color(nData, :), 'facealpha', 0.5, 'edgecolor', 'none');
    plot(xi, fc/f_*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(quantile(tmp_s, [.25 .75]), [nData, nData], 'linewid', 2, 'color', [0.5 0.5 0.5])
    plot(median(tmp_s), nData, 'o', 'linewid', 2, 'markerfacecolor', line_color(nData, :), 'color', [0.5 0.5 0.5])
    plot(quantile(tmp_c, [.25 .75]), [nData, nData], 'linewid', 2, 'color', line_color(nData, :))
    plot(median(tmp_c), nData, 'o', 'linewid', 2, 'markerfacecolor', 'w', 'color', line_color(nData, :))
end
ylim([0.9, 4.1])
xlim([0 0.32])
set(gca,'Ytick', 1.5:3.5)
plot([performanceMat.multi], 3:-1:1, 'k+', 'linewid', 1)
plot([refEphys.multi, refEphys.multi, refEphysFast.multi], 3:-1:1, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 3)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp_s = [analysisMat.peakiness];
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c = [analysisMat.peakiness];
    
    xi = 0:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    fc = ksdensity(tmp_c, xi);
    f_ = max([fs, fc]);
    
    fill([xi, xi(1)], [fs, fs(1)]/f_*0.97+nData, nData, 'facecolor', line_color(nData, :), 'facealpha', 0.5, 'edgecolor', 'none');
    plot(xi, fc/f_*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(quantile(tmp_s, [.25 .75]), [nData, nData], 'linewid', 2, 'color', [0.5 0.5 0.5])
    plot(median(tmp_s), nData, 'o', 'linewid', 2, 'markerfacecolor', line_color(nData, :), 'color', [0.5 0.5 0.5])
    plot(quantile(tmp_c, [.25 .75]), [nData, nData], 'linewid', 2, 'color', line_color(nData, :))
    plot(median(tmp_c), nData, 'o', 'linewid', 2, 'markerfacecolor', 'w', 'color', line_color(nData, :))
end
ylim([0.9, 4.1])
xlim([0.29 1.51])
set(gca,'Ytick', 1.5:3.5)
plot([performanceMat.peak], 3:-1:1, 'k+', 'linewid', 1)
plot([refEphys.peak, refEphys.peak, refEphysFast.peak], 3:-1:1, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 4)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp_s       = nan(1000, 1);
    for n_      = 1:1000
        PCAVar  = analysisMat(n_).PCAVar;
        tmp_s(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        PCAVar  = analysisMat(n_).PCAVar;
        tmp_c(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end
    
    xi = -0.4:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    fc = ksdensity(tmp_c, xi);
    f_ = max([fs, fc]);
    
    fill([xi, xi(1)], [fs, fs(1)]/f_*0.97+nData, nData, 'facecolor', line_color(nData, :), 'facealpha', 0.5, 'edgecolor', 'none');
    plot(xi, fc/f_*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(quantile(tmp_s, [.25 .75]), [nData, nData], 'linewid', 2, 'color', [0.5 0.5 0.5])
    plot(median(tmp_s), nData, 'o', 'linewid', 2, 'markerfacecolor', line_color(nData, :), 'color', [0.5 0.5 0.5])
    plot(quantile(tmp_c, [.25 .75]), [nData, nData], 'linewid', 2, 'color', line_color(nData, :))
    plot(median(tmp_c), nData, 'o', 'linewid', 2, 'markerfacecolor', 'w', 'color', line_color(nData, :))
    
    refPCAImage(nData) = 1 - performanceMat(nData).pca(2)/sum(performanceMat(nData).pca);
end
ylim([0.9, 4.1])
xlim([0 1.01])
set(gca,'Ytick', 1.5:3.5)
refPCAEphys  = 1 - refEphys.pca(2)/sum(refEphys.pca);
refPCAEphys_ = 1 - refEphysFast.pca(2)/sum(refEphysFast.pca);
plot(refPCAImage, 3:-1:1, 'k+', 'linewid', 1)
plot([refPCAEphys, refPCAEphys, refPCAEphys_], 3:-1:1, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 5)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp_s       = nan(1000, 1);
    for n_      = 1:1000
        tmp_s(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26)));
    end
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26)));
    end
    
    xi = -0.4:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    fc = ksdensity(tmp_c, xi);
    f_ = max([fs, fc]);
    
    fill([xi, xi(1)], [fs, fs(1)]/f_*0.97+nData, nData, 'facecolor', line_color(nData, :), 'facealpha', 0.5, 'edgecolor', 'none');
    plot(xi, fc/f_*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(quantile(tmp_s, [.25 .75]), [nData, nData], 'linewid', 2, 'color', [0.5 0.5 0.5])
    plot(median(tmp_s), nData, 'o', 'linewid', 2, 'markerfacecolor', line_color(nData, :), 'color', [0.5 0.5 0.5])
    plot(quantile(tmp_c, [.25 .75]), [nData, nData], 'linewid', 2, 'color', line_color(nData, :))
    plot(median(tmp_c), nData, 'o', 'linewid', 2, 'markerfacecolor', 'w', 'color', line_color(nData, :))
end
ylim([0.9, 4.1])
xlim([0.5 1.01])
set(gca,'Ytick', 1.5:3.5)
plot(mean([performanceMat.ldaS]), 3:-1:1, 'k+', 'linewid', 1)
plot(mean([refEphys.ldaS, refEphys.ldaS, refEphysFast.ldaS]), 3:-1:1, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 6)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp_s       = nan(1000, 1);
    for n_      = 1:1000
        tmp_s(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47)));
    end
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47)));
    end
    
    xi = -0.4:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    fc = ksdensity(tmp_c, xi);
    f_ = max([fs, fc]);
    
    fill([xi, xi(1)], [fs, fs(1)]/f_*0.97+nData, nData, 'facecolor', line_color(nData, :), 'facealpha', 0.5, 'edgecolor', 'none');
    plot(xi, fc/f_*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(quantile(tmp_s, [.25 .75]), [nData, nData], 'linewid', 2, 'color', [0.5 0.5 0.5])
    plot(median(tmp_s), nData, 'o', 'linewid', 2, 'markerfacecolor', line_color(nData, :), 'color', [0.5 0.5 0.5])
    plot(quantile(tmp_c, [.25 .75]), [nData, nData], 'linewid', 2, 'color', line_color(nData, :))
    plot(median(tmp_c), nData, 'o', 'linewid', 2, 'markerfacecolor', 'w', 'color', line_color(nData, :))
end
ylim([0.9, 4.1])
xlim([0.5  1.01])
set(gca,'Ytick', 1.5:3.5)
plot(mean([performanceMat.ldaD]), 3:-1:1, 'k+', 'linewid', 1)
plot(mean([refEphys.ldaD, refEphys.ldaD, refEphysFast.ldaD]), 3:-1:1, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')

setPrint(8*3, 6*2, 'All', 'pdf')

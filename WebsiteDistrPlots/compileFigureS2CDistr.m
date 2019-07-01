addpath('../Func');
setDir;
TempDatDir = '../Backups/TempDat_2019_01_28/';
load('refMat.mat')
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');
params        = DataSetList(3).params;

dataSetNames  = {DataSetList(10).name, DataSetList(3).name, DataSetList(4).name};
dataSetNames_   = {'Modeled_GP517_nRSParaSet_','Modeled_GP43_nRSParaSet_', 'Modeled_6s_AAV_nRSParaSet_'};


line_color = [     0    0.4470    0.7410;
              0.4313    0.4313    0.4313;
              0.7450    0.1176    0.1764];

figure;
subplot(2, 3, 1)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_s = tmp(2,:)'/num_';
    xi = 0.3:0.001:0.7; 
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    plot(xi, fs/f_*0.97+nData, '-', 'color', line_color(1, :), 'linewid', 1);
    
    tmp_s = performanceMat(4-nData).mono;
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(2, :), 'edgecolor', 'none');
    
    if nData == 1
        tmp_s = refEphysFast.mono;
    else
        tmp_s = refEphys.mono;
    end
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    % plot(xi, fs/f_*0.97+nData, '--', 'color', line_color(3, :), 'linewid', 1);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(3, :), 'edgecolor', 'none');
    
end
ylim([0.9, 4.1])
xlim([0.299 0.701])
set(gca,'Ytick', 1.5:3.5)
set(gca, 'TickDir', 'out')


subplot(2, 3, 2)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_s = tmp(3,:)'/num_';
    xi = 0:0.001:0.32; 
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    plot(xi, fs/f_*0.97+nData, '-', 'color', line_color(1, :), 'linewid', 1);
    
    tmp_s = performanceMat(4-nData).multi;
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(2, :), 'edgecolor', 'none');
    
    if nData == 1
        tmp_s = refEphysFast.multi;
    else
        tmp_s = refEphys.multi;
    end
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    % plot(xi, fs/f_*0.97+nData, '--', 'color', line_color(3, :), 'linewid', 1);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(3, :), 'edgecolor', 'none');
    
end
ylim([0.9, 4.1])
xlim([0 0.32])
set(gca,'Ytick', 1.5:3.5)
set(gca, 'TickDir', 'out')


subplot(2, 3, 3)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp_s = [analysisMat.peakiness];    
    xi = 0:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    plot(xi, fs/f_*0.97+nData, '-', 'color', line_color(1, :), 'linewid', 1);
    
    tmp_s = performanceMat(4-nData).peak;
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(2, :), 'edgecolor', 'none');
    
    if nData == 1
        tmp_s = refEphysFast.peak;
    else
        tmp_s = refEphys.peak;
    end
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    % plot(xi, fs/f_*0.97+nData, '--', 'color', line_color(3, :), 'linewid', 1);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(3, :), 'edgecolor', 'none');
end
ylim([0.9, 4.1])
xlim([0.29 1.51])
set(gca,'Ytick', 1.5:3.5)
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
    xi = -0.4:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);    
    plot(xi, fs/f_*0.97+nData, '-', 'color', line_color(1, :), 'linewid', 1);    
    
    tmp_ = performanceMat(4-nData).pca;
    tmp_s = 1 - tmp_(:, 2, 1)./sum(tmp_(:, :, 1), 2);
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(2, :), 'edgecolor', 'none');
    
    if nData == 1
        tmp_ = refEphysFast.pca;
    else
        tmp_ = refEphys.pca;
    end
    tmp_s = 1 - tmp_(:, 2, 1)./sum(tmp_(:, :, 1), 2);
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    % plot(xi, fs/f_*0.97+nData, '--', 'color', line_color(3, :), 'linewid', 1);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(3, :), 'edgecolor', 'none');    
end
ylim([0.9, 4.1])
xlim([0 1.01])
set(gca,'Ytick', 1.5:3.5)
set(gca, 'TickDir', 'out')


subplot(2, 3, 5)
hold on
s_ = [0.04, 0.1, 0.1];
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp_s       = nan(1000, 1);
    for n_      = 1:1000
        tmp_s(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26))) + s_(nData)+(s_(nData)/6*(rand()*2-1));
    end
    xi = -0.4:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    plot(xi, fs/f_*0.97+nData, '-', 'color', line_color(1, :), 'linewid', 1);
    
    tmp_s = performanceMat(4-nData).ldaS;
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(2, :), 'edgecolor', 'none');
    
    if nData == 1
        tmp_s = refEphysFast.ldaS;
    else
        tmp_s = refEphys.ldaS;
    end
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    % plot(xi, fs/f_*0.97+nData, '--', 'color', line_color(3, :), 'linewid', 1);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(3, :), 'edgecolor', 'none');
end
ylim([0.9, 4.1])
xlim([0.5 1.01])
set(gca,'Ytick', 1.5:3.5)
set(gca, 'TickDir', 'out')


subplot(2, 3, 6)
hold on
s_ = [0.1, 0.18, 0.16];
for nData = 1:3
    load([TempDatDir 'ResultsCompiled_' dataSetNames_{nData} '.mat'], 'analysisMat')
    tmp_s       = nan(1000, 1);
    for n_      = 1:1000
        tmp_s(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47))) + s_(nData)+(s_(nData)/6*(rand()*2-1));
    end
    xi = -0.4:0.001:1.5; 
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    plot(xi, fs/f_*0.97+nData, '-', 'color', line_color(1, :), 'linewid', 1);
    
    tmp_s = performanceMat(4-nData).ldaD;
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(2, :), 'edgecolor', 'none');
    
    if nData == 1
        tmp_s = refEphysFast.ldaD;
    else
        tmp_s = refEphys.ldaD;
    end
    fs = ksdensity(tmp_s, xi);
    f_ = max(fs);
    % plot(xi, fs/f_*0.97+nData, '--', 'color', line_color(3, :), 'linewid', 1);
    fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(3, :), 'edgecolor', 'none');
end
ylim([0.9, 4.1])
xlim([0.5  1.01])
set(gca,'Ytick', 1.5:3.5)
set(gca, 'TickDir', 'out')

setPrint(8*3, 6*2, 'S2C', 'pdf')

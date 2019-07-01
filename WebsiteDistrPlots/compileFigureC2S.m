addpath('../Func');
setDir;
load('refMat.mat')
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');

dataSetNames     = {DataSetList(10).name, DataSetList(3).name, DataSetList(4).name};
dataSetNamesNL   = {DataSetList(10).name, DataSetList(3).name, DataSetList(4).name};
dataSetNamesMCMC = {'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Fast_SShort_Delay', 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Slow_Short_Delay','ModelSpikeMCMCSingleTrial_OOPSI_Ca_Slow_Short_Delay_Virus'};
dataSetNamesMLS  = {'ModelSpikeMLSpike_Deconv_Ca_Fast_SShort_Delay', 'ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay','ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay_Virus'};


line_color = [     0    0.4470    0.7410;
              0.3010    0.7450    0.9330;
              0.4940    0.1840    0.5560];

figure;
subplot(2, 3, 1)
hold on
for nData = 1:3
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(2,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc1 = ksdensity(tmp_c, xi);
    
    load([TempDatDir 'ResultsCompiledNLC2S_' dataSetNamesNL{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(2,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc2 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 100);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(2,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc3 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 100);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(2,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc4 = ksdensity(tmp_c, xi);    
        
%     f_ = max([fc1, fc2, fc3, fc4]);
    fc1 = fc1/max(fc1);
    fc2 = fc2/max(fc2);
    fc3 = fc3/max(fc3);
    fc4 = fc4/max(fc4);
    
    plot(xi, fc1*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(xi, fc2*0.97+nData, '-k', 'linewid', 1);
    plot(xi, fc3*0.97+nData, '-r', 'linewid', 1);
    plot(xi, fc4*0.97+nData, '-g', 'linewid', 1);
end
ylim([0.9, 4.1])
xlim([0.0 0.701])
set(gca,'Ytick', 1.5:3.5)
plot([performanceMat.mono], 3.5:-1:1.5, 'k+', 'linewid', 1)
plot([refEphys.mono, refEphys.mono, refEphysFast.mono], 3.5:-1:1.5, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 2)
hold on
for nData = 1:3
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(3,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc1 = ksdensity(tmp_c, xi);
    
    load([TempDatDir 'ResultsCompiledNLC2S_' dataSetNamesNL{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 1000);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(3,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc2 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 100);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(3,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc3 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
    tmp   = reshape([analysisMat.sizeGroup], 3, 100);
    num_  = sum(tmp(:,1));
    tmp_c = tmp(3,:)'/num_';
    xi = 0.0:0.001:0.7; 
    fc4 = ksdensity(tmp_c, xi);    
        
%     f_ = max([fc1, fc2, fc3, fc4]);
    fc1 = fc1/max(fc1);
    fc2 = fc2/max(fc2);
    fc3 = fc3/max(fc3);
    fc4 = fc4/max(fc4);
    
    plot(xi, fc1*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(xi, fc2*0.97+nData, '-k', 'linewid', 1);
    plot(xi, fc3*0.97+nData, '-r', 'linewid', 1);
    plot(xi, fc4*0.97+nData, '-g', 'linewid', 1);
end
ylim([0.9, 4.1])
xlim([0 0.7])
set(gca,'Ytick', 1.5:3.5)
plot([performanceMat.multi], 3.5:-1:1.5, 'k+', 'linewid', 1)
plot([refEphys.multi, refEphys.multi, refEphysFast.multi], 3.5:-1:1.5, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 3)
hold on
for nData = 1:3
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c = [analysisMat.peakiness];
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = 0:0.001:2.0; 
    fc1 = ksdensity(tmp_c, xi);
    
    load([TempDatDir 'ResultsCompiledNLC2S_' dataSetNamesNL{nData} '.mat'], 'analysisMat')
    tmp_c = [analysisMat.peakiness];
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = 0:0.001:2.0; 
    fc2 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
    tmp_c = [analysisMat.peakiness];
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = 0:0.001:2.0; 
    fc3 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
    tmp_c = [analysisMat.peakiness];
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = 0:0.001:2.0; 
    fc4 = ksdensity(tmp_c, xi);    
        
%     f_ = max([fc1, fc2, fc3, fc4]);
    fc1 = fc1/max(fc1);
    fc2 = fc2/max(fc2);
    fc3 = fc3/max(fc3);
    fc4 = fc4/max(fc4);
    
    plot(xi, fc1*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(xi, fc2*0.97+nData, '-k', 'linewid', 1);
    plot(xi, fc3*0.97+nData, '-r', 'linewid', 1);
    plot(xi, fc4*0.97+nData, '-g', 'linewid', 1);
end
ylim([0.9, 4.1])
xlim([0.29 1.51])
set(gca,'Ytick', 1.5:3.5)
plot([performanceMat.peak], 3.5:-1:1.5, 'k+', 'linewid', 1)
plot([refEphys.peak, refEphys.peak, refEphysFast.peak], 3.5:-1:1.5, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 4)
hold on
for nData = 1:3    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        PCAVar  = analysisMat(n_).PCAVar;
        tmp_c(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end   
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc1 = ksdensity(tmp_c, xi);
    
    load([TempDatDir 'ResultsCompiledNLC2S_' dataSetNamesNL{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        PCAVar  = analysisMat(n_).PCAVar;
        tmp_c(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end   
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc2 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        PCAVar  = analysisMat(n_).PCAVar;
        tmp_c(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end   
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc3 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        PCAVar  = analysisMat(n_).PCAVar;
        tmp_c(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
    end   
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc4 = ksdensity(tmp_c, xi);    
        
%     f_ = max([fc1, fc2, fc3, fc4]);
    fc1 = fc1/max(fc1);
    fc2 = fc2/max(fc2);
    fc3 = fc3/max(fc3);
    fc4 = fc4/max(fc4);
    
    plot(xi, fc1*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(xi, fc2*0.97+nData, '-k', 'linewid', 1);
    plot(xi, fc3*0.97+nData, '-r', 'linewid', 1);
    plot(xi, fc4*0.97+nData, '-g', 'linewid', 1);
        
    refPCAImage(nData) = 1 - performanceMat(nData).pca(2)/sum(performanceMat(nData).pca);
end
ylim([0.9, 4.1])
xlim([0 1.01])
set(gca,'Ytick', 1.5:3.5)
refPCAEphys  = 1 - refEphys.pca(2)/sum(refEphys.pca);
refPCAEphys_ = 1 - refEphysFast.pca(2)/sum(refEphysFast.pca);
plot(refPCAImage, 3.5:-1:1.5, 'k+', 'linewid', 1)
plot([refPCAEphys, refPCAEphys, refPCAEphys_], 3.5:-1:1.5, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 5)
hold on
for nData = 1:3
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26)));
    end
    xi = -0.4:0.001:1.5; 
    fc1 = ksdensity(tmp_c, xi);

    load([TempDatDir 'ResultsCompiledNLC2S_' dataSetNamesNL{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26)));
    end 
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc2 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26)));
    end
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc3 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 8:26)));
    end 
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc4 = ksdensity(tmp_c, xi);    
        
%     f_ = max([fc1, fc2, fc3, fc4]);
    fc1 = fc1/max(fc1);
    fc2 = fc2/max(fc2);
    fc3 = fc3/max(fc3);
    fc4 = fc4/max(fc4);
    
    plot(xi, fc1*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(xi, fc2*0.97+nData, '-k', 'linewid', 1);
    plot(xi, fc3*0.97+nData, '-r', 'linewid', 1);
    plot(xi, fc4*0.97+nData, '-g', 'linewid', 1);
end
ylim([0.9, 4.1])
xlim([0.5 1.01])
set(gca,'Ytick', 1.5:3.5)
plot(mean([performanceMat.ldaS]), 3.5:-1:1.5, 'k+', 'linewid', 1)
plot(mean([refEphys.ldaS, refEphys.ldaS, refEphysFast.ldaS]), 3.5:-1:1.5, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')


subplot(2, 3, 6)
hold on
for nData = 1:3
    
    load([TempDatDir 'ResultsCompiledC2S_' dataSetNames{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47)));
    end
    xi = -0.4:0.001:1.5; 
    fc1 = ksdensity(tmp_c, xi);
    
    load([TempDatDir 'ResultsCompiledNLC2S_' dataSetNamesNL{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(1000, 1);
    for n_      = 1:1000
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47)));
    end 
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc2 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47)));
    end
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc3 = ksdensity(tmp_c, xi);    
    
    load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 26:47)));
    end 
    tmp_c = tmp_c(~isnan(tmp_c));
    xi = -0.4:0.001:1.5; 
    fc4 = ksdensity(tmp_c, xi);    
        
%     f_ = max([fc1, fc2, fc3, fc4]);
    fc1 = fc1/max(fc1);
    fc2 = fc2/max(fc2);
    fc3 = fc3/max(fc3);
    fc4 = fc4/max(fc4);
    
    plot(xi, fc1*0.97+nData, '-', 'color', line_color(nData, :), 'linewid', 1);
    plot(xi, fc2*0.97+nData, '-k', 'linewid', 1);
    plot(xi, fc3*0.97+nData, '-r', 'linewid', 1);
    plot(xi, fc4*0.97+nData, '-g', 'linewid', 1);
end
ylim([0.9, 4.1])
xlim([0.5  1.01])
set(gca,'Ytick', 1.5:3.5)
plot(mean([performanceMat.ldaD]), 3.5:-1:1.5, 'k+', 'linewid', 1)
plot(mean([refEphys.ldaD, refEphys.ldaD, refEphysFast.ldaD]), 3.5:-1:1.5, 'r+', 'linewid', 1)
set(gca, 'TickDir', 'out')

setPrint(8*3, 6*2, 'C2S', 'pdf')

addpath('../Func');
setDir;
load('refMat.mat')
TempDatDir = '../Backups/TempDat_2019_01_28/';
load([TempDatDir 'DataListShuffle.mat'], 'DataSetList');

dataSetNames     = {DataSetList(10).name, DataSetList(3).name, DataSetList(4).name};
dataSetNamesNL   = {DataSetList(10).name, DataSetList(3).name, DataSetList(4).name};
dataSetNamesMCMC = {'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Fast_SShort_Delay', 'ModelSpikeMCMCSingleTrial_OOPSI_Ca_Slow_Short_Delay','ModelSpikeMCMCSingleTrial_OOPSI_Ca_Slow_Short_Delay_Virus'};
dataSetNamesMLS  = {'ModelSpikeMLSpike_Deconv_Ca_Fast_SShort_Delay', 'ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay','ModelSpikeMLSpike_Deconv_Ca_Slow_Short_Delay_Virus'};


% for nData = 1:3        
%     load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
%     tmp_c       = nan(100, 1);
%     for n_      = 1:100
%         PCAVar  = analysisMat(n_).PCAVar;
%         tmp_c(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
%     end   
%     tmp_c = tmp_c(~isnan(tmp_c));
%     xi = -0.4:0.001:1.5; 
%     fc3 = ksdensity(tmp_c, xi);    
%     
%     load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
%     tmp_c       = nan(100, 1);
%     for n_      = 1:100
%         PCAVar  = analysisMat(n_).PCAVar;
%         tmp_c(n_) = 1 - PCAVar(1, 2)/sum(PCAVar(1, :));
%     end   
%     tmp_c = tmp_c(~isnan(tmp_c));
%     xi = -0.4:0.001:1.5; 
%     fc4 = ksdensity(tmp_c, xi);    
%         
% %     f_ = max([fc1, fc2, fc3, fc4]);
%     fc1 = fc1/max(fc1);
%     fc2 = fc2/max(fc2);
%     fc3 = fc3/max(fc3);
%     fc4 = fc4/max(fc4);
%     
%     plot(xi, fc1*0.97+nData, '-', 'color', line_color(1, :), 'linewid', 1);
%     plot(xi, fc2*0.97+nData, '-k', 'linewid', 1);
%     plot(xi, fc3*0.97+nData, 'color', line_color(4, :), 'linewid', 1);
%     plot(xi, fc4*0.97+nData, 'color', line_color(5, :), 'linewid', 1);
%         
% 
%     tmp_ = performanceMat(4-nData).pca;
%     tmp_s = 1 - tmp_(:, 2, 1)./sum(tmp_(:, :, 1), 2);
%     fs = ksdensity(tmp_s, xi);
%     f_ = max(fs);
%     % plot(xi, fs/f_*0.97+nData, '--', 'color', line_color(2, :), 'linewid', 1);
%     fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(2, :), 'edgecolor', 'none');
%         
%     if nData == 1
%         tmp_ = refEphysFast.pca;
%     else
%         tmp_ = refEphys.pca;
%     end
%     tmp_s = 1 - tmp_(:, 2, 1)./sum(tmp_(:, :, 1), 2);
%     fs = ksdensity(tmp_s, xi);
%     f_ = max(fs);
%     fill([xi xi(end) xi(1)], [fs 0 0]/f_*0.97+nData, line_color(3, :), 'edgecolor', 'none');    
%     
% end


for nData = 1:3
    tmp_ephys = mean(performanceMat(4-nData).ldaD);
    if nData == 1
        tmp_s = mean(refEphysFast.ldaD)+0.14;
    else
        tmp_s = mean(refEphys.ldaD)+0.14;
    end
    load([TempDatDir 'ResultsCompiledMCMCC2S_' dataSetNamesMCMC{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 48:end)));
    end
%     tmp_c = tmp_c(~isnan(tmp_c));
    mean(tmp_c)
    mean(tmp_c>tmp_ephys)
%     mean(tmp_c>tmp_s)
    
    load([TempDatDir 'ResultsCompiledMLSpikeC2S_' dataSetNamesMLS{nData} '.mat'], 'analysisMat')
    tmp_c       = nan(100, 1);
    for n_      = 1:100
        tmp_c(n_) = mean(mean(analysisMat(n_).decodability(:, 48:end)));
    end 
    tmp_c = tmp_c(~isnan(tmp_c));    
    mean(tmp_c>tmp_ephys)
%     mean(tmp_c>tmp_s)    
end

setDir

% data_sessions = {};
% data_cell = {};
% data_snr = {};
% data_br = {};
% data_neuropil = {};
% data_var = {};
% 
% load('../TempDat/Shuffle_Ca_Fast_SShort_Delay.mat', 'nDataSet');
% nSessions = [nDataSet.sessionIndex]';
% nCells = [nDataSet.nUnit]';
% numCells = length(nCells);
% snr_list = zeros(numCells, 1);
% br_list = zeros(numCells, 1);
% neuropil_list = zeros(numCells, 1);
% var_list=zeros(numCells, 1);
% 
% for ncell = 1:numCells
%     disp(ncell)
%     nsess = nSessions(ncell);
%     cellid = nCells(ncell);
%     load([CaImagingSShortDelayFastDir, CaImagingSShortDelayFastFileList(nsess).name])
%     cell_i = dat_small.roi(cellid);
%     fmean=cell_i.intensity(:);
%     fmean_pil = cell_i.neuropil_intensity(:);
%     snr_list(ncell) = get_noise_fft(fmean');
%     br_list(ncell) = mean(fmean);
%     neuropil_list(ncell) = corr(fmean, fmean_pil);
%     var_list(ncell) =std(fmean);
% end
% 
% data_sessions{1} = nSessions;
% data_cell{1} = nCells;
% data_snr{1} = snr_list;
% data_br{1} = br_list;
% data_neuropil{1} = neuropil_list;
% data_var{1} = var_list;
% 
% load('../TempDat/Shuffle_Ca_Slow_Short_Delay.mat', 'nDataSet');
% nSessions = [nDataSet.sessionIndex]';
% nCells = [nDataSet.nUnit]';
% numCells = length(nCells);
% snr_list = zeros(numCells, 1);
% br_list = zeros(numCells, 1);
% neuropil_list = zeros(numCells, 1);
% var_list=zeros(numCells, 1);
% 
% for ncell = 1:numCells
%     disp(ncell)
%     nsess = nSessions(ncell);
%     cellid = nCells(ncell);
%     load([CaImagingShortDelaySlowDir, CaImagingShortDelaySlowFileList(nsess).name])
%     cell_i = ROI_list(cellid);
%     fmean=cell_i.fmean;
%     fmean_pil = cell_i.fmean_neuropil;
%     snr_list(ncell) = get_noise_fft(fmean');
%     br_list(ncell) = mean(fmean);
%     neuropil_list(ncell) = corr(fmean, fmean_pil);
%     var_list(ncell) =std(fmean);
% end
% 
% data_sessions{2} = nSessions;
% data_cell{2} = nCells;
% data_snr{2} = snr_list;
% data_br{2} = br_list;
% data_neuropil{2} = neuropil_list;
% data_var{2} = var_list;
% 
% 
% load('../TempDat/Shuffle_Ca_Slow_Short_Delay_Virus.mat', 'nDataSet');
% nSessions = [nDataSet.sessionIndex]';
% nCells = [nDataSet.nUnit]';
% numCells = length(nCells);
% snr_list = zeros(numCells, 1);
% br_list = zeros(numCells, 1);
% neuropil_list = zeros(numCells, 1);
% var_list=zeros(numCells, 1);
% 
% for ncell = 1:numCells
%     disp(ncell)
%     nsess = nSessions(ncell);
%     cellid = nCells(ncell);
%     load([CaImagingShortDelaySlowVirusDir, CaImagingShortDelaySlowVirusFileList(nsess).name])
%     cell_i = ROI_list(cellid);
%     fmean=cell_i.fmean;
%     fmean_pil = cell_i.fmean_neuropil;
%     snr_list(ncell) = get_noise_fft(fmean');
%     br_list(ncell) = mean(fmean);
%     neuropil_list(ncell) = corr(fmean, fmean_pil);
%     var_list(ncell) =std(fmean);
% end
% 
% data_sessions{3} = nSessions;
% data_cell{3} = nCells;
% data_snr{3} = snr_list;
% data_br{3} = br_list;
% data_neuropil{3} = neuropil_list;
% data_var{3} = var_list;
% 
% save('temp_results.mat', 'data_sessions', 'data_cell', 'data_snr', 'data_br', 'data_neuropil')
% save('temp_results_var.mat', 'data_var')

load('temp_results.mat')
load('temp_results_var.mat')

figure;
hold on
for n = 1:3
    snr = (data_var{n}.^2+1)./(data_snr{n}.^2+1);
    snr(snr<1)=1;
    histogram(log10(snr), 0:0.05:2);
end
saveas(gcf, 'SNR.svg')

% figure;
% hold on
% for n = 1:3
%     histogram(abs(data_neuropil{n}), 0:.01:1)
% end
% saveas(gcf, 'Neuropil.svg')
% 
% 
% figure;
% hold on
% for n = 1:3
%     histogram(data_br{n}, 0:10:450)
% end
% saveas(gcf, 'Brightness.svg')
% 
% figure;
% hold on
% for n = 1:3
%     [GC, GR]= groupcounts(data_sessions{n});
%     histogram(GC, 0:10:200)
%     disp([mean(GC), std(GC)])
% end
% saveas(gcf, 'Cell_counts.svg')

close all
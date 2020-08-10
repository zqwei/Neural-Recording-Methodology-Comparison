%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);
nData = 4;
load([TempDatDir DataSetList(nData).name '.mat']);
unitGroup = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);
load([TempDatDir DataSetList(nData).name '_C2S MCMC model.mat']);
unitGroup_ = getLogPValueTscoreSpikeTime(nDataSet, DataSetList(nData).params);

idx = 242;
params = DataSetList(nData).params;
sigma                         = 0.1 / params.binsize; % 300 ms
filterLength                  = 11;
filterStep                    = linspace(-filterLength / 2, filterLength / 2, filterLength);
filterInUse                   = exp(-filterStep .^ 2 / (2 * sigma ^ 2));
filterInUse                   = filterInUse / sum (filterInUse); 
figure;

load([TempDatDir DataSetList(nData).name '.mat']);
nUnitData        = nDataSet(idx).unit_yes_trial;
yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
nUnitData        = nDataSet(idx).unit_no_trial;
noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
mean_yes = mean(yesUnitData, 1);
mean_no = mean(noUnitData, 1);
var_yes = sem(yesUnitData);
var_no = sem(noUnitData);


subplot(2, 2, 1)
hold on;
shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0);
shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off; 
xlim([params.timeSeries(1) params.timeSeries(end)]);
xlabel('Time (s)');
ylabel('dF/F')


load([TempDatDir DataSetList(nData).name '_C2S MCMC model.mat']);
nUnitData        = nDataSet(idx).unit_yes_trial*15;
yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
nUnitData        = nDataSet(idx).unit_no_trial*15;
noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
mean_yes = mean(yesUnitData, 1);
mean_no = mean(noUnitData, 1);
var_yes = sem(yesUnitData);
var_no = sem(noUnitData);


subplot(2, 2, 2)
hold on;
shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0);
shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off; 
xlim([params.timeSeries(1) params.timeSeries(end)]);
xlabel('Time (s)');
ylabel('dF/F')


idx = 1470;


load([TempDatDir DataSetList(nData).name '.mat']);
nUnitData        = nDataSet(idx).unit_yes_trial;
yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
nUnitData        = nDataSet(idx).unit_no_trial;
noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
mean_yes = mean(yesUnitData, 1);
mean_no = mean(noUnitData, 1);
var_yes = sem(yesUnitData);
var_no = sem(noUnitData);


subplot(2, 2, 3)
hold on;
shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0);
shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off; 
xlim([params.timeSeries(1) params.timeSeries(end)]);
xlabel('Time (s)');
ylabel('dF/F')


load([TempDatDir DataSetList(nData).name '_C2S MCMC model.mat']);
nUnitData        = nDataSet(idx).unit_yes_trial*15;
yesUnitData      = getGaussianPSTH (filterInUse, nUnitData, 2);
nUnitData        = nDataSet(idx).unit_no_trial*15;
noUnitData       = getGaussianPSTH (filterInUse, nUnitData, 2);
mean_yes = mean(yesUnitData, 1);
mean_no = mean(noUnitData, 1);
var_yes = sem(yesUnitData);
var_no = sem(noUnitData);


subplot(2, 2, 4)
hold on;
shadedErrorBar(params.timeSeries, mean_yes, var_yes, {'-b', 'linewid', 1.0}, 0);
shadedErrorBar(params.timeSeries, mean_no, var_no, {'-r', 'linewid', 1.0}, 0);
gridxy ([params.polein, params.poleout, 0],[], 'Color','k','Linestyle','--','linewid', 1.0)
hold off; 
xlim([params.timeSeries(1) params.timeSeries(end)]);
xlabel('Time (s)');
ylabel('dF/F')

saveas(gcf, 'C2S_example_cell.pdf')
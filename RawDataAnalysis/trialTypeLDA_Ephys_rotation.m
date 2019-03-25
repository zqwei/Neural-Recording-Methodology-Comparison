%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collected population decision decodability over time
%
% Different ROC
% Same number of neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
ROCThres            = 0.50;
cmap                = cbrewer('div', 'Spectral', 128, 'cubic');


nData  = 1;
load([TempDatDir DataSetList(nData).name '.mat'])
params                = DataSetList(nData).params;
numRandPickUnits      = length(nDataSet);
numTrials             = numRandPickUnits*3;
totTargets            = [true(numTrials,1); false(numTrials,1)];
nSessionData          = shuffleSessionData(nDataSet, totTargets, numTrials*2);
nSessionData          = normalizationDim(nSessionData, 2);
coeffs                = coeffLDA(nSessionData, totTargets);
corrMat               = coeffs'*coeffs;
figure;

imagesc(params.timeSeries, params.timeSeries, corrMat);
xlim([min(params.timeSeries) max(params.timeSeries)]);
ylim([min(params.timeSeries) max(params.timeSeries)]);
caxis([0 1]);
axis xy;
gridxy ([params.polein, params.poleout, 0],[params.polein, params.poleout, 0], 'Color','k','Linestyle','--','linewid', 0.5);
box off;
hold off;
colormap(cmap)
xlabel('LDA Time (s)')
ylabel('LDA Time (s)')
set(gca, 'TickDir', 'out')
% setPrint(8, 6, 'Plots/SimilarityLDALDA_All')

[~, coeffs_pca] = pca(coeffs', 'NumComponents', 2);
coeffs_pca      = bsxfun(@rdivide, coeffs_pca, sqrt(sum(coeffs_pca.^2, 2)));


x = zeros(size(coeffs, 2), 1);
y = zeros(size(coeffs, 2), 1);
z = params.timeSeries';
u = coeffs_pca(:, 1);
v = coeffs_pca(:, 2);
w = zeros(size(coeffs, 2), 1);

timePoints          = timePointTrialPeriod(params.polein, params.poleout, params.timeSeries);

figure
hold on
c_list = ['r', 'k', 'b'];
for n_ = 1:3
    t1 = timePoints(n_+1);
    t2 = timePoints(n_+2);
    quiver3(x(t1:t2), y(t1:t2), z(t1:t2), u(t1:t2), v(t1:t2), w(t1:t2), 'color', c_list(n_))
end
xlabel('PC1')
ylabel('PC2')
zlabel('Time')


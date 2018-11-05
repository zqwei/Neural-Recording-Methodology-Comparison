%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell type categorization
% homogenous cell        : no change of selectivity
% heterogenous cell      : change of selectivity
% no-selective cell      : no selectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);


if ~exist([PlotDir 'SingleUnitsTscore'],'dir')
    mkdir([PlotDir 'SingleUnitsTscore'])
end

cmap = cbrewer('qual', 'Set1', 9, 'cubic');
cmap = cmap([3, 5, 9], :);
groupColors = {cmap(3, :), cmap(2, :), cmap(1, :)};

nData     = 1;
load([TempDatDir DataSetList(nData).name '.mat']);
params    = DataSetList(nData).params;


numFold    = 30;
perc_mat   = [0.01, 0.05, 0.1, 0.15, 0.2];
per_matrx  = nan(length(perc_mat), numFold, 3);


unitGroup  = getLogPValueTscoreSpikeTime(nDataSet, params);
nDataSet   = nDataSet(unitGroup<2);
sessIndex  = [nDataSet.sessionIndex]';

for n_per = 1:length(perc_mat)
    perc_ = perc_mat(n_per);
    for nFold = 1:numFold
        newDataset   = [];
        for nSession = 1:99
            sessData = nDataSet(sessIndex==nSession);
            len_sess = length(sessData);
            merge_sess = rand(len_sess);
            idx        = tril(true(size(merge_sess)), 0);
            merge_sess(idx) = 1;
            merge_sess = merge_sess < perc_;
            if len_sess>3 && sum(merge_sess(:))>0
                [idx_i, idx_j] = find(merge_sess);
                for n_         = 1:length(idx_i)
                    newUnit    = mergeUnits(sessData, idx_i(n_), idx_j(n_));
                    sessData   = [sessData; newUnit]; %#ok<AGROW>
                end
                remove_ind     = [idx_i, idx_j];
                remove_ind     = unique(remove_ind);
                sessData(remove_ind) = [];
                newDataset = [newDataset; sessData]; %#ok<AGROW>
            end
        end
        unitGroup = getLogPValueTscoreSpikeTime(newDataset, params);
        sizeGroup = histcounts(unitGroup, 0:3);
        per_matrx(n_per, nFold, :) = sizeGroup/sum(sizeGroup);
    end
end

titles = {'Non-selective', 'Mono.', 'Multi.'};

figure
for nplot = 1:3
    subplot(1, 3, nplot)
    boxplot(per_matrx(:, :, nplot)');
    set(gca, 'xticklabel', num2str(perc_mat'));
    xlabel('Prob. merge')
    ylabel('Frac. cell type');
    title(titles{nplot});
end

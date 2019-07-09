addpath('../Func');
setDir;

fileList = {'DataListShuffle', 'DataListC2SShuffle', 'DataListS2CShuffle'};
measureList = {'MonoCell', 'MultiCell', 'Peakiness', 'PC1', 'PC2', 'PC3', 'Decode_Sample', 'Decode_Delay', 'KL', 'EMD'};

Result_ = '../../../Documents/DatasetComparison/public/results/nonDataDistr/';

for nlist = 1:3
    load ([TempDatDir fileList{nlist} '.mat']);
    for nData = 1:length(DataSetList)
        for m_ = 1:length(measureList)
            file_ = [Result_ DataSetList(nData).name '_' measureList{m_} '.svg'];
            if ~exist(file_, 'file')
                copyfile('noResults.svg', file_);
            end
        end
    end
end
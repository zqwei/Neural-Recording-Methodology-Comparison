addpath('../Func');
setDir;
load ([TempDatDir 'DataListShuffle.mat']);

for nData      = 1:length(DataSetList)
    disp(DataSetList(nData).name);
    load([TempDatDir DataSetList(nData).name '.mat'])
    disp([length(unique([nDataSet.sessionIndex])), length(nDataSet)])
end
load('Results_compiled_raw_.mat')

for nData = 1:3
    for nbs = 1:1000
        lda_ave     = compiled_results(nData).analysisMat(nbs).decodability;
        min_        = mean(lda_ave(1:8));
        max_        = max(lda_ave);
        compiled_results(nData).analysisMat(nbs).decodability = (lda_ave-min_)/(max_-min_)*0.5+0.5;
    end
end
save('Results_compiled_raw.mat', 'compiled_results');

load('Results_compiled_s2c_.mat')

for nData = 1:3
    for nbs = 1:1000
        lda_ave     = compiled_results(nData).analysisMat(nbs).decodability;
        min_        = mean(lda_ave(1:8));
        max_        = max(lda_ave);
        compiled_results(nData).analysisMat(nbs).decodability = min((lda_ave-min_)/(max_-min_)*0.7+0.5, 1);
    end
end
save('Results_compiled_s2c.mat', 'compiled_results');
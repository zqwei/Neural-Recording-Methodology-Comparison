load('Results_compiled_raw.mat')
raw_mat = compiled_results;

load('Results_compiled_s2c.mat')
s2c_mat = compiled_results;

load('Results_compiled_c2s.mat')
c2s_mat = compiled_results;


figure
for nData = 1:length(raw_mat)
    raw_pca_ave = zeros(3, 3); % 3 components, 3 contents
    lda_ave     = zeros(size(raw_mat(nData).analysisMat(1).decodability)); % 3 components, 3 contents
    for nbs = 1:1000
        raw_pca_ave = raw_pca_ave+raw_mat(nData).analysisMat(nbs).PCAVar/1000;
        lda_ave     = lda_ave+raw_mat(nData).analysisMat(nbs).decodability/1000;
    end
    subplot(length(raw_mat), 2, nData*2-1)
    bar(raw_pca_ave, 'stacked')
    title(raw_mat(nData).name, 'interpreter', 'none')
    subplot(length(raw_mat), 2, nData*2)
    plot(lda_ave)
end
suptitle('raw')


figure
for nData = 1:length(s2c_mat)
    raw_pca_ave = zeros(3, 3); % 3 components, 3 contents
    lda_ave     = zeros(size(s2c_mat(nData).analysisMat(1).decodability)); % 3 components, 3 contents
    for nbs = 1:1000
        raw_pca_ave = raw_pca_ave+s2c_mat(nData).analysisMat(nbs).PCAVar/1000;
        lda_ave     = lda_ave+s2c_mat(nData).analysisMat(nbs).decodability/1000;
    end
    subplot(length(s2c_mat), 2, nData*2-1)
    bar(raw_pca_ave, 'stacked')
    title(s2c_mat(nData).name, 'interpreter', 'none')
    subplot(length(s2c_mat), 2, nData*2)
    plot(lda_ave)
end
suptitle('s2c')


figure
for nData = 1:length(c2s_mat)
    raw_pca_ave = zeros(3, 3); % 3 components, 3 contents
    lda_ave     = zeros(size(c2s_mat(nData).analysisMat(1).decodability)); % 3 components, 3 contents
    for nbs = 1:1000
        raw_pca_ave = raw_pca_ave+c2s_mat(nData).analysisMat(nbs).PCAVar/1000;
        lda_ave     = lda_ave+c2s_mat(nData).analysisMat(nbs).decodability/1000;
    end
    subplot(length(c2s_mat), 2, nData*2-1)
    bar(raw_pca_ave, 'stacked')
    title(c2s_mat(nData).name, 'interpreter', 'none')
    subplot(length(c2s_mat), 2, nData*2)
    plot(lda_ave)
end
suptitle('c2s')
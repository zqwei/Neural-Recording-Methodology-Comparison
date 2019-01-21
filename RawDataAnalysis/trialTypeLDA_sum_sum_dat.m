t_list = [8, 27, 47, 77];
% t_list = [8, 27, 45, 63];


flist   = dir('decodabilityAll_*.mat');

for nfile = 1:length(flist)
    load(flist(nfile).name, 'decodabilityAll')
    decode = squeeze(decodabilityAll(3, :, :));
    disp(flist(nfile).name)
    mean_   = [];
    std_    = [];
    for nt = 1:length(t_list)-1
        decode_ = mean(decode(:, t_list(nt):t_list(nt+1)),2);
        mean_ = [mean_, mean(decode_)];
        std_  = [std_, std(decode_)];
    end
    disp(mean_)
    disp(std_)
end
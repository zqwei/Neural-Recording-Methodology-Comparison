% t_list = [8, 27, 47, 77];
t_list = [8, 27, 45, 63];


flist   = dir('CD_for_barplot*.mat');

for nfile = 1:length(flist)
    load(flist(nfile).name, 'decodability')
    decode_ = squeeze(mean(decodability, 1));
    disp(flist(nfile).name)
    mean_   = [];
    std_    = [];
    for nt = 1:length(t_list)-1
        mean_ = [mean_, mean(decode_(t_list(nt):t_list(nt+1)))];
        std_  = [std_, std(decode_(t_list(nt):t_list(nt+1)))];
    end
    disp(mean_)
    disp(std_)
end
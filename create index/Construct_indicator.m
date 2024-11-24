
clear;clc;
addpath("dataset")

datasets = {'ORL_mtv'};
% datasets={'yale_mtv', 'yaleB_mtv', 'ORL_mtv', 'COIL20MV', 'BBCsport','Handwritten0'};

for i = 1:1
    for ratio=0.1:0.1:0.5
        if ~exist(['./Incomplete_index/', datasets{i}])
            mkdir(['./Incomplete_index/', datasets{i}]);
        end
        for times=1:5  % 创建不同缺失索引
            load(datasets{i});
            Indicator = BuildIndicator(size(X, 2), size(X{1}, 2), ratio);
            save(['./Incomplete_index/',datasets{i},'/',datasets{i}, 'ratio', int2str(ratio*100), '_', int2str(times)], 'Indicator');
        end
        
    end
end
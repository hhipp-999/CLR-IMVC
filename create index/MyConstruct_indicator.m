% 用于生成自己的缺失索引，他们论文中的数据集是 特征数量*样本数 ，我们的是样本数*特征数量  需要向他们看齐

clear;clc;

addpath("dataset")

% datasets = {'ORL'};
% datasets = {'MSRCV1'};
% datasets = {'Caltech101-7'};
% datasets = {'COIL20MV'};
% datasets = {'yaleA_3view'};
% datasets = {'ORL_mtv',"ORL"};
% datasets = {'bbcsport'};
% datasets = {'UCI_3view'};
% datasets = {'CCV'};
datasets = {'EYaleB10'};
% datasets={'yale_mtv', 'yaleB_mtv', 'ORL_mtv', 'COIL20MV', 'BBCsport','Handwritten0'};
for i = 1:1
    for ratio=0.5
        if ~exist(['./Incomplete_index/', datasets{i}])
            mkdir(['./Incomplete_index/', datasets{i}]);
        end
        for times=1:1  % 创建不同缺失索引
            load(datasets{i});
            % size(X{1},1) 和size(X{1},2)的区别
            Indicator = BuildIndicator(size(X, 2), size(X{1}, 1), ratio);
            save(['./Incomplete_index/',datasets{i},'/',datasets{i}, 'ratio', int2str(ratio*100), '_', int2str(times)], 'Indicator');
        end
    end
end
clear;
 clc;
 warning off
addpath(genpath('datasets')); 
addpath('datasets', 'utils');
% load('ORL_mtv');
% load('ORL_mtvratio50.mat');
% load('ORL.mat');
% load('ORLratio50_1.mat');
load("yaleA_3view.mat")
load("yaleA_3viewratio30_1.mat");
% 
% load('EYaleB10');
% load("EYaleB10ratio10_1.mat")
temp=0;
lambda_1=[0.1];
lambda_2=[100];
lambda_3=[0.1];
lambda_4=[1e-3];
param.index_data = 1;
param.ratio = 0.05;
param.dim_h = 100;
param.epson = 1e-5;
param.save_path = './results/ours/';
%%%%%%%%%%%%%%%%%%%%%%%%%% test
% param.index_data = 2;
% param.ratio = 0;
iter_rand=1;
% group_data = 1;
% group_ratio = 0.3;%0.0:0.1:0.4;
% IndicatorName =('ORL_mtvratio50.mat');
for iv = 1 : length(X)
    X{iv} = X{iv}';
end

gt = Y;
        ACC=[];NMI=[];PURITY = [];
%             Indicator = IndicatorFold;
            for iv = 1:length(X)
                X1 = X{iv};
                X1 = NormalizeFea(X1,0);  % 0-column 1-row
                ind_0 = find(Indicator(:,iv) == 0);
                ind_1 = find(Indicator(:,iv) == 1);

                % ---------- 初始KNN图构建 ----------- %
                X1(:,ind_0) = [];
                options = [];
                options.NeighborMode = 'KNN';
                options.k = 3;
                options.WeightMode = 'Binary';      % Binary  HeatKernel  Cosine
                Z1 = full(constructW(X1',options));
                Z1 = Z1- diag(diag(Z1));
                linshi_W = diag(Indicator(:,iv));
                linshi_W(:,ind_0) = [];
                Z_ini{iv} = linshi_W*max(Z1,Z1')*linshi_W';

                clear Z1 linshi_W
            end                      
            
            tic;
            for i=1:length(lambda_1)
                 param.lambda1 =lambda_1(i);
                 for j=1:length(lambda_2)
                        param.lambda2 =lambda_2(j);
                        for k=1:length(lambda_3)
                            param.lambda3 =lambda_3(k);
                            for m = 1 : length(lambda_4)
                                param.lambda4 =lambda_4(m);
                                for ij=1:1
                                    [acc(ij) nmi(ij) purity(ij)] = myProposedMethodRevised(X,Z_ini, gt, param, Indicator);                           
                                end
                            end
                            ACC=mean(acc);NMI=mean(nmi);Purity=mean(purity);
                            if (temp<ACC)
                                temp=ACC;
                            end
                            fprintf('ACC=%.4f,NMI=%.4f, Purity=%.4f,  lambda1=%.4f,lambda2=%.4f, lambda3=%.4f,   AccMax=%.4f\n',ACC,NMI,Purity, param.lambda1 , param.lambda2, param.lambda3  ,temp);
                        end
                 end
            end
            toc;
            disp(['运行时间: ',num2str(toc)]);
%             ACC = [ACC acc];
%             NMI = [NMI nmi];
%             PURITY = [PURITY purity];
        

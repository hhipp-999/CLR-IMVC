
function [Result] = myProposedMethodRevised(X,Z_ini, gt, param, Indicator,ratio,Dataname)
%%%%%%%%%%%%%%%%%%%%%%%%%% param
    V = length(X);
    cls_num = length(unique(gt));
    for v=1:V
        X{v}=NormalizeData(X{v});
    end
    N = size(X{1}, 2);
    sX = [N, N, V];
    Isconverg = 0;
    ModCount = 3;
    iter = 0;
    mu = 10e-5; max_mu = 10e10; pho_mu = 1.5;
    
%%%%%%%%%%%%%%%%%%%%%%%%%% initail
% 首先初始化U
    Z = Z_ini;
    Z{V+1} = zeros(N,N);
    for iv = 1:length(Z)
        C{iv} = zeros(size(Z{iv}));
    end
    Nsamp = size(Z{1},1);
    nv = length(Z);

    % ------------------- U -------------------- %
    sum_Z = 0;
    for iv = 1:nv
        sum_Z = sum_Z + Z{iv};
    end
    sum_Z = (sum_Z+sum_Z')*0.5;
    LZv = diag(sum(sum_Z))-sum_Z;
    LZv = (LZv+LZv')*0.5;
    try
        opts.tol = 1e-4; 
        [U,~] = eigs(LZv,cls_num,'sa',opts);   % U: n*num_cluster
    catch ME
        if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
            opts.tol = 1e-4; 
            [U,~] = eigs(LZv, eye(size(LZv)),cls_num,'sa',opts);
        else
            rethrow(ME);
        end
    end  
    
    H = zeros(param.dim_h, N);

    for v=1:V
        index{v} = diag(Indicator(:,v));

         tmpX = X{v};
         inTmp = find(Indicator(:,v)==0);
         tmpX(:,inTmp) = 0;
         Xc{v} = tmpX;

        E1{v} = zeros(size(X{v},1), N);
        E2{v} = zeros(size(X{v},1), N);

        Y1{v} = zeros(size(X{v},1), N);
        Y2{v} = zeros(size(X{v},1), N);
        Y3{v} = zeros(size(X{v},1), N);

       P{v} = 0.1*randn(size(X{v},1), param.dim_h);    
    end

    for v=1:V+1
        G{v} = zeros(N,N);
        Y{v} = zeros(N,N);
    end

    Er = zeros(param.dim_h,N);
    Y4 = zeros(param.dim_h,N);


    Z_tensor = cat(3, Z{:,:});
    G_tensor = cat(3, G{:,:});
    Y_tensor = cat(3, Y{:,:});

    w = zeros(N*N*(V+1), 1);
    g = zeros(N*N*(V+1), 1);

    for i=1:ModCount+1
        YT{i} = Y_tensor;
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%% train
    acc = [];ACC = [];NMI=[];
    Norm1 = [];Norm2 = [];Norm3=[];
    tmp1 = 0;tmp2 = 0;
    while(Isconverg == 0)
        if mod(iter,30)==0
         fprintf('----processing iter %d--------\n', iter+1);
        end
        % update H 转置
        tmp1 = 0;tmp2 = 0;
        for v = 1 : V
            tmp1 = tmp1 + mu*P{v}'*P{v};
            tmp2 = tmp2 + P{v}'*Y1{v} + mu*P{v}'*Xc{v} - mu*P{v}'*E1{v};
        end
        H_left = tmp1 + mu * eye(size(P{1},2),size(P{1},2));
        H_right = mu*(Z{V+1}*Z{V+1}'-Z{V+1}-Z{V+1}');
        H_all = tmp2 + Y4*(Z{V+1}-eye(N)) + mu * Er * (eye(N)-Z{V+1}');        
        H = lyap(H_left,H_right,H_all);

        for v=1:V 
            % update P  ok
            P{v} = (Xc{v}+Y1{v}/mu-E1{v})*H'*inv(H*H');

            % update Xc ok
            tmp1 = eye(N) + index{v} + (eye(N)-Z{v})*(eye(N)-Z{v}');
            tmp2 = E1{v}+P{v}*H + X{v}*index{v} + E2{v}*(eye(N)-Z{v}') - (Y1{v}+Y2{v}*index{v}+Y3{v}*(eye(N)-Z{v}'))/mu ;
            Xc{v} = tmp2/tmp1;
               
            % update Z  ok
            Q = EuDist2(U,U,0);
            tmp1 = Xc{v}'*Xc{v} + eye(N);
            tmp2 = (Xc{v}'*Y3{v}-Y{v}-0.5*param.lambda3*Q')/mu + Xc{v}'*Xc{v} + G{v} - Xc{v}'*E2{v};
            tmpZ = tmp1\tmp2;
            
            Z1 = zeros(N,N);
            for is = 1:N
                ind_c = 1:N;
                ind_c(is) = [];
                Z1(is,ind_c) = EProjSimplex_new(tmpZ(is,ind_c));
            end
            Z{v} = Z1;
            clear tmpZ Z1;

            % update E1 ok
            E1 = UpdateE1(Xc, P, H, Y1, mu, param.lambda1, V);
       
            % update E2
            E2 = UpdateE2(Xc, Z, Y3, mu, param.lambda2, V);

            % update Y1-Y3
            Y1{v} = Y1{v} + mu*(Xc{v} - P{v}*H - E1{v});
            Y2{v} = Y2{v} + mu*(Xc{v}*index{v} - X{v}*index{v});
            Y3{v} = Y3{v} + mu*(Xc{v} - Xc{v}*Z{v} - E2{v});

        end

        % update Er
        Er = UpdateEr(H,Z{V+1},Er,Y4,mu,param.lambda4);
        % update Y4
        Y4 = Y4 + mu*(H- H*Z{V+1}-Er);


        Z_tensor = cat(3, Z{:,:});
        Y_tensor = cat(3, Y{:,:});
        z = Z_tensor(:);
        y = Y_tensor(:);

        % update G
        sX = [N, N, V+1];

        [g, objV] = wshrinkObj(z + 1/mu*y,1/mu,sX,1,3);
        G_tensor = reshape(g,sX);

        y = y + mu*(z-g);
        
            % ------------------- U -------------------- % ok
        sum_Z = 0;
        for iv = 1:nv
            sum_Z = sum_Z + Z{iv};
        end
        sum_Z = (sum_Z+sum_Z')*0.5;
        LZv = diag(sum(sum_Z))-sum_Z;
        LZv = (LZv+LZv')*0.5;
        try
            opts.tol = 1e-4; 
            [U,~] = eigs(LZv,cls_num,'sa',opts);   % U: n*num_cluster
        catch ME
            if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
                opts.tol = 1e-4; 
                [U,~] = eigs(LZv, eye(size(LZv)),cls_num,'sa',opts);
            else
                rethrow(ME);
            end
        end 
% CONVERGENCE-------------------- 
        Isconverg = 1;
        for v=1:V 
            if (norm(Xc{v} - Xc{v}*Z{v} - E2{v}, inf)>param.epson)
                Isconverg = 0;
            end
            G{v} = G_tensor(:,:,v);
            Y_tensor = cat(3, YT{:,:});
            Y{v} = Y_tensor(:,:,v);  

            if (norm(Z{v}-G{v}, inf)>param.epson)
                Isconverg = 0;
            end
            if (norm(Xc{v}-P{v}*H - E1{v}, inf)>param.epson)
                Isconverg = 0;
            end
        end

        if (norm(H-H*Z{V+1}-Er)>param.epson)
            Isconverg = 0;
        end

        if (iter>100)
            Isconverg = 1;
        end
        iter = iter + 1;
        mu = min(mu*pho_mu, max_mu);
        
        Norm1 = [Norm1 norm(Xc{1}-P{1}*H - E1{1}, inf)];
        Norm2 = [Norm2 norm(Xc{1} - Xc{1}*Z{1} - E2{1}, inf)];
        Norm3 = [Norm3 norm(Z{1}-G{1}, inf)]; 
    end

    Fng =  NormalizeFea(U,1);
    for iter_c = 1 : 10 
        pre_labels = kmeans(real(Fng),cls_num,'maxiter',1000,'replicates',20,'EmptyAction','singleton');
        result_cluster = ClusteringMeasure(gt, pre_labels)*100;
        ACC(iter_c) = result_cluster(1); 
        NMI(iter_c) = result_cluster(2); 
        PUR(iter_c) = result_cluster(3);  
        result = compute_f(gt, pre_labels);
        Fscore(iter_c)=result(1)*100;
        Precision(iter_c)=result(2)*100;
        Recall(iter_c)=result(3)*100;
    end
    Result(1,:) = ACC;
    Result(2,:) = NMI;
    Result(3,:) = PUR;
    Result(4,:) = Fscore;
    Result(5,:) = Precision;
    Result(6,:) = Recall;
    Result(7,1) = mean(ACC);
    Result(7,2) = mean(NMI);
    Result(7,3) = mean(PUR);
    Result(7,4) = mean(Fscore);
    Result(7,5) = mean(Precision);
    Result(7,6) = mean(Recall);
    Result(8,1) = std(ACC);
    Result(8,2) = std(NMI);
    Result(8,3) = std(PUR);
    Result(8,4) = std(Fscore);
    Result(8,5) = std(Precision);
    Result(8,6) = std(Recall);
%     save([char(Dataname),'_result.mat'],'Result');
    filename = strcat(Dataname,'_' ,ratio, '.mat');
    save(filename,'Result')
    fprintf('ACC=%.4f,NMI=%.4f, Purity=%.4f,\n',mean(ACC),mean(NMI),mean(PUR));
end
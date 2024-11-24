% function [acc nmi] = demo_add_E1E3(X, gt, param)
function [acc nmi] = demo_add_E1E3(X, gt, param, Indicator)
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
    H = zeros(param.dim_h, N);
%     index = eye(N);
%     Indicator = BuildIndicator(V, N, param.ratio);
    for v=1:V
        index{v} = diag(Indicator(:,v));
        
        Z{v} = zeros(N,N);
        G{v} = zeros(N,N);
        Y4{v} = zeros(N,N);

        Xc{v} = 0.1*randn(size(X{v},1), N);

        E1{v} = zeros(size(X{v},1), N);
        E3{v} = zeros(size(X{v},1), N);

        Y1{v} = zeros(size(X{v},1), N);
        Y2{v} = zeros(size(X{v},1), N);
        Y3{v} = zeros(size(X{v},1), N);

        P{v} = 0.1*randn(size(X{v},1), param.dim_h);
    end

    Z_tensor = cat(3, Z{:,:});
    G_tensor = cat(3, G{:,:});
    Y_tensor = cat(3, Y4{:,:});
    w = zeros(N*N*V, 1);
    g = zeros(N*N*V, 1);

    for i=1:ModCount
        YT{i} = Y_tensor;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%% train
    acc = [];ACC = [];NMI=[];
    Norm1 = [];Norm2 = [];Norm3=[];
    tmp1 = 0;tmp2 = 0;
    while(Isconverg == 0)
        fprintf('----processing iter %d--------\n', iter+1);
        % update H  公式正确
        tmp1 = 0;tmp2 = 0;
        for v=1:V 
            tmp1 = tmp1 + P{v}'*P{v};
            tmp2 = tmp2 + P{v}'*(Xc{v}+Y1{v}/mu - E1{v}); 
        end
        H = inv(tmp1)*tmp2;

        for v=1:V 
            % update P  公式正确
            P{v} = (Xc{v}+Y1{v}/mu-E1{v})*H'*inv(H*H');

            % update Xc 公式正确
            tmp1 = eye(N) + index{v} + (eye(N)-Z{v})*(eye(N)-Z{v}');
            tmp2 = E1{v}+P{v}*H + X{v}*index{v} + E3{v}*(eye(N)-Z{v}') - (Y1{v}+Y2{v}*index{v}+Y3{v}*(eye(N)-Z{v}'))/mu ;
            Xc{v} = tmp2*inv(tmp1);

            % update Z 公式正确
            tmp1 = Xc{v}'*Xc{v} + eye(N);
            tmp2 = (Xc{v}'*Y3{v}-Y4{v})/mu + Xc{v}'*Xc{v} + G{v} - Xc{v}'*E3{v};
            Z{v} = inv(tmp1)*tmp2;

            % update E1
            E1 = UpdateE1(Xc, P, H, Y1, mu, param.lambda1, V);
            
            
            % update E3
            E3 = UpdateE3(Xc, Z, Y3, mu, param.lambda3, V);

            % update Y1-Y3
            Y1{v} = Y1{v} + mu*(Xc{v} - P{v}*H - E1{v});
            Y2{v} = Y2{v} + mu*(Xc{v}*index{v} - X{v}*index{v});
            Y3{v} = Y3{v} + mu*(Xc{v} - Xc{v}*Z{v} - E3{v});

        end


        Z_tensor = cat(3, Z{:,:});
        Y_tensor = cat(3, Y4{:,:});
        z = Z_tensor(:);
        y = Y_tensor(:);

        % update G

        [g, objV] = wshrinkObj(z + 1/mu*y,1/mu,sX,1,3);
        G_tensor = reshape(g,sX);

        y = y + mu*(z-g);

        Isconverg = 1;
        for v=1:V 
            if (norm(Xc{v} - Xc{v}*Z{v} - E3{v}, inf)>param.epson)
                Isconverg = 0;
    %             norm(Xc{v} - Xc{v}*Z{v} - E3{v}, inf)
            end
            % update G
            G{v} = G_tensor(:,:,v);
            Y_tensor = cat(3, YT{:,:});
            Y{v} = Y_tensor(:,:,v);  
            % udate Y4

            % Y4{v} = Y4{v} + mu*(Z{v} - G{v});

            if (norm(Z{v}-G{v}, inf)>param.epson)
                Isconverg = 0;
    %             norm(Z{v}-G{v}, inf)
            end
            if (norm(Xc{v}-P{v}*H - E1{v}, inf)>param.epson)
                Isconverg = 0;
    %             norm(Xc{v}-P{v}*H - E1{v}, inf)
            end
        end
        if (iter>100)
            Isconverg = 1;
        end
        iter = iter + 1;
        mu = min(mu*pho_mu, max_mu);
        
%         norm(Xc{1}-P{1}*H - E1{1}, inf)
%         norm(Xc{1}*index{1} - X{1}*index{1}, inf)
%         norm(Xc{1} - Xc{1}*Z{1} - E3{1}, inf)
%         norm(Z{1}-G{1}, inf)
        Norm1 = [Norm1 norm(Xc{1}-P{1}*H - E1{1}, inf)];
        Norm2 = [Norm2 norm(Xc{1} - Xc{1}*Z{1} - E3{1}, inf)];
        Norm3 = [Norm3 norm(Z{1}-G{1}, inf)];
        
        S = 0;
        for v=1:V 
            S = S + abs(Z{v}) + abs(Z{v}');
        end
        C = SpectralClustering(S,cls_num);
        ACC = [ACC Accuracy(C,gt)];
        [A nmi avgent] = compute_nmi(gt,C);
        NMI = [NMI nmi];
%         groups = kmeans(H',cls_num,'maxiter',1000,'replicates',20,'EmptyAction','singleton');
%         acc = [acc Accuracy(groups,gt)];
        
%         if Accuracy(C,gt)>0.95
%             break;
%         end
    end
    save('converge.mat','ACC','NMI','Norm1','Norm2','Norm3')
    % x_axis = 1:size(ACC,2);
    % plot(x_axis, ACC,'-*r', x_axis, NMI,'-ob');
    % legend('Z','H');
    S = 0;
    for v=1:V 
        S = S + abs(Z{v}) + abs(Z{v}');
    end
    C = SpectralClustering(S,cls_num);
    acc = Accuracy(C,gt);
    [A nmi avgent] = compute_nmi(gt,C);
%     [acc nmi]
    % groups = kmeans(H',cls_num,'maxiter',1000,'replicates',20,'EmptyAction','singleton');
    %         acc = [acc Accuracy(groups,gt)];
end
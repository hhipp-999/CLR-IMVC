function [E1] = UpdateE1(Xc, P, H, Y1, mu, lambda1, V)
        F = [];
        for k=1:V
             F = [F;Xc{k}-P{k}*H+Y1{k}/mu];
        end
%         F = [X{1}-X{1}*Z{1}+Y{1}/mu;X{2}-X{2}*Z{2}+Y{2}/mu;X{3}-X{3}*Z{3}+Y{3}/mu];
        [Econcat] = solve_l1l2(F,lambda1/mu);
        
        beg_ind = 0;
        end_ind = 0;
        for k=1:V
            if(k>1)
                beg_ind = beg_ind+size(Xc{k-1},1);
            else
                beg_ind = 1;
            end
            end_ind = end_ind+size(Xc{k},1);
            E1{k} = Econcat(beg_ind:end_ind,:);
        end
end
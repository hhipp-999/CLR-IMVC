% function [indicator] = BuildIndicator(V, N, ratio)
%     ratio = 1-ratio;
%     indicator = zeros(1, N*V);
%     new_index = BinaryRandom(indicator, ratio);
%     while(check(indicator, new_index, N, V) ~= 1)
%         new_index = BinaryRandom(indicator, ratio);
%     end
%     indicator(new_index) = 0;
%     indicator = reshape(indicator, [N, V]);
% %     list_class = unique(gt);
% %     indicator = ones(1,N);
% %     if ratio ~= 0
% %         for i=1:size(list_class, 1)
% %             classes = list_class(i);
% %             index = find(gt==classes);
% % %             index = reshape(index, size(index, 1));
% %             indicator_index = BinaryRandom(index, ratio);
% % %             indicator_index = reshape(indicator_index, size(indicator_index, 1));
% %             indicator(index(indicator_index)) = 0;
% %         end
% %     end
% %     indicator = diag(indicator);
% end
% function [flag] = check(indicator, index, N, V)
%     flag = 1;
%     indicator(index)=1;
%     if (any(sum(reshape(indicator, [N, V]), 2) == 0) == 1)
%         flag = 0;
%     end
% end
% function [array] = BinaryRandom(index, ratio)
%     new_index = randperm(size(index, 2)); % œ¬±Í÷ÿ≈≈–Ú
%     array = new_index(1:int32(ratio*size(index, 2)));
% end
function [indicator] = BuildIndicator(V, N, ratio)
    if (ratio < 0 || ratio > (V-1)/V)
        fprintf('error: missingg rate out of the range');
    end
    indicator = zeros(N, V);
    matrix_retain = randi(V, N, 1);
    matrix_other = remain(matrix_retain, V, N, ratio);
    indicator = merge(indicator, matrix_retain, matrix_other, N, V, ratio);
end
function [matrix] = remain(matrix_retain, V, N, ratio)
    matrix = zeros(N, (V-1));
    for i=1:N
       view = 1;
       for j=1:V-1
           if (matrix_retain(i,1) == view)
               view = view + 1;
           end
           matrix(i, j) = view;
           view = view + 1;
       end
    end
end
function [indicator] = merge(indicator, matrix_retain, matrix_other, N, V, ratio)
    rand_index = randperm(N*(V-1));
    rand_index = rand_index(1:int32(ratio*N*V));
    binary_index = ones(N*(V-1),1);
    binary_index(rand_index) = 0;
    binary_index = reshape(binary_index,[N, V-1]);
    matrix_other = matrix_other.*binary_index;
    for i=1:N
        miss_index = find(matrix_other(i,:)~=0);
        indicator(i, matrix_retain(i, 1)) = 1;
        for j=1:size(miss_index, 2)
            indicator(i, matrix_other(i,miss_index(j))) = 1;
        end
    end
end

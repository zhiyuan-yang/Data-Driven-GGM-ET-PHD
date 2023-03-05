function [Zkp] = split(Zkp,model)
lambda = model.lambda;
num_p = length(Zkp);

for i=1:1:num_p
    num_w = Zkp(i).N_W; %number of cell in current partition
    curr_p = Zkp(i);    %current partition
    split_idx = [];     %idx of cell to be splitted
    for j = 1:1:num_w
        num_meas = length(curr_p.P(j).idx);
        [~ , max_n] = max(poisspdf(num_meas,lambda*[1:5]));
        if max_n > 1
            temp_meas = Zkp(i).P(j).W;
            temp_idx = Zkp(i).P(j).idx;
            idx = kmeans(temp_meas',max_n);
            new_w = struct('W',temp_meas(:,find(idx == 1)),'idx',temp_idx(find(idx == 1)));
            for k = 2:1:max_n
                new_w(k).W = temp_meas(:,find(idx == k));
                new_w(k).idx = temp_idx(find(idx == k));
            end
            new_w = new_w + max_n - 1;
        else
            Zkp(i).P(new_w).W = Zkp(i).P(j).W;
            Zkp(i).P(new_w).idx = Zkp(i).P(j).idx;
            new_w = new_w + 1;
        end
    end
    Zkp(i).N_W = new_w;
end


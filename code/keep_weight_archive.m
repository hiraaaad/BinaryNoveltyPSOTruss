function [weight_archive,total_eval,best_feas,GenerationData] = keep_weight_archive(pos,weight_archive,best_feas,benchmark,GenerationData,gen)
ps =  size(pos,1);
dim = numel(pos(1,:));
for i = 1: size(pos,1)
    kp = key_p(pos(i,:));
    weight_archive.pos{end+1} = kp;
    if benchmark ~= 68
        idx = bin2dec(kp)+1;
    else
        N=68;
        strbin = kp; %char(randi(2,1,N)+'0'-1);
        pows2 = 2.^(N-1:-1:0);
        idx=pows2*(strbin-'0')';
    end
    M_upper = pos(i,:);
    [best_feas_design,final_eval,total_eval] = run_subsolver(benchmark, M_upper);
    weight_archive.weight(end+1) = best_feas_design(end);
    weight_archive.final_eval(end+1) = final_eval;
    best_feas(kp)=best_feas_design;
    GenerationData(end+1,:)=[idx,gen,best_feas_design(end)];
end
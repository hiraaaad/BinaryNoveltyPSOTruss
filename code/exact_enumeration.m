function weight_archive = exact_enumeration(benchmark)

switch benchmark
    case 10
        Dimension=10;
        load('testproblem_data/10bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
    case 25
     Dimension=8;
     load('testproblem_data/25bar_info.mat','GSCP','basic_node','N_member','N_node','X_const_val','sym_M','DOF_constrained')
     is_constrained_DOF=zeros(1,2*N_node);
     is_constrained_DOF(DOF_constrained)=1;
    case 52
     Dimension = 12;
     load('testproblem_data/52bar_info.mat','GSCP','basic_node','N_member','N_node','X_const_val','sym_M','DOF_constrained')
     is_constrained_DOF=zeros(1,2*N_node);
     is_constrained_DOF(DOF_constrained)=1;
     
end

total = 2^Dimension;
pos = dec2bin(0:total-1)' - '0';
pos = pos';

weight_archive = ones(total,1).*1e150;

for i = 2: total
    M_upper = pos(i,:);
    [best_feas_design] = run_subsolver(benchmark, M_upper);
    weight_archive(i,1) = best_feas_design(end);
end
     
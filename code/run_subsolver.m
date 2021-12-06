function [best_feas_design,final_eval,total_eval] = run_subsolver(benchmark, M_upper)

switch benchmark
    case 10
        [best_feas_design,final_eval,total_eval] = bar10(M_upper);
    case 15
        [best_feas_design,final_eval,total_eval] = bar15(M_upper);
    case 25
        [best_feas_design,final_eval,total_eval] = bar25(M_upper);
    case 47
        [best_feas_design,final_eval,total_eval] = bar47(M_upper);
    case 68
        [best_feas_design,final_eval,total_eval] = bar68(M_upper);
    case 52
        [best_feas_design,final_eval,total_eval] = bar52(M_upper);
    case 72
        [best_feas_design,final_eval,total_eval] = bar72(M_upper);
    case 200
        [best_feas_design,final_eval,total_eval] = bar200(M_upper);
    case 224
        [best_feas_design,final_eval,total_eval] = bar224(M_upper);  
end


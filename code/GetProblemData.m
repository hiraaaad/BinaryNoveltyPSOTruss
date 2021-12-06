function [opt,D,N_node,N_member,GSCP,sym_M,sym_A,X_indep_ind,basic_node,basic_member,X_const_ind,X_const_val,DOF_constrained,Fext,M_elasticity,Density,dis_all,SigT_all, SigC_all,kappa,Fy,section_no,Min_A,Max_A,D_X, U_X, minL,maxL]=GetProblemData(opt)

if strcmp(opt.ProblemName,'10barDiscrete')
    load('testproblem_data/10bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 10; % 10 bar truss discrete
elseif strcmp(opt.ProblemName,'15barDiscrete')
    load('testproblem_data/15bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 15; % 15 bar truss discrete
elseif strcmp(opt.ProblemName,'25barDiscrete')
    load('testproblem_data/25bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 25; % 25 bar truss discrete
elseif strcmp(opt.ProblemName,'25barDiscrete2')
    load('testproblem_data/26bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 26; % 25 bar truss discrete 2
elseif strcmp(opt.ProblemName,'52barDiscrete')
    load('testproblem_data/52bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 52; % 52 bar truss discrete
elseif strcmp(opt.ProblemName,'72barDiscrete')
    load('testproblem_data/72bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 72; % 52 bar truss discrete
elseif strcmp(opt.ProblemName,'72barDiscrete2')
    load('testproblem_data/72bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 73; % 52 bar truss discrete
elseif strcmp(opt.ProblemName,'200barDiscrete')
    load('testproblem_data/200bar_info.mat')
    Fy = 1e12;
    minL=0;
    maxL=1e10;
    section_no = 2; % 200 bar truss discrete : sections 2
end

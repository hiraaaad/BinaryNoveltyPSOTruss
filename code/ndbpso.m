function T = ndbpso(benchmark)
% this script is for novelty-driven binary PSO
% function output = ndbspo(id,benchmark,Particle_Number,Max_Gen,radius_novelty)
infeasible_upper = containers.Map; %visited infeasible topologies
feasible_upper = containers.Map; %visited feasible topologies
best_feas = containers.Map;
archive_novelty = [];
%
id = 1; % this id stems from the .sh file to perform multiple experiments with different rng.
Particle_Number=30;
Max_Gen=1000;
radius_novelty=3;
% random number generator
RandSeedNo=randi([1,100000],1,100000);
RandSeedNo = RandSeedNo(id); % random seed number
rng(RandSeedNo)
% 
% testproblems
switch benchmark
    case 10
        Dimension=10;
        load('testproblem_data/10bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
    case 15
     Dimension=15;
     load('testproblem_data/15bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
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
    case 72
        Dimension = 16;
     load('testproblem_data/72bar_info.mat','GSCP','basic_node','N_member','N_node','X_const_val','sym_M','DOF_constrained')
     is_constrained_DOF=zeros(1,3*N_node);
     is_constrained_DOF(DOF_constrained)=1;
    case 68
        Dimension = 68;
        load('testproblem_data/68bar_info.mat','GSCP','basic_node','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        is_constrained_DOF=zeros(1,2*N_node);
        is_constrained_DOF(DOF_constrained)=1;
    case 47
        Dimension = 27;
        load('testproblem_data/47bar_info.mat','GSCP','basic_node','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        is_constrained_DOF=zeros(1,2*N_node);
        is_constrained_DOF(DOF_constrained)=1;
    case 224
        Dimension = 32;
        load('testproblem_data/224bar_info.mat','GSCP','basic_node','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        is_constrained_DOF=zeros(1,3*N_node);
        is_constrained_DOF(DOF_constrained)=1;
    case 200
     Dimension=29;
     load('testproblem_data/200bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
end
% information storage
weight_archive = struct();
weight_archive.pos = [];
weight_archive.weight = [];
weight_archive.final_eval = [];
TotalEval = 0;
GenerationData = [];
% binary PSO
% keep an archive
VRmin=range(1);
VRmax=range(2);

me=Max_Gen;
ps=Particle_Number;
D=Dimension;
cc=[1 1];   %acceleration constants
iwt=0.9-(1:me).*(0.5./me);

phi_max = 5;
phi_min = 1;
phi = phi_max - [1:1:me].*((phi_max-phi_min)./me);


%
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end

mv=0.5*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin = -6;
Vmax = -Vmin;
%
pos = generate_initial_population(Dimension,ps,is_constrained_DOF,GSCP, basic_node, DOF_constrained, N_node,sym_M,benchmark); % Initial population

for k=1:ps
        if size(archive_novelty,1)>0
            if ~ismember(pos(k,:),archive_novelty,'rows')
                if rand <= 1
                    archive_novelty(end+1,:) = pos(k,:);
                end
            end
        else
            archive_novelty(end+1,:) = pos(k,:);
        end
end

fitcount=0;

result = evaluate_novelty(pos,archive_novelty,radius_novelty);

gen = 0;

[weight_archive,total_eval,best_feas,GenerationData] = keep_weight_archive(pos,weight_archive,best_feas,benchmark,GenerationData,gen);

for k=1:ps
    kp = key_p(pos(k,:));
    feasible_upper(kp)=1;     
end

fitcount=fitcount+ps;
TotalEval = 0;

vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=result; %initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1);
g_res(1)=gbestval;

i=1;

while i<=me
    
    for k=1:ps
        
        gen = i;
        % update velocity and position - binary variables
        v = vel(k,:);
        p = pos(k,:);
        % velocity boundary
        v(v>Vmax) = Vmax;
        v(v<Vmin) = Vmin;
                
        transform_vel = sigmoid(v,phi(i)); 
        r = rand(1,D);
        
        p(r<transform_vel) = 1;
        p(r>=transform_vel) = 0;
        
        offspring_pass = false;
        
        % decide for new offspring
        f_x = fitness_topology(p,is_constrained_DOF,GSCP, basic_node, DOF_constrained, N_node,sym_M,benchmark);
        X = p;
        ctr = 0;
        while offspring_pass == false
            ctr = ctr + 1;
           kp = key_p(X);
            
            if sum(f_x == 0) == 3 % feasible
                % visited before : no calculation
                if ~isKey(feasible_upper,kp)
                      [weight_archive,total_eval,best_feas,GenerationData] = keep_weight_archive(X,weight_archive,best_feas,benchmark,GenerationData,gen); 
                    fitcount  = fitcount+ 1;
                    TotalEval = 0;
                    feasible_upper(kp)=1; 
                    % add it to feasible 
                    % else to infeasible
                % new : truss weight calculation
                % else: nothing new 
                end
                offspring_pass = true;
            else % infeasible
                if ~isKey(infeasible_upper,kp) % not visited before, add it to rejection list
                    infeasible_upper(kp)=1;
                    
                    if benchmark ~= 68
                       idx = bin2dec(kp)+1;
                    else
                        N=68;
                        strbin = kp; %char(randi(2,1,N)+'0'-1);
                        pows2 = 2.^(N-1:-1:0);
                        idx=pows2*(strbin-'0')';
                    end
                    
                    GenerationData(end+1,:)=[idx,gen,1e12];
                    if ~isempty(archive_novelty)
                            if ~ismember(X,archive_novelty,'rows')
                                archive_novelty(end+1,:) = X;
                            end
                    end
                end
% repair
[Y,ky] = ea_1p1_bitflip(X);
f_y = fitness_topology(Y,is_constrained_DOF,GSCP, basic_node, DOF_constrained, N_node, sym_M, benchmark);
if sum(f_y == 0) ~= 3 %infeasible
    if ~isKey(infeasible_upper,ky)
        infeasible_upper(ky)=1;
        if ~isempty(archive_novelty)
            if ~ismember(Y,archive_novelty,'rows')
                archive_novelty(end+1,:) = Y;
                
                if benchmark ~= 68
                    idx = bin2dec(ky)+1;
                else
                    N=68;
                    strbin = ky; %char(randi(2,1,N)+'0'-1);
                    pows2 = 2.^(N-1:-1:0);
                    idx=pows2*(strbin-'0')';
                end
                GenerationData(end+1,:)=[idx,gen,0];
            end
        end
    end
else
    if ~isKey(feasible_upper,ky)
        feasible_upper(ky)=1;
%                         [weight_archive,GenerationData] = keep_weight_archive_raw(Y,Fitness,weight_archive,id,GenerationData,gen);
        [weight_archive,total_eval,best_feas,GenerationData] = keep_weight_archive(Y,weight_archive,best_feas,benchmark,GenerationData,gen);
        fitcount  = fitcount+ 1;
        TotalEval = 0;
    end
end
        if lex_compare(f_y,f_x)
            X = Y;
            f_x = f_y;
        end
    end
end
p = X;

% put them back
vel(k,:) = v;
pos(k,:) = p;

end
    
    result=evaluate_novelty(pos,archive_novelty,radius_novelty);
        
    for k=1:ps
        
        recalculate_pbest = evaluate_novelty(pos(k,:),archive_novelty,radius_novelty,pbest(k,:)); % recalculate novelty score of pbest
        
           tmp=(recalculate_pbest<result(k));
           temp=repmat(tmp,1,D);
           pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
           pbestval(k)=tmp.*recalculate_pbest+(1-tmp).*result(k); %update the pbest
    end
        
    [pbestval_sorted, idx_sorted ] = sort(pbestval);    
    tied = sum(pbestval_sorted(2:end) == pbestval_sorted(1));
    
    if tied == 0
%       there is only one top best pbestval
        idx = idx_sorted(1);
            gbest=pbest(idx,:);
            gbestval=pbestval(idx);
            gbestrep=repmat(gbest,ps,1);%update the gbest
%            end
    else
        % there is a tie in top best pbestval
        top_candidates_quota = tied + 1; % how many tied top best pval are avaliable including first index
        top_candidates = pbest(idx_sorted(1:top_candidates_quota),:);
                
            if sum(sum(gbest == top_candidates,2) == numel(gbest)) == 0
            %if it was true then gbest remains the same
                % choose randomly 
                perm = randi([1,top_candidates_quota]);
                idx = idx_sorted(perm);
                
                 if pbestval(k)<gbestval
                    gbest=pbest(k,:);
                    gbestval=pbestval(k);
                    gbestrep=repmat(gbest,ps,1);%update the gbest
                 end
            end
            
    end   

% adding to archive
    for k=1:ps
        if size(archive_novelty,1)>0
            if ~ismember(pos(k,:),archive_novelty,'rows')
                if rand <= 1
                    archive_novelty(end+1,:) = pos(k,:);
                end
            end
        else
            archive_novelty(end+1,:) = pos(k,:);
        end
    end

    for k = 1: ps
    
        aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbestrep(k,:)-pos(k,:)); % second and third term of position update
        vel(k,:)=iwt(i).*vel(k,:)+aa(k,:); % velocity update
    
    end
    i=i+1;
    if fitcount>=1e6 % Max_FES controller for computational expense of the problems.
        break;
    end
end

position=gbest;
value=gbestval;
iteration=i;
num_FES=fitcount;
T=struct2table(weight_archive);
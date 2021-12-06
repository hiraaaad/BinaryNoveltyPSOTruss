function T = pso_novelty_truss_optimizer(id,benchmark)
% parameters
Max_FES=1000;
range=[0,1];
%
archive_novelty = [];
radius_novelty = 5; %k nearest neighbors

infeasible_upper = containers.Map; %visited infeasible topologies
feasible_upper = containers.Map; %visited feasible topologies
best_feas = containers.Map;
% start a pop from a predefined list
switch benchmark
    case 10
        Dimension=10;
%         load('valid_topologies_10bar.mat')
        load('10bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        Max_Gen=10;
        Particle_Number=10;
        str_save = 'run_10bar_novelty_iter_';
        DD = 2;
    case 11
        Dimension=10;
%         load('valid_topologies_10bar.mat')
        load('10bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        Max_Gen=10;
        Particle_Number=10;
        str_save = 'run_11bar_novelty_iter_';
        DD = 2;
    case 15
        Dimension=15;
%         load('valid_topologies_15bar.mat')
        load('15bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        Max_Gen=10;
        Particle_Number=10;
        str_save = 'run_15bar_novelty_iter_';
        DD = 2;
    case 52
        Dimension=12;
%         load('valid_topologies_52bar.mat')
        load('52bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        Max_Gen=100;
        Particle_Number=30;
        str_save = 'run_52bar_novelty_iter_';
        DD = 2;
    case 72
        Dimension=16;
%         load('valid_topologies_72bar.mat')
        load('72bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        Max_Gen=100;
        Particle_Number=30;
        str_save = 'run_72bar_novelty_iter_';
        DD = 3;
    case 73
        Dimension=16;
%         load('valid_topologies_72bar.mat')
        load('72bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        Max_Gen=100;
        Particle_Number=30;
        str_save = 'run_73bar_novelty_iter_';
        DD = 3;
    case 200
        Dimension=29;
%         load('valid_topologies_200bar.mat')
        load('200bar_info.mat','GSCP','basic_node','is_constrained_DOF','N_member','N_node','X_const_val','sym_M','DOF_constrained')
        Max_Gen=100;
        Particle_Number=30;
        str_save = 'run_200bar_novelty_iter_';
        DD = 2;
end

% load('run_10bar_fitness.mat','F');

RandSeedNo=randi([1,100000],1,100000);
RandSeedNo = RandSeedNo(id); % random seed number
rng(RandSeedNo)

weight_archive = struct();
weight_archive.pos = [];
weight_archive.weight = [];
weight_archive.final_eval = [];
TotalEval = 0;
% repair an individual to be fit in a predefined list
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
pos = generate_feasible_initial_population(benchmark,ps) % Initial population

fitcount=0;

result = evaluate_novelty(pos,archive_novelty,radius_novelty);

[weight_archive,total_eval,best_feas] = keep_weight_archive(pos,weight_archive,best_feas,benchmark);

for k=1:ps
    kp = key_p(pos(k,:));
    feasible_upper(kp)=1;     
end

fitcount=fitcount+ps;
TotalEval = TotalEval+ total_eval;

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
        f_x = fitness_topology(p,is_constrained_DOF,GSCP, basic_node, DOF_constrained, DD, N_node,sym_M,benchmark);
        X = p;
        ctr = 0;
        while offspring_pass == false
            ctr = ctr + 1;
           kp = key_p(X);
            
            if sum(f_x == 0) == 3 % feasible
                % visited before : no calculation
                if ~isKey(feasible_upper,kp)
                    [weight_archive,total_eval,best_feas] = keep_weight_archive(X,weight_archive,best_feas,benchmark)
                    fitcount  = fitcount+ 1;
                    TotalEval = TotalEval+ total_eval;
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
                    % else : it has been seen before and already in
                    % rejection listx
                end
                % resample
%                 r = rand(1,D);
%                 p(r<transform_vel) = 1;
%                 p(r>=transform_vel) = 0;
%                 f_x = fitness_topology(p,is_constrained_DOF,GSCP, basic_node, DOF_constrained, 2, N_node,sym_M,benchmark);
                % repair
                [Y,ky] = ea_1p1_bitflip(X);
                f_y = fitness_topology(Y,is_constrained_DOF,GSCP, basic_node, DOF_constrained, DD, N_node, sym_M, benchmark);
                if sum(f_y == 0) ~= 3 %infeasible
                    if ~isKey(infeasible_upper,ky)
                        infeasible_upper(ky)=1;
                    end
                else
                    if ~isKey(feasible_upper,ky)
                        feasible_upper(ky)=1;
                        [weight_archive,total_eval,best_feas] = keep_weight_archive(Y,weight_archive,best_feas,benchmark)
                        fitcount  = fitcount+ 1;
                        TotalEval = TotalEval+ total_eval;
                    end
                end
                if lex_compare(f_y,f_x)
                    X = Y;
                    f_x = f_y;
                end
            end
        end
            
        
        
%         % repair p
%         
%         
%         
%         
%         counter = 0;
%         X = p;
%         while sum(f_x == 0) ~= 3
%             r = rand(1,D);
%             p(r<transform_vel) = 1;
%             p(r>=transform_vel) = 0;
%             X = p;
%             f_x = fitness_topology(X,is_constrained_DOF,GSCP, basic_node, DOF_constrained, 2, N_node,sym_M,benchmark)
%             counter = counter+1
%         end

%             Y = ea_rls_bitflip(X);
%             f_y = fitness_topology(Y,is_constrained_DOF,GSCP, basic_node, DOF_constrained, 2, N_node, sym_M, benchmark);
%             if lex_compare(f_y,f_x)
%                 X = Y;
%                 f_x = f_y;
%             end
%             if counter > 1000
%                 counter = 0;
%                 index = randperm(size(valid_topologies,1),Particle_Number);
%                 X = valid_topologies(index(1),:);
%                 f_x = fitness_topology(X,is_constrained_DOF,GSCP, basic_node, DOF_constrained, 2, N_node,sym_M,benchmark);
%             end
%         end
        p = X;
        
        % put them back
        vel(k,:) = v;
        pos(k,:) = p;
     
    end
%         if (sum(pos(k,:)>VRmax(k,:))+sum(pos(k,:)<VRmin(k,:)))==0
        
        result=evaluate_novelty(pos,archive_novelty,radius_novelty);
        
%         [weight_archive,total_eval,best_feas] = keep_weight_archive(pos,weight_archive,best_feas,valid_topologies,benchmark);

        
%         fitcount=fitcount+ ps;
%         TotalEval = TotalEval+ total_eval;

%         result_data(fitcount)=result(k);
%         if fitcount>=Max_FES
%             break;
%         end
        
    for k=1:ps
        
        recalculate_pbest = evaluate_novelty(pos,archive_novelty,radius_novelty,pbest(k,:)); % recalculate novelty score of pbest
        
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
%         if pbestval(idx)<gbestval
            gbest=pbest(idx,:);
            gbestval=pbestval(idx);
            gbestrep=repmat(gbest,ps,1);%update the gbest
%            end
    else
        % there is a tie in top best pbestval
        top_candidates_quota = tied + 1; % how many tied top best pval are avaliable including first index
        top_candidates = pbest(idx_sorted(1:top_candidates_quota),:);
        
%          if pbestval(idx)<gbestval
        
            if sum(sum(gbest == top_candidates,2) == numel(gbest)) == 0
            %if it was true then gbest remains the same
                % choose randomly 
                perm = randi([1,top_candidates_quota]);
                idx = idx_sorted(perm);
                
                 gbestval
                 if pbestval(k)<gbestval
                    gbest=pbest(k,:)
                    gbestval=pbestval(k);
                    gbestrep=repmat(gbest,ps,1);%update the gbest
                 end
            end
            
    end   
    
    

%        
%     end
%     end

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
   
    
%     tmp1=abs(repmat(gbest,ps,1)-pos)+abs(pbest-pos);
%     temp1=ones(ps,1);
%     temp2=ones(ps,1);

    for k = 1: ps
    
        aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbestrep(k,:)-pos(k,:)); % second and third term of position update
        vel(k,:)=iwt(i).*vel(k,:)+aa(k,:); % velocity update
    
    end
    
    i=i+1;
    
    if fitcount>=Max_FES
        break;
    end
    
    
%     if (i==me)&&(fitcount<Max_FES)
%         i=i-1;
%     end
    
%     X(end+1)=1./gbestval;
%     XX(end+1)= sum(gbest);
end

position=gbest;
value=gbestval;
iteration=i;
num_FES=fitcount;
% found_topology = find(sum(valid_topologies_10bar == gbest,2)==10);
T=struct2table(weight_archive);
fitcount

save([str_save,num2str(id),'.mat'])

% end

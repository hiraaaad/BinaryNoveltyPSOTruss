function [best_feas_design,final_eval,total_eval] = bar200(M_upper)
kk = 1;
load('testproblem_data/200bar_info.mat')
%Define Mechanical Properties
M_elasticity=30000;  %modulus of elasticity
SigT_all=10; % tensile strength
SigC_all=10; % comprehessive strength
Density=0.283; % density of truss material
kappa=3.96; % The parameter in the elastic buckling equation ( sigma_c<=kappa*EI/l^2)
dis_all=1e12; % allowable displacement along each DOF. If there is no displacement, set dis_all=1e12 
% Now specify the list number for the given cross section. Store the cross sections in m-file sections.m. Use a number greater than 100 when AISC-ASD spcifications are employed.
section_no=2; % The avilable section list for this problem. (see sections.m)
Min_A=0.1;
Max_A=3.2;
% *************************************************************************
opt.initial_M=1; % initial values for recombinant design for topology variables
opt.topology_efficiency_check=0; % check the necessary condition for efficiency of the sampled topology (deactivated)
opt.max_infeas=.5; % the threshold for ratio of constraint violating members. 
opt.WriteInterval=5; % Interval for writing the optimization process
opt.lambda_coeff=3; % coefficient for the population size
opt.mulambdaratio=.3; % ratio of parents to offspring
opt.maxiter_coeff=10; % coefficient for the maximum number of iterations
opt.Tccoeff=1; % Coefficient for the learning rate of the scaling factors
opt.Tscoeff=1; % Coefficient for the learning rate of the global step size
opt.ucr_red=.05; % reduction rate of the critical displacement constraint
opt.feastol=1e-7; % tolerance for constraint violation
% ********** No more modifications are required  from this line ***********
SF_buck = 1e12; %h

RandSeedNo = randi([1,1e9]);

rng(RandSeedNo)
if section_no>0
    availsec=sections(section_no);
    Min_A=min(availsec(:,1));
    Max_A=max(availsec(:,1));
end
N_loadcase=size(Fext,2);
is_constrained_DOF=zeros(1,D*N_node);
is_constrained_DOF(DOF_constrained)=1;
DOFonNAP=sum(reshape(is_constrained_DOF,D,N_node)); % used for checking necessary condition of stability

if opt.topology_efficiency_check==1 
     F0=max(abs(Fext)>0,[],2);
     DOFonNAP=DOFonNAP+max(reshape(F0,D,N_node))-1;
end

%Initial valuess of design variables
Mmean=opt.initial_M*ones(1,N_member); % topology
Amean=(Max_A+Min_A)/2*ones(1,N_member); % cross-section areas
% Find independent variables. Consider only variables that are independent
Xvar_ind=setdiff(X_indep_ind, X_const_ind);   % no shape analysis
% Avar_ind=setdiff(1:N_member,reshape(sym_A(2:end,:),1,numel(sym_A(2:end,:))));
Avar_ind = sym_A(1,:);
Mvar_ind = sym_M(1,:);
N_Avar=numel(Avar_ind); % Number of independent size variables
N_Mvar=numel(Mvar_ind); 
YMEAN=([Amean(Avar_ind) Mmean(Mvar_ind)]);
U_Y=[Max_A*ones(1,N_Avar) ones(1,N_Mvar)]; 
D_Y=[Min_A*ones(1,N_Avar) zeros(1,N_Mvar)];

STR=(U_Y-D_Y); % This is vector of scaling factors
Max_Kcond=1e12; % The limit that specifies whether the stiffness matrix is singular

SMEAN=.25; % global step size

% Calculating effective number of design parameters
N_Xvar = 0; % no shape optimization

N_var_all=( sqrt(N_Mvar)+sqrt(N_Xvar)+sqrt(N_Avar))^2;
N_DefC=(N_node*D-numel(DOF_constrained))*(dis_all<1e10);
N_constraint_all=(sqrt(N_DefC)+sqrt(N_member))^2*sqrt(N_loadcase);
N_var_eff=(N_var_all)* sqrt(1+N_constraint_all/N_var_all); %efective number of variables
% Setting lambda, mu, maxiter
lambda=round(opt.lambda_coeff*sqrt(N_var_eff)); 

maxiter=round(opt.maxiter_coeff*sqrt(N_var_eff));

mu=max(round((lambda*opt.mulambdaratio)),1); 
Tc=(1+opt.Tccoeff*N_var_eff*(N_var_eff+1)/mu)^(-1);   % learning rate of STR
Ts=1/sqrt(2*opt.Tscoeff*N_var_eff);   % learning rate of the global step size

weights=log(1+mu)-log(1:mu);
weights=weights/sum(weights); % weights for the recombination of the selected parents
% Initial values 
iter=0; %iteration number
bestf=1e180; % best penalized function value
Pcoeff=ones(1,N_member); % penalty coefficients for constriant violation of members 

best_feas_design=zeros(1,2*numel(Mmean)+D*N_node+1); % best feasible solution and weight so far
best_feas_design(end)=1e150;
Max_Max_Avar=(Max_A/Min_A)^1; % maximum value of the move limit ratio

if section_no>0
    Min_Max_Avar=max(availsec(2:end,1)./availsec(1:end-1,1)); % minimum value of the move limit ratio
else
    Min_Max_Avar=1.01;
end
Max_Avar=Max_Max_Avar; % move limit ratio
minPcoeff=1;  % minimum value of penalty coefficients
infeas_ratio=opt.max_infeas*ones(1,N_member); % maximum fraction of members that may violate a constraint
counteval=0; % evaluation count

stopping = false;
stopping_record = [];
stopping_evals = [];

while iter<maxiter %&& stopping == false % main optimization loop starts here
    iter=iter+1; 
    % presetting values for speed boost
    keep_M=zeros(2*lambda,N_member);
    Y=zeros(2*lambda,numel(STR)); 
    Y_W=Y;
    f=1e180*ones(1,2*lambda);
    constraint_vio=zeros(2*lambda,N_member);
    Z=zeros(2*lambda,numel(STR));
    S=SMEAN*ones(1,2*lambda); 
    for k=1:lambda %start sampling solutions
        stable=0;
        while ~stable %hirad: checks external stability c% Stay in this loop untill a solution that satisfies necessary conditions for stability (and efficiency, if applied) is generated
            S(k)=SMEAN*exp(Ts*randn); % mutate the global step size
            % now use the truncated normal distribution to sample a new design between D_Y and U_Y
            mincdf=normcdf(D_Y,YMEAN,S(k)*STR);
            maxcdf=normcdf(U_Y,YMEAN,S(k)*STR);
            r=max(1e-8,min(1-1e-8,rand(1,numel(U_Y)))); % we use 1e-8 bound to avoid infinity error
            Y(k,:)=norminv(mincdf+(maxcdf-mincdf).*r,YMEAN,S(k)*STR); %  the k-th candidate design with continuous values
            A=Amean;M=Mmean;
            M(Mvar_ind)=(Y(k,N_Xvar+N_Avar+1:end)>rand(1,N_Mvar)); % the topology is defined in M
            M(basic_member)=1; % enforce presence of basic members, if any
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1); % deteremine absence/presence of topologically coupled members            
            [NAP,NAPconnector]=find_NAP(GSCP,M); % find which nodes are active, and the number of the members connected to them
            NAPconnector_all=NAPconnector+DOFonNAP; % The number of members connected to each node + the number of reactions on the nodes
            
            chk1=(min(NAP(basic_node)==1)>.5); % necessary condition: Are all basic nodes active? 
            chk2=(sum(M) >= (sum(NAP)*D-numel(DOF_constrained))); % Necessary condition: Does the candidate truss have the minimum required number of members?
            chk3= min( max([ NAPconnector_all>=D ; NAP==0])      ); % Necessary condition: Are sufficient number of members connected to each node? 
            
            if (chk1 && chk2 && chk3)
                stable=1; % The candidate design is accepted if it satisfies all the necessary conditions for stability
            end
        end % The k-th solution that satisfies the necessary conditions for stability is generated
%         X(Xvar_ind)=Y(k,1:N_Xvar); % shape of the design
%         X=enforce_shape_sym(X,N_member,opt);% Enforce shape symmetry, if any

         M = zeros(1,N_member); %upper
         M(sym_M(1,:)) = M_upper;
         M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
            
        
        NAP=find_NAP(GSCP,M);

        A(Avar_ind)=Y(k,N_Xvar+1:N_Xvar+N_Avar); % size values
        [A(Avar_ind),~,~]=round_A_W(A(Avar_ind),M(Avar_ind),section_no,0.5*ones(1,numel(Avar_ind))); % Stochastically round the continuous size values to the closest upper/lower value in the available set of sections   
        A(sym_A(2:end,:))= repmat(A(sym_A(1,:)),size(sym_A,1)-1,1); % Assign the size values of dependent members
       % find out which members and shape variables are active in the sampled design   (independent and dependent)       
        A_active_ind=find(M); % active members 
        A_passive_ind=setdiff(1:N_member,A_active_ind);  % passive members (independent and dependent)
        X_active_ind=D*find(NAP);
        if D==3
            X_active_ind=sort([X_active_ind X_active_ind-1 X_active_ind-2]);
        elseif D==2
            X_active_ind=sort([X_active_ind X_active_ind-1]);
        end
        X_passive_ind=setdiff(1:N_node*D,X_active_ind); % passive shape variables (independent and dependent)
        keep_M(k,:)=M; % store the topology of the design
         [~,q1]=ismember(intersect(X_active_ind,Xvar_ind),Xvar_ind) ; % find index of independent active coordinates
         [~,q2]=ismember(intersect(A_active_ind,Avar_ind),Avar_ind);
        Y_active_ind=[ q1 q2+N_Xvar  (1:N_Mvar)+N_Xvar+N_Avar]; % index of active variables in the k-th solution (Y(k,:))
        Y_W(k,:)=Y(k,:)*0;
        Y_W(k,Y_active_ind)=1; % store which variables in the k-th solution were active (required for excluding the effect of the passive variables in recombination)
        %   Now update the the k-th solution
%         Y(k,1:N_Xvar)=X(Xvar_ind);
        Y(k,N_Xvar+1:N_Xvar+N_Avar)=A(Avar_ind);
        Y(k,N_Xvar+N_Avar+1:end)=M(Mvar_ind);
        Z(k,:)=(Y(k,:)-YMEAN)/S(k); % perturbation vector 
        
        X = X_const_val; % no shape analysis
        
        %  analyze the design 
        [Length_M,Vol,displ_ratio,buck_ratio,   stress_ratio,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_simp(NAP,X,DOF_constrained,A,   M,GSCP,Fext,M_elasticity,SigT_all,SigC_all,dis_all,Max_Kcond,D); 
        counteval=counteval+1; 
        
        %  assign the highest constraint violation to all the coupled variables
        buck_ratio = zeros(1,N_member); % 10 bar
        stress_ratio(sym_A)= repmat(max(stress_ratio(sym_A)),size(sym_A,1),1);
        %  Estimate the required increase in the cross sections such that all constraints are satisfied
        Agoal0=A;
        Agoal2=Agoal0.*max(1,max([stress_ratio;buck_ratio])); % Agoal2>=Agoal0    
        Agoal2=round_A_W(Agoal2,M,section_no,ones(1,N_member));
        Agoal3=Agoal0.*max(displ_ratio); % proportional increase for satisfaction of displacement constriants
        Agoal= max([Agoal0+(Agoal2-Agoal0).*Pcoeff;Agoal3]); % The estimated required increase + the initial cross section area
        f(k)=Density*(  Vol+sum((Agoal-A).*M.*Length_M)   ); % the value of the objective function (penalized weight)
        constraint_vio(k,:)=M.*max(0,max([stress_ratio;buck_ratio;max(displ_ratio)*ones(1,N_member)])-1); % contraint violation 
        if (max(constraint_vio(k,:))<=opt.feastol)  && (f(k)<best_feas_design(end)) % update the best feasible solution found so far
            best_feas_design=[X A M f(k)];
        end
        % Now resize. Pre-allocate first
        Y(k+lambda,:)=Y(k,:);  
        A_resized=A;  
        S(k+lambda)=S(k); 
        Y_W(k+lambda,:)=Y_W(k,:);  
        keep_M(k+lambda,:)=M;
        keep_M(k+lambda,:)=XXX; % Upper level
        if max(abs(stress_ratio))==0 % the k-th design was kinematically unstable, no resizing is possible
            f(k+lambda)=f(k);
            Z(k+lambda,:)= Z(k,:);
            constraint_vio(k+lambda,:)= constraint_vio(k,:);
        else % the k-th design was kinematically stable, perform resizing
            A_resized=           resize_FSDII_simp(Length_M,dis_all,buck_ratio,stress_ratio,f_int_ext,displacement,M,A,      M_elasticity,Pcoeff,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,opt.ucr_red);
            Y(k+lambda,N_Xvar+1:N_Xvar+N_Avar)=A_resized(Avar_ind); % Update the resized solution
            Z(k+lambda,:)=(Y(k+lambda,:)-YMEAN)/S(k+lambda);   % Update the corresponding perturbation vector 
            [Length_M,Vol,displ_ratio,buck_ratio,   stress_ratio,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_simp(NAP,X,DOF_constrained,A_resized,M,GSCP,Fext,M_elasticity,SigT_all,SigC_all,dis_all,Max_Kcond,D); 
            counteval=counteval+1;
            % assign the most critical constraint to all coupled sections    
            buck_ratio= sqrt(SF_buck*buck_ratio);
            buck_ratio(sym_A)= repmat(max(buck_ratio(sym_A)),size(sym_A,1),1);
            stress_ratio(sym_A)= repmat(max(stress_ratio(sym_A)),size(sym_A,1),1);
            % estimate the required increase in cross sections such that all constraints are satisfied
            Agoal0=A_resized;
            Agoal2=Agoal0.*max(1,max([stress_ratio;buck_ratio])); % Agoal2>=Agoal0    
            Agoal2=round_A_W(Agoal2,M,section_no,ones(1,N_member));
            Agoal3=Agoal0.*max(displ_ratio); % proportional increase for satisfaction of displacement constriants
            Agoal= max([Agoal0+(Agoal2-Agoal0).*Pcoeff;Agoal3]); % The estimated required increase + the initial cross section area
            f(k+lambda)=Density*(  Vol+sum((Agoal-A_resized).*M.*Length_M)   ); % the value of the objective function (penalized weight)
            constraint_vio(k+lambda,:)=M.*max(0,max([stress_ratio;buck_ratio;max(displ_ratio)*ones(1,N_member)])-1); % constraint violation of the resized solution
            if (max(constraint_vio(k+lambda,:))<=opt.feastol)  && (f(k+lambda)<best_feas_design(end)) % update the best identified feasible solution
                best_feas_design=[X A_resized M f(k+lambda)];
            end
        end % resizing completed
     end %  lambda+lambda solutions were created and evaluated
    % Now perform recombination and update strategy parameters
%     best_feas_design
    resize_eff=mean(sign(1-f(lambda+1:2*lambda) ./f(1:lambda))); %  resizing efficiency
    feas_rat=mean(max(constraint_vio')<=opt.feastol);%  ratio of feasible solutions
    Max_Avar=Max_Avar*exp(resize_eff*sqrt(Ts)); % update the move limit ratio
    Max_Avar=max(min(Max_Max_Avar,Max_Avar),Min_Max_Avar); % make sure the move limit ratio remains within the predefined range  
    % Recombination of design parameters
    [f,ind]=sort(f); % sort solutions
    activeness=weights*Y_W(ind(1:mu),:); % weighted fraction of parents in which an arbitrary variable is active
    YMEAN_update1=YMEAN+(weights*(Y(ind(1:mu),:).*Y_W(ind(1:mu),:))-mean(Y.*Y_W));  % Update based on comparing a variable in the whole population and the selected parents. This type of update is used for topology variables
    YMEAN_update2=(1-activeness).*YMEAN+weights*(Y(ind(1:mu),:).*Y_W(ind(1:mu),:)); % Update based on direct values in the parents. This type of update is used for shape and size variables
    YMEAN(1:N_Xvar+N_Avar)=YMEAN_update2(1:N_Xvar+N_Avar);
    YMEAN(1+N_Xvar+N_Avar:end)=YMEAN_update1(1+N_Xvar+N_Avar:end); % updated recombinant design
    YMEAN=min(U_Y,max(D_Y,YMEAN)); % make sure the recombinant design remains inside the bounds
    % Recombination of the scaling factors 
    SMEAN=exp(weights*(log(S(ind(1:mu)))')); % update the global step size
    goodZ=Y_W(ind(1:mu),:).*(Z(ind(1:mu),:)); 
    suggD=weights*(goodZ.^2); % weighted average of the scaling factors in parental population
    randomD=mean((Y_W.*Z).^2); % average of the scaling factors in the whole population
    STR_update1=STR.^2+Tc*(suggD-randomD); % update rule for scaling factors of topology variables
    STR_update2=(1-Tc*weights*(Y_W(ind(1:mu),:).^2)).*STR.^2+Tc*suggD; % update rule for scaling factors of shape and size variables
    STR(1:N_Xvar+N_Avar)=STR_update2(1:N_Xvar+N_Avar);
    STR(1+N_Xvar+N_Avar:end)=STR_update1(1+N_Xvar+N_Avar:end); % updated scaling factors
    STR=sqrt(max(1e-14,STR)); % make sure the scaling factors remain positive
    Xmean(Xvar_ind)=YMEAN(1:N_Xvar);
    Amean(Avar_ind)=YMEAN(1+N_Xvar:N_Xvar+N_Avar);
    Mmean(Mvar_ind)=YMEAN(N_Xvar+N_Avar+1:end);
    keep_history(iter,:)=[iter counteval best_feas_design(end) f(1)  max(constraint_vio(ind(1),:)) max(Pcoeff) min(Pcoeff) resize_eff]; % history of convergence 
    % update penalty coefficients (Pcoeff)
    old_infeas_ratio=infeas_ratio; % the fraction of constriant violatring maembers in the previous iteration
    exclude_it=    [ keep_M(ind(1:mu),:)==0]  | [[f(ind(1:mu))'*ones(1,N_member)]>(1e99/Density) ];  % absent members or unstable designs, they must not have any effect on the update of penalty coefficients
    infeas_ratio=    weights*(       (1-exclude_it).* (constraint_vio(ind(1:mu),:)  >opt.feastol)   + opt.max_infeas*exclude_it   ); % The fraction of parents in which a specific member has violated a constraint
    update_Pcoeff=((infeas_ratio>=old_infeas_ratio) | (infeas_ratio<=opt.max_infeas)); % the condition for updating Pcoeff
    Pcoeff=Pcoeff.*exp( (infeas_ratio-opt.max_infeas)*sqrt(Ts) ).*update_Pcoeff + Pcoeff .* (1-update_Pcoeff); % update Pcoeff
    Pcoeff= max(minPcoeff,Pcoeff); % make sure the penalty coefficients remain larger than the minimum value

 if numel(stopping_record) == 0
     stopping_record(end+1) = best_feas_design(end);
     stopping_evals(end+1) = counteval;
 else
     if stopping_record(end) == best_feas_design(end)
         stopping_record(end+1) = best_feas_design(end);
         stopping_evals(end+1) = counteval;
     else
         stopping_record = [];
         stopping_evals = [];
         stopping_record(end+1) = best_feas_design(end);
         stopping_evals(end+1) = counteval;
     end
 end
 
end %run completed
bestX=best_feas_design(1:D*N_node);
bestA=best_feas_design(D*N_node+(1:N_member));
bestM=best_feas_design(D*N_node+N_member+(1:N_member));
bestNAP=find_NAP(GSCP,bestM);
final_eval = stopping_evals(1);
total_eval = counteval;
function [best_feas_design,final_eval,total_eval] = bar72(M_upper)
opt.ProblemName='72barDiscrete';
kk = 1;
D=3; %plannar or spatial
voroodi=[
16	2	2	2	1	0	0.3	1	1.5	2	1	0.5	100000
];
m = 1;
opt.resize_budget=voroodi(m,2);
opt.lambda_coeff=voroodi(m,3);
opt.Tccoeff=voroodi(m,4);
opt.initial_M=voroodi(m,5);
opt.ExploitSymmetry=voroodi(m,6);
opt.mulambdaratio=voroodi(m,7);
opt.iniSTRMult=voroodi(m,8:10);
opt.FSDBasedPenalty=voroodi(m,11);
opt.max_infeas=voroodi(m,12);
opt.BiasIniAroundTavan=voroodi(m,13);
opt.incAmeanTsTavan=1000000;
opt.TToplearn=1;
opt.DisConstMult=1;
opt.TerminateResizeNoChange=1;
opt.ratiosUnsucc(1:3)=[1 1 0];
opt.TruncatedCompensation=0;
opt.iniSMEAN=.25;
 % opt.applyDisConst=1;
opt.updateSTR=[2 2 1];
opt.resizeStrategy=2;
opt.GiveTimePcoeff=1;
opt.AdjustMoveLimit=1;
opt.UpdateLambda=0;
opt.UpdateTau=[0 0];
opt.updateZ=1;
opt.Tscoeff=1;
opt.maxeval_coeff=200;%200/opt.lambda_coeff;
opt.tavan=2;
opt.topology_efficiency_check=0;
opt.WriteInterval=120; % writing interval in seconds
opt.inconly=0;
opt.feastol=1e-7;
opt.c_sec_red=.05;
opt.specification = 'simplified';
opt.stopping_interval = 10;
 % ********************* Enter problem specific data here % **********************
[opt,D,N_node,N_member,GSCP,sym_M,sym_A,X_indep_ind,basic_node,basic_member,X_const_ind,X_const_val,DOF_constrained,Fext,M_elasticity,Density,dis_all,SigT_all, SigC_all,kappa,Fy,section_no,Min_A,Max_A,D_X, U_X, minL,maxL]=GetProblemData(opt);
% generating required data
if strcmp(opt.specification,'simplified')
    SF_buck=pi^2/kappa;
else
    SF_buck=1;
end

if strcmp(opt.ProblemName,'10barDiscrete') || strcmp(opt.ProblemName,'15barDiscrete') || strcmp(opt.ProblemName,'72barDiscrete') || strcmp(opt.ProblemName,'52barDiscrete') || strcmp(opt.ProblemName,'200barDiscrete')
    SF_buck=1e12;
end

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

% Find independent variables. Consider only variables that are independent
Xvar_ind=setdiff(X_indep_ind, X_const_ind); 
Avar_ind=setdiff(1:N_member,reshape(sym_A(2:end,:),1,numel(sym_A(2:end,:))));
Avar_ind=union(Avar_ind,sym_A(1,:)); 
Mvar_ind=setdiff(1:N_member,[basic_member,reshape(sym_A,1,numel(sym_A))]);
Mvar_ind=union(Mvar_ind,setdiff(sym_A(1,:),basic_member)); % independent topology variables
N_Xvar=numel(Xvar_ind); % Number of independent shape variables
N_Avar=numel(Avar_ind); % Number of independent size variables
N_Mvar=numel(Mvar_ind);  % Number of independent topology variables
U_Y=[U_X(Xvar_ind) Max_A*ones(1,N_Avar) ones(1,N_Mvar)]; 
D_Y=[D_X(Xvar_ind) Min_A*ones(1,N_Avar) zeros(1,N_Mvar)];
Max_Kcond=1e12; % The limit that specifies whether the stiffness matrix is singular
% Calculating effective number of design parameters
N_var_all=( sqrt(N_Mvar)+sqrt(N_Xvar)+sqrt(N_Avar))^2;
N_DefC=(N_node*D-numel(DOF_constrained))*(dis_all<1e10);
N_constraint_all=(sqrt(N_DefC)+sqrt(N_member))^2*sqrt(N_loadcase);

if 0
N_var_eff=(N_var_all)* sqrt(1+N_constraint_all/N_var_all); %efective number of variables
else
N_var_eff=N_var_all* sqrt(1+ sqrt(N_loadcase)*(N_member+D*N_node)/N_var_all   ); %efective number of variables
end

% Setting lambda, mu, maxiter
maxeval=round(opt.maxeval_coeff*N_var_eff);
target_eval_per_iter=round(2*opt.lambda_coeff*sqrt(N_var_eff));

Max_Max_Avar=(Max_A/Min_A)^1; % maximum value of the move limit ratio
if section_no>0
    Min_Max_Avar=max(availsec(2:end,1)./availsec(1:end-1,1)); % minimum value of the move limit ratio
else
    Min_Max_Avar=1.01;
end
minPcoeff=1;  % minimum value of penalty coefficients
infeas_ratio=opt.max_infeas*ones(1,N_member); % maximum fraction of members that may violate a constraint

    Xmean=(D_X+U_X)/2;
    Xmean(X_const_ind)=X_const_val;

        %%% UPPER LEVEL %%%%
    if strcmp(opt.ProblemName,'39barDiscrete')
        Mmean = M_upper;
        Mmean(sym_M(2:end,:))= repmat(Mmean(sym_M(1,:)),size(sym_M,1)-1,1);
    else
        Mmean = zeros(1,N_member);
        Mmean(sym_M(1,:)) = M_upper;
        Mmean(sym_M(2:end,:))= repmat(Mmean(sym_M(1,:)),size(sym_M,1)-1,1);
    end
    
    Amean=(Max_A+Min_A)/2*ones(1,N_member);
    YMEAN=([ Xmean(Xvar_ind) Amean(Avar_ind) Mmean(Mvar_ind)]); %Y is the set of all independent variables
    SMEAN=opt.iniSMEAN; % global step size
    STR=(U_Y-D_Y); % This is vector of scaling factors
    
    STR(1:N_Xvar)=STR(1:N_Xvar)*opt.iniSTRMult(1);
    STR((1:N_Avar)+N_Xvar)=STR((1:N_Avar)+N_Xvar)*opt.iniSTRMult(2);
    STR((1:N_Mvar)+N_Avar+N_Xvar)=STR((1:N_Mvar)+N_Avar+N_Xvar)*opt.iniSTRMult(3);
    
    iter=0; %iteration number
    Pcoeff=ones(1,N_member); % penalty coefficients for constriant violation of members 
    best_feas_design=zeros(1,2*numel(Mmean)+D*N_node+1); % best feasible solution and weight so far
    best_feas_design(end)=1e150;
    Max_Avar=Max_Max_Avar; % move limit ratio
    counteval=0; % evaluation count    
    lambda=round(target_eval_per_iter/(opt.resize_budget+1)+.5);
    keep_best=[];
% end
    mu=max(round((lambda*opt.mulambdaratio)),1);
    weights=log(1+mu)-log(1:mu);
    weights=weights/sum(weights); % weights for the recombination of the selected parents
    mueff=sum(weights)^2/sum(weights.^2);
    Tc=(1+opt.Tccoeff*N_var_eff/mueff)^(-1);   % learning rate of STR
    Ts=1/sqrt(2*opt.Tscoeff*N_var_eff);   % learning rate of the global step size
zarib5=1;
format shortg
tic
TimeInterval=0;

if 0
zind=find(Xmean(3:3:end)>0);
xy=intersect(find(Xmean(1:3:end)>=19) ,find(Xmean(2:3:end)>=15.5))
nindep=intersect(xy,zind)

javab=zeros(4,numel(nindep))
javab(1,:)=nindep;
for k=nindep
    p=setdiff(zind,k)
    for m=p
        if (Xmean(3*k-2)==(Xmean(3*m-2)) & Xmean(3*k-1)==(31-Xmean(3*m-1)))  |   (Xmean(3*k-2)==(38-Xmean(3*m-2)) & Xmean(3*k-1)==(31-Xmean(3*m-1)))  |   (Xmean(3*k-2)==(38-Xmean(3*m-2)) & Xmean(3*k-1)==(Xmean(3*m-1)))
            
            
            
            
            ind=find(javab(1,:)==k);
            bood=sum(javab(:,ind)>0);
            javab(bood+1,ind)=m;
        end
    end
 
end
javab
 
   
end

stopping = false;
stopping_record = [];
stopping_evals = [];

while counteval<maxeval && stopping == false % main optimization loop starts here
    iter=iter+1; old_counteval=counteval;
    
    % presetting values for speed boost
 

    if opt.resizeStrategy==2
        keep_M=zeros(2*lambda,N_member);
        Y=zeros(2*lambda,numel(STR)); 
        Y_W=Y;
        f=1e180*ones(1,2*lambda);
        constraint_vio=zeros(2*lambda,N_member);
        Z=zeros(2*lambda,numel(STR));
        S=SMEAN*ones(1,2*lambda); 
    else
        keep_M=zeros(1*lambda,N_member);
        Y=zeros(1*lambda,numel(STR)); 
        Y_W=Y;
        f=1e180*ones(1,1*lambda);
        constraint_vio=zeros(1*lambda,N_member);
        Z=zeros(1*lambda,numel(STR));
        S=SMEAN*ones(1,1*lambda);
        suggMax_Avar=Max_Avar*ones(1,lambda);
    end
    for k=1:lambda %start sampling solutions
        stable=0;
        while ~stable % Stay in this loop untill a solution that satisfies necessary conditions for stability (and efficiency, if applied) is generated
            S(k)=SMEAN*exp(Ts*randn); % mutate the global step size
            % now use the truncated normal distribution to sample a new design between D_Y and U_Y
            mincdf=normcdf(D_Y,YMEAN,S(k)*STR);
            maxcdf=normcdf(U_Y,YMEAN,S(k)*STR);
            r=max(1e-8,min(1-1e-8,rand(1,numel(U_Y)))); % we use 1e-8 bound to avoid infinity error
            Y(k,:)=norminv(mincdf+(maxcdf-mincdf).*r,YMEAN,S(k)*STR); %  the k-th candidate design with continuous values
            
            X=Xmean;A=Amean;M=Mmean;
            
            M(Mvar_ind)=(Y(k,N_Xvar+N_Avar+1:end)>rand(1,N_Mvar)); % the topology is defined in M
            M(basic_member)=1; % enforce presence of basic members, if any
            X(Xvar_ind)=Y(k,1:N_Xvar);
            X=enforce_shape_sym(X,N_member,opt);%using shape symmetry
            L_Link=calc_L(M,X,GSCP,D);
            M=M.*(L_Link>=minL).*(L_Link<=maxL);
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1); % deteremine absence/presence of topologically coupled members
          
            %%% UPPER LEVEL %%%%
            M = zeros(1,N_member);
            M(sym_M(1,:)) = M_upper;
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
                
            
            
            
            
            [NAP,NAPconnector]=find_NAP(GSCP,M); % find which nodes are active, and the number of the members connected to them
            NAPconnector_all=NAPconnector+DOFonNAP; % The number of members connected to each node + the number of reactions on the nodes
            chk1=(min(NAP(basic_node)==1)>.5); % necessary condition: Are all basic nodes active? 
            chk2=(sum(M) >= (sum(NAP)*D-numel(DOF_constrained))); % Necessary condition: Does the candidate truss have the minimum required number of members?
            chk3= min( max([ NAPconnector_all>=D ; NAP==0])      ); % Necessary condition: Are sufficient number of members connected to each node? 
              % [iter k chk1   chk2   chk3]
             if (chk1 && chk2 && chk3)
                stable=1; % The candidate design is accepted if it satisfies all the necessary conditions for stability
                if 0
                    for l2=1:size(javab,2)
                    l2
                    NAP=ones(1,N_node);
                    NAP(javab(:,l2))=0;
                    figure(l2)
                    plot_truss(GSCP,X,NAP,M,A,D);
                    end
                    fffff
                end  
                             
 
             end
        end % The k-th solution that satisfies the necessary conditions for stability is generated
       
        A(Avar_ind)=Y(k,N_Xvar+1:N_Xvar+N_Avar); % size values
        Ri=zeros(1,N_member); % Radii of gyration of sections, presetting
        [A(Avar_ind),Ri(Avar_ind),~]=round_A_W(A(Avar_ind),M(Avar_ind),section_no,(1-0.5*exp(1-Pcoeff(Avar_ind).^opt.BiasIniAroundTavan)).*ones(1,N_Avar)); % Stochastically round the continuous size values to the closest upper/lower value in the available set of sections   
        A(sym_A(2:end,:))= repmat(A(sym_A(1,:)),size(sym_A,1)-1,1); % Assign the size values of dependent members
        Ri(sym_A(2:end,:))= repmat(Ri(sym_A(1,:)),size(sym_A,1)-1,1); %Assign the radius of gyration of dependent memebrs
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
        X(X_const_ind)=X_const_val; % set the fixed coordinates to the given value
        
        %%% UPPER LEVEL %%%%
        if strcmp(opt.ProblemName,'39barDiscrete')
            M = M_upper;
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
        else
            M = zeros(1,N_member);
            M(sym_M(1,:)) = M_upper;
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
        end
            
            
        keep_M(k,:)=M; % store the topology of the design
        
        %  find out which independent shape and size variables are active in the k-th solution (Y(k,:))  
        [~,q1]=ismember(intersect(X_active_ind,Xvar_ind),Xvar_ind) ; % find index of independent active coordinates
        [~,q2]=ismember(intersect(A_active_ind,Avar_ind),Avar_ind);
        Y_active_ind=[ q1 q2+N_Xvar  (1:N_Mvar)+N_Xvar+N_Avar]; % index of active variables in the k-th solution (Y(k,:))
        Y_W(k,:)=Y(k,:)*0;
        Y_W(k,Y_active_ind)=1; % store which variables in the k-th solution were active (required for excluding the effect of the passive variables in recombination)
        %   Now update the the k-th solution

        resized_iter=0;
        suggMax_Avar(k)=Max_Avar;
       
        while (opt.resizeStrategy==1 && resized_iter<=opt.resize_budget) | resized_iter==0
             if strcmp(opt.specification,'AISC-ASD')
                [Length_M,Vol,displ_ratio,slender_ratio,stress_ratio,f_int_ext,displacement,f_int_unit]= FE_solve_2D3D_ASD(NAP,X,DOF_constrained,A,Ri,M,GSCP,Fext,M_elasticity,Fy,dis_all,Max_Kcond,D); 
            elseif strcmp(opt.specification,'simplified')
                [Length_M,Vol,displ_ratio,slender_ratio,stress_ratio,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_simp(NAP,X,DOF_constrained,A,   M,GSCP,Fext,M_elasticity,SigT_all,SigC_all,dis_all,Max_Kcond,D); 
                slender_ratio= sqrt(SF_buck*slender_ratio);
            else
                'Error, correct option for specification was not provided. Aborting ....'
                ABORT
            end
            counteval=counteval+1; 
           %  assign the highest constraint violation to all the coupled variables
           if opt.FSDBasedPenalty==1
                slender_ratio(sym_A)= repmat(max(slender_ratio(sym_A)),size(sym_A,1),1);
                stress_ratio(sym_A)= repmat(max(stress_ratio(sym_A)),size(sym_A,1),1);
                %  Estimate the required increase in the cross sections such that all constraints are satisfied
                slender_ratio = zeros(1,N_member); % No buckling
                Agoal0=A;
                Rigoal0=Ri;
                if strcmp(opt.specification,'AISC-ASD')
                    f_int_ext_increased=f_int_ext;%f_int_ext.*repmat(max(1,stress_ratio.^(opt.tavanAgoal+Pcoeff.^opt.TavanApplyPcoeffInAgoal-1)),N_loadcase,1)';
                 
                    Agoal2= find_Agoal_ASD(Length_M,f_int_ext_increased,M,Agoal0,Rigoal0,Fy,M_elasticity,section_no,sym_A,Avar_ind,opt); % individual section area increase for satisfaction of member-based constraints 
                   % Agoal2=find_Agoal_ASD2(Length_M,f_int_ext_increased,M,Agoal0,Rigoal0,Fy,M_elasticity,section_no,sym_A,Avar_ind,max([slender_ratio;stress_ratio]),opt); % individual section area increase for satisfaction of member-based constraints 
                 
                
                elseif strcmp(opt.specification,'simplified')
                    Agoal2=Agoal0.*max(1,max([stress_ratio;slender_ratio])); % Agoal2>=Agoal0    
                    Agoal2=max([ Agoal2;round_A_W(Agoal2,M,section_no,ones(1,N_member))]);
                end
                Agoal3=Agoal0.*max(displ_ratio); % proportional increase for satisfaction of displacement constriants
                Agoal= max([Agoal0+Pcoeff.*(Agoal2-Agoal0);Agoal3]); % The estimated required increase + the initial cross section area
                f(k)=Density*(  Vol+sum((Agoal-A).*M.*Length_M)   ); % the value of the objective function (penalized weight)
           else
                jarime=sum((stress_ratio.^2-1).*(stress_ratio>1))  +  sum((slender_ratio.^2-1).*(slender_ratio>1))  +  sum((displ_ratio.^2-1).*(displ_ratio>1));  
                f(k)=Density*(Vol)*(1+jarime); % the value of the objective function (penalized weight)
           end
           
                constraint_vio(k,:)=M.*max(0,max([stress_ratio.^opt.tavan;slender_ratio;max(displ_ratio).^opt.tavan*ones(1,N_member)])-1); % contraint violation 
            if max(abs(stress_ratio))==0 
                break
            end
            if (max(constraint_vio(k,:))<=opt.feastol)  && (f(k)<best_feas_design(end)) % update the best feasible solution found so far
                best_feas_design=[X A M f(k)];
                keep_best=[keep_best; [counteval best_feas_design(end)]];
            end
            if resized_iter>=opt.resize_budget | opt.resizeStrategy==2
                break
            end
            if resized_iter>0
                if oldf(k)<=f(k)
                    A=oldA;
                    Ri=oldRi;
                    Vol=oldVol;
                    displ_ratio=olddispl_ratio;
                    slender_ratio=oldslender_ratio;
                    stress_ratio=oldstress_ratio;
                    f_int_ext=oldf_int_ext;
                    displacement=olddisplacement;
                    f_int_unit=olff_int_unit;
                    suggMax_Avar(k)=suggMax_Avar(k)^(1-Ts);
                    f=oldf;
                    constraint_vio(k,:)=M.*max(0,max([stress_ratio.^opt.tavan;slender_ratio;max(displ_ratio).^opt.tavan*ones(1,N_member)])-1); % contraint violation 
                else
                    suggMax_Avar(k)=suggMax_Avar(k)^(1+Ts);
                end
                suggMax_Avar(k)=max(min(Max_Max_Avar,suggMax_Avar(k)),Min_Max_Avar); % make sure the move limit ratio remains within the predefined range  
                %[k resized_iter oldf(k) f(k) max(constraint_vio(k,:)) suggMax_Avar(k)]
            end
            if (opt.inconly==0) | (resized_iter==1)
                 if strcmp(opt.specification,'AISC-ASD')  
                     [A_resized,Ri_resized]=resize_FSDII_ASD(Length_M,dis_all,f_int_ext,displacement,M,A,Ri,Fy,M_elasticity,Pcoeff,suggMax_Avar(k),section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);  % returns the resized section
                 elseif strcmp(opt.specification,'simplified')
                     A_resized=            resize_FSDII_simp(Length_M,dis_all,slender_ratio,stress_ratio,f_int_ext,displacement,M,A,M_elasticity,Pcoeff,suggMax_Avar(k),section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);
                 end
            elseif opt.inconly==1 
                 if strcmp(opt.specification,'AISC-ASD')  
                     [A_resized,Ri_resized]=resize_FSDII_ASD(Length_M,dis_all,f_int_ext,displacement,M,A,Ri,Fy,M_elasticity,Pcoeff*0+25,suggMax_Avar(k),section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);  % returns the resized section
                 elseif strcmp(opt.specification,'simplified')
                     A_resized=            resize_FSDII_simp(Length_M,dis_all,slender_ratio,stress_ratio,f_int_ext,displacement,M,A,M_elasticity,Pcoeff*0+25,suggMax_Avar(k),section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);
                 end
            else
                 'Unsupported option for the resizing. Aborting the optimization process ...'
                 abort
            end
            resized_iter=resized_iter+1;
            if max(abs(A_resized-A))==0
                break
            else
            oldA=A;
            oldRi=Ri;
            oldVol=Vol;
            olddispl_ratio=displ_ratio;
            oldslender_ratio=slender_ratio;
            oldstress_ratio=stress_ratio;
            oldf_int_ext= f_int_ext;
            olddisplacement=displacement;
            olff_int_unit=f_int_unit;
            oldf=f;  
            A=A_resized;
            Ri=Ri_resized;
            end
        end
        if max(constraint_vio(k,:))>0 & max(A)<Max_A & opt.resize_budget>3 & opt.resizeStrategy==1
                  [k resized_iter oldf(k) f(k) max(constraint_vio(k,:)) suggMax_Avar(k)]
                  'warning, we could not find a fesible solution?'
        end   
        
        
        
        Y(k,1:N_Xvar)=X(Xvar_ind);
        Y(k,N_Xvar+1:N_Xvar+N_Avar)=A(Avar_ind);
        Y(k,N_Xvar+N_Avar+1:end)=M(Mvar_ind);
        
        if opt.updateZ==1
            thisYk=Y(k,:);
        end
                thisYk=Y(k,:);
        truncatadedVar=calc_var(D_Y,U_Y,YMEAN,S(k)*STR);
        if opt.TruncatedCompensation==0
            Z(k,:)=(thisYk-YMEAN)/S(k); % perturbation vector 
        elseif opt.TruncatedCompensation==1;
            Z(k,:)=(thisYk-YMEAN)./sqrt(truncatadedVar);
        end
        
         %%% UPPER LEVEL %%%%
        if strcmp(opt.ProblemName,'39barDiscrete')
            M = M_upper;
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
        else
            M = zeros(1,N_member);
            M(sym_M(1,:)) = M_upper;
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1);
        end
      
        %  analyze the design 
               % Now resize. Pre-allocate first
        if opt.resizeStrategy==2      
        Y(k+lambda,:)=Y(k,:);  
        A_resized=A;  
        S(k+lambda)=S(k); 
        Y_W(k+lambda,:)=Y_W(k,:);  
        keep_M(k+lambda,:)=M;
        end
        if max(abs(stress_ratio))==0 & opt.resizeStrategy==2   | (opt.resize_budget==0 & opt.resizeStrategy==2)    % the k-th design was kinematically unstable, no resizing is possible
            f(k+lambda)=f(k);
            Z(k+lambda,:)= Z(k,:);
            constraint_vio(k+lambda,:)= constraint_vio(k,:);
        elseif opt.resizeStrategy==2  % the k-th design was kinematically stable, perform resizing
            A_resized=A;
            Ri_resized=Ri;
            resized_iter=1;
     
            while resized_iter<=opt.resize_budget && (max(abs(stress_ratio))>0)
                 old_A_resized=A_resized;
                 if (opt.inconly==0) | (resized_iter==1)
                     if strcmp(opt.specification,'AISC-ASD')  
                         [A_resized,Ri_resized]=resize_FSDII_ASD(Length_M,dis_all,f_int_ext,displacement,M,A_resized,Ri_resized,Fy,M_elasticity,Pcoeff,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);  % returns the resized section
                     elseif strcmp(opt.specification,'simplified')
                          A_resized=            resize_FSDII_simp(Length_M,dis_all,slender_ratio,stress_ratio,f_int_ext,displacement,M,A_resized,M_elasticity,Pcoeff,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);
                          
                     end
                 elseif opt.inconly==1 
                     if strcmp(opt.specification,'AISC-ASD')  
                         [A_resized,Ri_resized]=resize_FSDII_ASD(Length_M,dis_all,f_int_ext,displacement,M,A_resized,Ri_resized,Fy,M_elasticity,Pcoeff*0+25,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);  % returns the resized section
                     elseif strcmp(opt.specification,'simplified')
                          A_resized=            resize_FSDII_simp(Length_M,dis_all,slender_ratio,stress_ratio,f_int_ext,displacement,M,A_resized,M_elasticity,Pcoeff*0+25,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,opt.c_sec_red);
                     end
                 else
                     'Unsupported option for the resizing. Aborting the optimization process ...'
                     abort
                 end
                 resized_iter=resized_iter+1;
                 if norm(A_resized-old_A_resized)<1e-12 & opt.TerminateResizeNoChange
                     break;
                 end
                 if strcmp(opt.specification,'AISC-ASD')
                      [Length_M,Vol,displ_ratio,slender_ratio,stress_ratio,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_ASD(NAP,X,DOF_constrained,A_resized,Ri_resized,M,GSCP,Fext,M_elasticity,Fy,dis_all,Max_Kcond,D); 
                 elseif strcmp(opt.specification,'simplified')
                      [Length_M,Vol,displ_ratio,slender_ratio,stress_ratio,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_simp(NAP,X,DOF_constrained,A_resized,   M,GSCP,Fext,M_elasticity,SigT_all,SigC_all,dis_all,Max_Kcond,D); 
                      slender_ratio= sqrt(SF_buck*slender_ratio);
                 end
                
                 counteval=counteval+1;
                
            end % resized
             if max([displ_ratio,slender_ratio,stress_ratio])>1 &  max(A_resized)<Max_A & opt.resize_budget>4
                     max([displ_ratio,slender_ratio,stress_ratio])
                     'Warning! why we could not reach a feasible design?'
                 end
               Y(k+lambda,N_Xvar+1:N_Xvar+N_Avar)=A_resized(Avar_ind); % Update the resized solution
                 if opt.updateZ==1
                     thisYk=Y(k+lambda,:);
                 end
                 if opt.TruncatedCompensation==0
                     Z(k+lambda,:)=(thisYk-YMEAN)/S(k); % perturbation vector 
                 elseif opt.TruncatedCompensation==1;
                     Z(k+lambda,:)=(thisYk-YMEAN)./sqrt(truncatadedVar);
                 end
                 if opt.FSDBasedPenalty==1
                     slender_ratio(sym_A)= repmat(max(slender_ratio(sym_A)),size(sym_A,1),1);
                     stress_ratio(sym_A)= repmat(max(stress_ratio(sym_A)),size(sym_A,1),1); %important for update of infeas_ratio    
                     Agoal0=A_resized;
                     Rigoal0=Ri_resized;
                     if strcmp(opt.specification,'AISC-ASD')
                         f_int_ext_increased=f_int_ext;%f_int_ext.*repmat(max(1,stress_ratio.^(opt.tavanAgoal+Pcoeff.^opt.TavanApplyPcoeffInAgoal-1)),N_loadcase,1)';
                         Agoal2=find_Agoal_ASD(Length_M,f_int_ext_increased,M,Agoal0,Rigoal0,Fy,M_elasticity,section_no,sym_A,Avar_ind,opt); % individual section area increase for satisfaction of member-based constraints 
  
                     elseif strcmp(opt.specification,'simplified')
                         Agoal2=Agoal0.*max(1,max([stress_ratio;slender_ratio])); % Agoal2>=Agoal0    
                         Agoal2=max([ Agoal2;round_A_W(Agoal2,M,section_no,ones(1,N_member))]);
                     end
                     Agoal3=Agoal0.*max(displ_ratio); % proportional increase for satisfaction of displacement constriants
                     Agoal= max([Agoal0+Pcoeff.*(Agoal2-Agoal0);Agoal3]); % The estimated required increase + the initial cross section area
                     f(k+lambda)=Density*(  Vol+sum((Agoal-A_resized).*M.*Length_M)   ); % the value of the objective function (penalized weight)
                 else
                     jarime=sum((stress_ratio.^2-1).*(stress_ratio>1))  +  sum((slender_ratio.^2-1).*(slender_ratio>1))  +  sum((displ_ratio.^2-1).*(displ_ratio>1));  
                     f(k+lambda)=Density*(Vol)*(1+jarime); % the value of the objective function (penalized weight)
                 end
           
                 
                 constraint_vio(k+lambda,:)=M.*max(0,max([stress_ratio.^opt.tavan;slender_ratio;max(displ_ratio).^opt.tavan*ones(1,N_member)])-1); % constraint violation of the resized solution
                 if (max(constraint_vio(k+lambda,:))<=opt.feastol)  && (f(k+lambda)<best_feas_design(end)) % update the best identified feasible solution
                     best_feas_design=[X A_resized M f(k+lambda)];
                     keep_best=[keep_best; [counteval best_feas_design(end)]];
                 end

             end % resizing completed
    end
    % Now perform recombination and update strategy parameters
    if opt.resizeStrategy==2
        resize_eff=mean(sign(1-f(lambda+1:2*lambda) ./f(1:lambda))); %  resizing efficiency
        feas_rat=mean(max(constraint_vio')<=opt.feastol);%  ratio of feasible solutions
        if opt.AdjustMoveLimit==1
            Max_Avar=Max_Avar*exp(resize_eff*sqrt(Ts)); % update the move limit ratio
            Max_Avar=max(min(Max_Max_Avar,Max_Avar),Min_Max_Avar); % make sure the move limit ratio remains within the predefined range  
        end
    else
         if opt.AdjustMoveLimit==1
            Max_Avar=prod(suggMax_Avar.^(1/lambda));
            Max_Avar=max(min(Max_Max_Avar,Max_Avar),Min_Max_Avar); % make sure the move limit ratio remains within the predefined range  
         end
    end
    % Recombination of design parameters
    oldf=f;
    [f,ind]=sort(f); % sort solutions
    activeness=weights*Y_W(ind(1:mu),:); % weighted fraction of parents in which an arbitrary variable is active
    YMEAN_update1=YMEAN+(weights*(Y(ind(1:mu),:).*Y_W(ind(1:mu),:))-mean(Y.*Y_W));  % Update based on comparing a variable in the whole population and the selected parents. This type of update is used for topology variables
    if ~(opt.TToplearn==1)
        ind4=N_Xvar+N_Avar+(1:N_Mvar);
        YMEAN_update1(ind4)=YMEAN(ind4)+opt.TToplearn*(weights*(Y(ind(1:mu),ind4).*Y_W(ind(1:mu),ind4))-mean(Y(:,ind4).*Y_W(:,ind4)));  % Update based on comparing a variable in the whole population and the selected parents. This type of update is used for topology variables
    end
    YMEAN_update2=(1-activeness).*YMEAN+weights*(Y(ind(1:mu),:).*Y_W(ind(1:mu),:)); % Update based on direct values in the parents. This type of update is used for shape and size variables
    YMEAN(1:N_Xvar+N_Avar)=YMEAN_update2(1:N_Xvar+N_Avar);
    YMEAN(1+N_Xvar+N_Avar:end)=YMEAN_update1(1+N_Xvar+N_Avar:end); % updated recombinant design
    YMEAN=min(U_Y,max(D_Y,YMEAN)); % make sure the recombinant design remains inside the bounds
    % Recombination of the scaling factors 
    SMEAN=exp(weights*(log(S(ind(1:mu)))')); % update the global step size
    goodZ=Y_W(ind(1:mu),:).*(Z(ind(1:mu),:)); 
    suggD=weights*(goodZ.^2); % weighted average of the scaling factors in parental population
    randomD=mean((Y_W.*Z).^2); % average of the scaling factors in the whole population
    oldSTR=STR;
    STR_update1=STR.^2+Tc*(suggD-randomD); % update rule for scaling factors of topology variables
    STR_update11=STR_update1;
    if ~(opt.TToplearn==1)
        STR_update1(ind4)=STR(ind4).^2+Tc*(suggD(ind4)-randomD(ind4)); % update rule for scaling factors of topology variables
    end
    STR_update2=(1-Tc*weights*(Y_W(ind(1:mu),:).^2)).*STR.^2+Tc*suggD; % update rule for scaling factors of shape and size variables
    
    if opt.updateSTR(1)==1
        STR(1:N_Xvar)=STR_update1(1:N_Xvar);
    end
    if opt.updateSTR(2)==1
        STR((1:N_Avar)+N_Xvar)=STR_update1((1:N_Avar)+N_Xvar);
    end
    if opt.updateSTR(3)==1
        STR((1:N_Mvar)+N_Xvar+N_Avar)=STR_update1((1:N_Mvar)+N_Xvar+N_Avar);
    end
    if opt.updateSTR(1)==2
        STR(1:N_Xvar)=STR_update2(1:N_Xvar);
    end
    if opt.updateSTR(2)==2
        STR((1:N_Avar)+N_Xvar)=STR_update2((1:N_Avar)+N_Xvar);
    end
    if opt.updateSTR(3)==2
        STR((1:N_Mvar)+N_Xvar+N_Avar)=STR_update2((1:N_Mvar)+N_Xvar+N_Avar);
    end
    STR=sqrt(max(1e-14,STR)); % make sure the scaling factors remain positive
    Xmean(Xvar_ind)=YMEAN(1:N_Xvar);
    Amean(Avar_ind)=YMEAN(1+N_Xvar:N_Xvar+N_Avar);
    Mmean(Mvar_ind)=YMEAN(N_Xvar+N_Avar+1:end);
      % update penalty coefficients (Pcoeff)
     
        old_infeas_ratio=infeas_ratio; % the fraction of constriant violatring maembers in the previous iteration
        exclude_it=    [ keep_M(ind(1:mu),:)==0]  | [[f(1:mu)'*ones(1,N_member)]>(1e99/Density) ];  % absent members or unstable designs, they must not have any effect on the update of penalty coefficients
        infeas_ratio=    weights*(       (1-exclude_it).* (constraint_vio(ind(1:mu),:)  >opt.feastol)   + opt.max_infeas*exclude_it   ); % The fraction of parents in which a specific member has violated a constraint

        if opt.GiveTimePcoeff==1
             update_Pcoeff=((infeas_ratio>=old_infeas_ratio) | (infeas_ratio<=opt.max_infeas)); % the condition for updating Pcoeff
             Pcoeff=Pcoeff.*exp( (infeas_ratio-opt.max_infeas)*sqrt(Ts) ).*update_Pcoeff + Pcoeff .* (1-update_Pcoeff); % update Pcoeff
        elseif opt.GiveTimePcoeff==0
             update_Pcoeff=ones(1,N_member);
             Pcoeff=Pcoeff.*exp( (infeas_ratio-opt.max_infeas)*sqrt(Ts) ).*update_Pcoeff + Pcoeff .* (1-update_Pcoeff); % update Pcoeff
        elseif opt.GiveTimePcoeff==.5
             update_Pcoeff=((infeas_ratio>=old_infeas_ratio) | (infeas_ratio<=opt.max_infeas)); % the condition for updating Pcoeff
             new_infeas_ratio=infeas_ratio;
             infeas_ratio=infeas_ratio.*update_Pcoeff + (1-update_Pcoeff).*opt.max_infeas;
             Pcoeff=Pcoeff.*exp( (infeas_ratio-opt.max_infeas)*(sqrt(Ts)) ) ;  
        end
        Pcoeff= max(minPcoeff,Pcoeff); % make sure the penalty coefficients remain larger than the minimum value
     
    
    Nvar_top_active=sum(var(Y(:,N_Xvar+N_Avar+(1:N_Mvar)))>0);
    Nvar_shape_active=sum(mean(Y_W(:,1:N_Xvar))); 
    Nvar_size_active=sum(mean(Y_W(:,(1:N_Avar)+N_Xvar))); 
    N_var_all_active=( sqrt(Nvar_top_active)+sqrt(Nvar_shape_active)+sqrt(Nvar_size_active))^2;
    N_var_eff=(N_var_all_active)* sqrt(1+N_constraint_all/N_var_all); %efective number of variables
    if opt.UpdateLambda==1
        lambda=round(opt.lambda_coeff*sqrt(N_var_eff));  
        mu=max(round((lambda*opt.mulambdaratio)),1); 
        weights=log(1+mu)-log(1:mu);
        weights=weights/sum(weights); % weights for the recombination of the selected parents
        mueff=sum(weights)^2/sum(weights.^2);
    end
    if opt.UpdateTau==1
        Tc=(1+opt.Tccoeff*N_var_eff/mueff)^(-1);   % learning rate of STR
        Ts=1/sqrt(2*opt.Tscoeff*N_var_eff);   % learning rate of the global step size 
    end
     zarib5=(infeas_ratio(Avar_ind)-opt.max_infeas)*(Ts^opt.incAmeanTsTavan);
     zarib5=max(0,zarib5);
     Amean(Avar_ind)=Amean(Avar_ind).*exp(zarib5);
     Amean(Avar_ind)=min(max(Min_A,Amean(Avar_ind)),Max_A);
     YMEAN(N_Xvar+(1:N_Avar))=Amean(Avar_ind);
    keep_history(iter,:)=[iter counteval best_feas_design(end) f(1)  weights*max(constraint_vio(ind(1:mu),:)')' max(Pcoeff) min(Pcoeff)  Max_Avar N_var_eff lambda maxeval]; % history of convergence 
    
    usedeval_this_iteration=counteval-old_counteval;
    lambda=round(target_eval_per_iter/usedeval_this_iteration*lambda+.5);
    mu=max(round((lambda*opt.mulambdaratio)),1);
    weights=log(1+mu)-log(1:mu);
    weights=weights/sum(weights); % weights for the recombination of the selected parents
    mueff=sum(weights)^2/sum(weights.^2);
    Tc=(1+opt.Tccoeff*N_var_eff/mueff)^(-1);   % learning rate of STR
    Ts=1/sqrt(2*opt.Tscoeff*N_var_eff);   % learning rate of the global step size
    
bestX=best_feas_design(1:D*N_node);
bestA=best_feas_design(D*N_node+(1:N_member));
bestM=best_feas_design(D*N_node+N_member+(1:N_member));
bestNAP=find_NAP(GSCP,bestM);



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
 
 if numel(stopping_record)==opt.stopping_interval
     stopping = true;
 end

end %run completed
bestX=best_feas_design(1:D*N_node);
bestA=best_feas_design(D*N_node+(1:N_member));
bestM=best_feas_design(D*N_node+N_member+(1:N_member));
bestNAP=find_NAP(GSCP,bestM);
final_eval = stopping_evals(1);
total_eval = counteval;
end


% Please study the readme.pdf file before using this code. 
% This file optimizes the 224-bar pyramid problem.
% AISC-ASD (9-th edition) specifications govern the constraints.

function [best_feas_design,final_eval,total_eval] = bar224(M_upper)
kk = 1;

% ********** Control parameters of FSD-ES II (Default values) *************
opt.initial_M=1; % initial values for recombinant design for topology variables
opt.topology_efficiency_check=0; % check the necessary condition for efficiency of the sampled topology (deactivated)
opt.max_infeas=.5; % the threshold for ratio of constraint violating members. 
opt.WriteInterval=5; % Interval for writing the optimization process
opt.lambda_coeff=2; % coefficient for the population size
opt.mulambdaratio=.3; % ratio of parents to offspring
opt.maxiter_coeff=1; % coefficient for the maximum number of iterations
opt.Tccoeff=1; % Coefficient for the learning rate of the scaling factors
opt.Tscoeff=1; % Coefficient for the learning rate of the global step size
opt.ucr_red=.05; % reduction rate of the critical displacement constraint
opt.feastol=1e-7; % tolerance for constraint violation

opt.ExploitSymmetry= 1;
% ********************* Enter problem specific data here % **********************
D=3; %planar or spatial truss
N_node=65; % number of nodes in the ground structure
N_member=224; % number of nodes in the ground structure
GSCP=[49	34	35	50	34	35	36	50	35	36	37	51	36	37	38	52	37	38	39	53	38	39	40	54	39	40	41	55	40	41	42	56	41	42	43	57	42	43	44	58	43	44	45	59	44	45	46	60	45	46	47	61	46	47	48	62	47	48	49	63	48	49	34	64	33	18	19	34	18	19	20	34	19	20	21	35	20	21	22	36	21	22	23	37	22	23	24	38	23	24	25	39	24	25	26	40	25	26	27	41	26	27	28	42	27	28	29	43	28	29	30	44	29	30	31	45	30	31	32	46	31	32	33	47	32	33	18	48	17	2	3	18	2	3	4	18	3	4	5	19	4	5	6	20	5	6	7	21	6	7	8	22	7	8	9	23	8	9	10	24	9	10	11	25	10	11	12	26	11	12	13	27	12	13	14	28	13	14	15	29	14	15	16	30	15	16	17	31	16	17	2	32	1	2	1	2	1	3	1	4	1	5	1	6	1	7	1	8	1	9	1	10	1	11	1	12	1	13	1	14	1	15	1	16
50	50	50	65	51	51	51	51	52	52	52	52	53	53	53	53	54	54	54	54	55	55	55	55	56	56	56	56	57	57	57	57	58	58	58	58	59	59	59	59	60	60	60	60	61	61	61	61	62	62	62	62	63	63	63	63	64	64	64	64	65	65	65	65	34	34	34	49	35	35	35	35	36	36	36	36	37	37	37	37	38	38	38	38	39	39	39	39	40	40	40	40	41	41	41	41	42	42	42	42	43	43	43	43	44	44	44	44	45	45	45	45	46	46	46	46	47	47	47	47	48	48	48	48	49	49	49	49	18	18	18	33	19	19	19	19	20	20	20	20	21	21	21	21	22	22	22	22	23	23	23	23	24	24	24	24	25	25	25	25	26	26	26	26	27	27	27	27	28	28	28	28	29	29	29	29	30	30	30	30	31	31	31	31	32	32	32	32	33	33	33	33	2	17	3	3	4	4	5	5	6	6	7	7	8	8	9	9	10	10	11	11	12	12	13	13	14	14	15	15	16	16	17	17
]'; % a matrix determining the connectivity plot in the ground structure (N_member rows, 2 columns)
%define grouped members in a matrix
sym_M=[
193	197	195	196	198	130	138	134	136	140	133	131	135	137	66	74	70	72	76	69	67	71	73	2	10	6	8	12	5	3	7	9
201	205	199	194	200	146	154	142	132	144	143	129	141	139	82	90	78	68	80	79	65	77	75	18	26	14	4	16	15	1	13	11
209	213	203	202	206	162	170	150	148	156	149	145	151	153	98	106	86	84	92	85	81	87	89	34	42	22	20	28	21	17	23	25
217	221	207	204	208	178	186	158	152	160	159	147	157	155	114	122	94	88	96	95	83	93	91	50	58	30	24	32	31	19	29	27
217	221	211	210	214	178	186	166	164	172	165	161	167	169	114	122	102	100	108	101	97	103	105	50	58	38	36	44	37	33	39	41
217	221	215	212	216	178	186	174	168	176	175	163	173	171	114	122	110	104	112	111	99	109	107	50	58	46	40	48	47	35	45	43
217	221	219	218	222	178	186	182	180	188	181	177	183	185	114	122	118	116	124	117	113	119	121	50	58	54	52	60	53	49	55	57
217	221	223	220	224	178	186	190	184	192	191	179	189	187	114	122	126	120	128	127	115	125	123	50	58	62	56	64	63	51	61	59
]; % topologically grouped members. 
sym_A=sym_M; % members with similar sections (In most test problems in the literature, sym_M=sym_A)
X_indep_ind=[4     6     7     8    11    52    54    55    56    59   100   102   103   104   107   148   151   152]; % coordinate variables that are not coupled to other coordinates
basic_node=[1 52:4:64]; %Basic nodes that must be present in the design
basic_member=[]; % members that must be active in design
X_const_ind=[ 1 2 3 52*3-[2 1 0]   150 153     (4:16:52)*3-2]; % indices of the fixed coordinates
X_const_val=[ 0 0 10  0 -10 0      0   0   zeros(1,4)]; % values of the fixed coordinates
DOF_constrained=[(50:65)*3-2 (50:65)*3-1 (50:65)*3]; % DOFs anchored by the supports
% Define external load(s) applied to the structure
Fext=zeros(1,3*N_node)';
Fext([1 2 3])=[500 500 -1000]*1000; % D*N_node
%Define Mechanical Properties
M_elasticity=2e11;  %Modulus of elasticity
Fy=248.21e6; % Yield strength according to AISC-ASD 9th edition
Density=7850; % Density of truss material
dis_all=.01; % allowable displacement along each DOF 
% Now specify the list number for the given cross section. Store the cross sections in m-file sections.m. Use a number greater than 100 when AISC-ASD spcifications are employed.
section_no=101; % The available section list for this problem. 
% Specify the search range of shape variables. Only those corresponding the independent shape variables are important
D_X=[0	0	10	-3.75	-3.75	5	-2.5	-3.75	5	-1.25	-3.75	5	0	-3.75	5	1.25	-3.75	5	1.25	-2.5	5	1.25	-1.25	5	1.25	0	5	1.25	1.25	5	0	1.25	5	-1.25	1.25	5	-2.5	1.25	5	-3.75	1.25	5	-3.75	0	5	-3.75	-1.25	5	-3.75	-2.5	5	-7.5	-7.5	2.5	-5	-7.5	2.5	-2.5	-7.5	2.5	0	-7.5	2.5	2.5	-7.5	2.5	2.5	-5	2.5	2.5	-2.5	2.5	2.5	0	2.5	2.5	2.5	2.5	0	2.5	2.5	-2.5	2.5	2.5	-5	2.5	2.5	-7.5	2.5	2.5	-7.5	0	2.5	-7.5	-2.5	2.5	-7.5	-5	2.5	-11.25	-11.25	0	-7.5	-11.25	0	-3.75	-11.25	0	0	-11.25	0	3.75	-11.25	0	3.75	-7.5	0	3.75	-3.75	0	3.75	0	0	3.75	3.75	0	0	3.75	0	-3.75	3.75	0	-7.5	3.75	0	-11.25	3.75	0	-11.25	0	0	-11.25	-3.75	0	-11.25	-7.5	0	-15	-15	0	-10	-15	0	-5	-15	0	0	-15	0	5	-15	0	5	-10	0	5	-5	0	5	0	0	5	5	0	0	5	0	-5	5	0	-10	5	0	-15	5	0	-15	0	0	-15	-5	0	-15	-10	0];
U_X=[0	0	10	-1.25	-1.25	10	0	-1.25	10	1.25	-1.25	10	2.5	-1.25	10	3.75	-1.25	10	3.75	0	10	3.75	1.25	10	3.75	2.5	10	3.75	3.75	10	2.5	3.75	10	1.25	3.75	10	0	3.75	10	-1.25	3.75	10	-1.25	2.5	10	-1.25	1.25	10	-1.25	0	10	-2.5	-2.5	7.5	0	-2.5	7.5	2.5	-2.5	7.5	5	-2.5	7.5	7.5	-2.5	7.5	7.5	0	7.5	7.5	2.5	7.5	7.5	5	7.5	7.5	7.5	7.5	5	7.5	7.5	2.5	7.5	7.5	0	7.5	7.5	-2.5	7.5	7.5	-2.5	5	7.5	-2.5	2.5	7.5	-2.5	0	7.5	-3.75	-3.75	5	0	-3.75	5	3.75	-3.75	5	7.5	-3.75	5	11.25	-3.75	5	11.25	0	5	11.25	3.75	5	11.25	7.5	5	11.25	11.25	5	7.5	11.25	5	3.75	11.25	5	0	11.25	5	-3.75	11.25	5	-3.75	7.5	5	-3.75	3.75	5	-3.75	0	5	-5	-5	0	0	-5	0	5	-5	0	10	-5	0	15	-5	0	15	0	0	15	5	0	15	10	0	15	15	0	10	15	0	5	15	0	0	15	0	-5	15	0	-5	10	0	-5	5	0	-5	0	0];

% *********** No more modifications are required to be applied from this line *********************


availsec=sections(section_no);

Min_A=min(availsec(:,1));
Max_A=max(availsec(:,1));
N_loadcase=size(Fext,2);
is_constrained_DOF=zeros(1,D*N_node);
is_constrained_DOF(DOF_constrained)=1;
DOFonNAP=sum(reshape(is_constrained_DOF,D,N_node)); % used for checking necessary condition of stability
if opt.topology_efficiency_check==1 
     F0=max(abs(Fext)>0,[],2);
     DOFonNAP=DOFonNAP+max(reshape(F0,D,N_node))-1;
end
%Initial valuess of design variables
Xmean=(D_X+U_X)/2;
Xmean(X_const_ind)=X_const_val;
Mmean=opt.initial_M*ones(1,N_member);
Amean=(Max_A+Min_A)/2*ones(1,N_member);
% Find independent variables. Consider only variables that are independent
Xvar_ind=setdiff(X_indep_ind, X_const_ind); 
Avar_ind=setdiff(1:N_member,reshape(sym_A(2:end,:),1,numel(sym_A(2:end,:))));
Avar_ind=union(Avar_ind,sym_A(1,:)); 
Mvar_ind=setdiff(1:N_member,[basic_member,reshape(sym_A,1,numel(sym_A))]);
Mvar_ind=union(Mvar_ind,setdiff(sym_A(1,:),basic_member)); % independent topology variables
N_Xvar=numel(Xvar_ind); % Number of independent shape variables
N_Avar=numel(Avar_ind); % Number of independent size variables
N_Mvar=numel(Mvar_ind);  % Number of independent topology variables
YMEAN=([ Xmean(Xvar_ind) Amean(Avar_ind) Mmean(Mvar_ind)]); %Y is the set of all independent variables
U_Y=[U_X(Xvar_ind) Max_A*ones(1,N_Avar) ones(1,N_Mvar)]; 
D_Y=[D_X(Xvar_ind) Min_A*ones(1,N_Avar) zeros(1,N_Mvar)];
STR=(U_Y-D_Y); % This is vector of scaling factors
Max_Kcond=1e12; % The limit that specifies whether the stiffness matrix is singular
SMEAN=.25; % global step size
% Calculating effective number of design parameters
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
Min_Max_Avar=max(availsec(2:end,1)./availsec(1:end-1,1)); % minimum value of the move limit ratio
Max_Avar=Max_Max_Avar; % move limit ratio
minPcoeff=1;  % minimum value of penalty coefficients
infeas_ratio=opt.max_infeas*ones(1,N_member); % maximum fraction of members that may violate a constraint
counteval=0; % evaluation count
while iter<maxiter % main optimization loop starts here
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
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1); % deteremine absence/presence of topologically coupled members
            [NAP,NAPconnector]=find_NAP(GSCP,M); % find which nodes are active, and the number of the members connected to them
            NAPconnector_all=NAPconnector+DOFonNAP; % The number of members connected to each node + the number of reactions on the nodes
            chk1=(min(NAP(basic_node)==1)>.5); % necessary condition: Are all basic nodes active? 
            chk2=(sum(M) >= (sum(NAP)*D-numel(DOF_constrained))); % Necessary condition: Does the candidate truss have the minimum required number of members?
            chk3= min( max([ NAPconnector_all>=D ; NAP==0])      ); % Necessary condition: Are sufficient number of members connected to each node? 
            if (chk1 && chk2 && chk3)
                stable=1; % The candidate design is accepted if it satisfies all the necessary conditions for stability
            end
        end
        %%% UPPER LEVEL %%%%
            M = zeros(1,N_member);
            M(sym_M(1,:)) = M_upper;
            M(sym_M(2:end,:))= repmat(M(sym_M(1,:)),size(sym_M,1)-1,1); 
            [NAP,NAPconnector]=find_NAP(GSCP,M); % find which nodes are active, and the number of the members connected to them
            NAPconnector_all=NAPconnector+DOFonNAP;
            
        X(Xvar_ind)=Y(k,1:N_Xvar); % shape of the design
        X=enforce_shape_sym(X,N_member,opt);% Enforce shape symmetry, if any
        A(Avar_ind)=Y(k,N_Xvar+1:N_Xvar+N_Avar); % size values
        Ri=zeros(1,N_member); % Radii of gyration of sections, presetting
        [A(Avar_ind),Ri(Avar_ind),~]=round_A_W(A(Avar_ind),M(Avar_ind),section_no,0.5*ones(1,numel(Avar_ind))); % Stochastically round the continuous size values to the closest upper/lower value in the available set of sections   
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
        keep_M(k,:)=M; % store the topology of the design
        %  find out which independent shape and size variables are active in the k-th solution (Y(k,:))  
        [~,q1]=ismember(intersect(X_active_ind,Xvar_ind),Xvar_ind) ; % find index of independent active coordinates
        [~,q2]=ismember(intersect(A_active_ind,Avar_ind),Avar_ind);
        Y_active_ind=[ q1 q2+N_Xvar  (1:N_Mvar)+N_Xvar+N_Avar]; % index of active variables in the k-th solution (Y(k,:))
        Y_W(k,:)=Y(k,:)*0;
        Y_W(k,Y_active_ind)=1; % store which variables in the k-th solution were active (required for excluding the effect of the passive variables in recombination)
        %   Now update the the k-th solution
        Y(k,1:N_Xvar)=X(Xvar_ind);
        Y(k,N_Xvar+1:N_Xvar+N_Avar)=A(Avar_ind);
        Y(k,N_Xvar+N_Avar+1:end)=M(Mvar_ind);
        Z(k,:)=(Y(k,:)-YMEAN)/S(k); % perturbation vector 
        %  analyze the design 
        [Length_M,Vol,displ_ratio,slender_ratio,stress_ratio,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_ASD(NAP,X,DOF_constrained,A,Ri,M,GSCP,Fext,M_elasticity,Fy,dis_all,Max_Kcond,D); 
        counteval=counteval+1; 
        %  assign the highest constraint violation to all the coupled variables
        slender_ratio(sym_A)= repmat(max(slender_ratio(sym_A)),size(sym_A,1),1);
        stress_ratio(sym_A)= repmat(max(stress_ratio(sym_A)),size(sym_A,1),1);
        %  Estimate the required increase in the cross sections such that all constraints are satisfied
        Agoal0=A;
        Agoal2=find_Agoal_ASD_224(Length_M,f_int_ext,M,A,Ri,Fy,M_elasticity,Pcoeff, Max_Avar, section_no,sym_M, Avar_ind); % individual section area increase for satisfaction of member-based constraints 
        
        Agoal3=Agoal0.*max(displ_ratio); % proportional increase for satisfaction of displacement constriants
        Agoal= max([Agoal0+(Agoal2-Agoal0).*Pcoeff;Agoal3]); % The estimated required increase + the initial cross section area
        f(k)=Density*(  Vol+sum((Agoal-A).*M.*Length_M)   ); % the value of the objective function (penalized weight)
        constraint_vio(k,:)=M.*max(0,max([stress_ratio;slender_ratio;max(displ_ratio)*ones(1,N_member)])-1); % contraint violation 
        if (max(constraint_vio(k,:))<=opt.feastol)  && (f(k)<best_feas_design(end)) % update the best feasible solution found so far
            best_feas_design=[X A M f(k)];
        end
        % Now resize. Pre-allocate first
        Y(k+lambda,:)=Y(k,:);  
        A_resized=A;  
        S(k+lambda)=S(k); 
        Y_W(k+lambda,:)=Y_W(k,:);  
        keep_M(k+lambda,:)=M;
        if max(abs(stress_ratio))==0 % the k-th design was kinematically unstable, no resizing is possible
            f(k+lambda)=f(k);
            Z(k+lambda,:)= Z(k,:);
            constraint_vio(k+lambda,:)= constraint_vio(k,:);
        else % the k-th design was kinematically stable, perform resizing
            [A_resized,Ri_resized]=resize_FSDII_ASD(Length_M,dis_all,f_int_ext,displacement,M,A,Ri,Fy,M_elasticity,Pcoeff,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,opt.ucr_red);  % returns the resized section
            Y(k+lambda,N_Xvar+1:N_Xvar+N_Avar)=A_resized(Avar_ind); % Update the resized solution
            Z(k+lambda,:)=(Y(k+lambda,:)-YMEAN)/S(k+lambda);   % Update the corresponding perturbation vector 
            % now analyze the resized solution
            [Length_M,Vol,displ_ratio,slender_ratio,stress_ratio,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_ASD(NAP,X,DOF_constrained,A_resized,Ri_resized,M,GSCP,Fext,M_elasticity,Fy,dis_all,Max_Kcond,D); 
            counteval=counteval+1;
            % assign the most critical constraint to all coupled sections    
            slender_ratio(sym_A)= repmat(max(slender_ratio(sym_A)),size(sym_A,1),1);
            stress_ratio(sym_A)= repmat(max(stress_ratio(sym_A)),size(sym_A,1),1);
            % estimate the required increase in cross sections such that all constraints are satisfied
            Agoal0=A_resized;
            Agoal2=find_Agoal_ASD_224(Length_M,f_int_ext,M,A_resized,Ri_resized,Fy,M_elasticity,Pcoeff, Max_Avar, section_no,sym_M, Avar_ind); % individual increase in cross sections for satisfaction of member-based constraints  
            Agoal3=Agoal0.*max(displ_ratio); % proportional increase for satisfaction of displacement constriants
            Agoal= max([Agoal0+(Agoal2-Agoal0).*Pcoeff;Agoal3]); % the estimated required increase + the initial cross section area
            f(k+lambda)=Density*(  Vol+sum((Agoal-A_resized).*M.*Length_M)   ); % the value of the objective function (penalized weight)
            constraint_vio(k+lambda,:)=M.*max(0,max([stress_ratio;slender_ratio;max(displ_ratio)*ones(1,N_member)])-1); % constraint violation of the resized solution
            if (max(constraint_vio(k+lambda,:))<=opt.feastol)  && (f(k+lambda)<best_feas_design(end)) % update the best identified feasible solution
                best_feas_design=[X A_resized M f(k+lambda)];
            end
        end % resizing completed
    end %  lambda+lambda solutions were created and evaluated
    % Now perform recombination and update strategy parameters
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
end %run completed
bestX=best_feas_design(1:D*N_node);2
bestA=best_feas_design(D*N_node+(1:N_member));
bestM=best_feas_design(D*N_node+N_member+(1:N_member)); % 169:278
bestNAP=find_NAP(GSCP,bestM);
final_eval = counteval;
total_eval = counteval;



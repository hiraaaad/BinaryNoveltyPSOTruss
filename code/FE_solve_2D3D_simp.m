% revised 11/12/2015 by Ali Ahrari
% This function performs FE analysis and returns the output required for calculating weight and constraints wrt simplified specifications.  
function [L_member,vol, Def_ratio_all, buck_ratio_all,  Sig_ratio_all,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_simp(NAP,coordinates,DOF_constrained,A,AAP,GSCP,Fext,M_elasticity,SigT_all,SigC_all,Dis_all,Max_Kcond,D)
%Max_Link: maximum number of members = the number of members in the ground structre
%L_member (1*Max_Link): A vector returning the length of the members
% Def_ratio (1*2N) or (1*3N): returns the absolute value of the ratio of deflection along a DOF  to the Aloable limit
% buck_ratio (1*max_link): returns the ratio of compresive force to the buckling limit (for tensile forces this value is zero)
% Sig_ratio (1*max_link): returns the ratio of member stress to the allowable stress (absolute value)
% NAP (1*N) : consists of 0/1 elements which shows whether a nod is active (1) or passive (0)
% nodes (1*2N) or (1*3N): stores coordinates of the nodes . (x1 y1 z1 x2 y2 z2 ... zn yn zn)
% DOF_constrained: The indii of degrees of freedom which are constrained.
% A (1*Max_Link): The member areas of all elements
% AAP (1*Max_Link) : specifies whether a member is active or passive
% GSCP (Max_Link*2) : the matrix that specifies the two nodes of each member
% Fext (2N*1) or (3N*1) : the vector of external load on nodes
% M_elasticity: E
%SigT_all,SigC_all: allowable stress limit (tensile and compressive)
%Def_all allowable deflection
%Max_Kcond the maximum condition number of global stiffness matrix
% GSCP
% A
% AAP
%pause
X=coordinates(1:D:end);
Y=coordinates(2:D:end);
N_loadcase=size(Fext,2);
N_member=numel(AAP);
N_node=numel(NAP);
if D==3
    Z=coordinates(3:D:end);
end
I=A.^2; %assume this and use Euler criteria. if I=kA^2, use k in expected safety factor for buckling
% for 160-bar truss we use specific buckling criteria
r_buckling = [0.47, 0.57, 0.67, 0.77, 0.87, 0.97, 0.97, 1.06, 1.16, 1.26,...
    1.15, 1.26, 1.36, 1.46, 1.35, 1.36, 1.45, 1.55, 1.75, 1.95, 1.74, 1.94,...
    2.16, 2.36, 2.57, 2.35, 2.56, 2.14, 2.33, 2.97, 2.54, 2.93, 2.94, 2.94,...
    2.92, 3.54, 3.96, 3.52, 3.51, 3.93, 3.92, 3.92]*1e-2;
%
L_member=zeros(1,N_member);
Kg=zeros(N_node*D,N_node*D); % initial values of global stiffness matrix
if D==3
    for r=1:N_member
        if AAP(r)==1 % if the member exist, apply its contribution to the global stiffness matrix
            inds=GSCP(r,:);
            m=inds(1); %1st node of member
       	    n=inds(2); %2nd node of  member
            L_member(r)=norm([X(n)-X(m) Y(n)-Y(m) Z(n)-Z(m)]);
      	    CTX=(X(n)-X(m))/L_member(r);
            CTY=(Y(n)-Y(m))/L_member(r);
            CTZ=(Z(n)-Z(m))/L_member(r);
            kk=M_elasticity*A(r)/L_member(r)*[1 -1;-1 1];
            L_mat=[CTX CTY CTZ 0 0 0;0 0 0 CTX CTY CTZ];    
            k=L_mat'*kk*L_mat;   
            indii=[3*m-2 3*m-1 3*m 3*n-2 3*n-1 3*n];
            for g1=1:6
                for g2=1:6
                     Kg(indii(g1),indii(g2))= Kg(indii(g1),indii(g2))+k(g1,g2);
                end
            end
        else
            L_mat=zeros(2,2*D);
        end
        keep_L_mat(2*r-1:2*r,:)=L_mat; %matrices converting nodal coordinates to global coordinates
    end
elseif D==2
    for r=1:N_member
        if AAP(r)==1
            inds=GSCP(r,:);
            m=inds(1);
            n=inds(2);
            L_member(r)=norm([X(m)-X(n) Y(m)-Y(n)]);
            teta=atan( (Y(n)-Y(m)) / (X(n)-X(m))  );
            if (X(n)-X(m))<0
                teta=teta+pi;
            end
            L_mat=[cos(teta) sin(teta) 0 0;0 0 cos(teta) sin(teta)];
            k=L_mat'*M_elasticity*A(r).*AAP(r)/L_member(r)*[1 -1;-1 1]*L_mat;
            indii=[2*m-1 2*m 2*n-1 2*n];
            for g1=1:4
                for g2=1:4
                    Kg(indii(g1),indii(g2))= Kg(indii(g1),indii(g2))+k(g1,g2);
                end
            end
        else
            L_mat=zeros(2,2*D);
        end
        keep_L_mat(2*r-1:2*r,:)=L_mat; %matrices converting nodal coordinates to global coordinates
    end
end
omit_DOF=zeros(1,D*N_node);
omit_DOF(1:numel(DOF_constrained))=DOF_constrained;
omit_index=numel(DOF_constrained);
for kkk=1:N_node
    if NAP(kkk)==0
        omit_DOF(omit_index+(1:D))= D*kkk-((D-1):-1:0);
        omit_index=omit_index+D;
    end
end
keep_DOF=setdiff(1:D*numel(NAP),omit_DOF);
Kg_red2=Kg(keep_DOF,:);
Kg_red=Kg_red2(:,keep_DOF);
Fext_red=Fext(keep_DOF,:);
stable=(rcond(Kg_red)>(1/Max_Kcond));
if ~ stable
    Def_ratio_all=zeros(1,D*numel(NAP));
    buck_ratio_all=zeros(1,N_member);
    Sig_ratio_all=buck_ratio_all;
    displacement=Fext*0;
    f_int_ext=zeros(N_member,size(Fext,2));
    req_AAP=D*sum(NAP)-numel(DOF_constrained);
    nulity=numel(keep_DOF)-rank(Kg_red);
    %penalt=nulity+(sum(AAP)-req_AAP)/req_AAP; % this one is better in 
    %directing towards stable topologies (224 bar with no symmetry, for
    %example)
    penalt=(nulity+sum(AAP)-req_AAP)/req_AAP;
    vol=1e100*(1+penalt);
    f_int_unit=[];
else
    vol=(A.*AAP)*L_member';
    U=Fext*0; % initial value of node displacement when the external force is applied 
    buck_ratio=zeros(N_member,N_loadcase);
    Sig_ratio=zeros(N_member,N_loadcase);
    Fext_unit=eye(numel(keep_DOF)); %unit load matrix
    if Dis_all<1e10 %There are disoplcement constriant, find the unit load response
        U_unit=zeros(numel(coordinates));%initial value of node displacement when the unit force is applied 
        allU=Kg_red\[Fext_unit Fext_red]; %Dispalcememnts of all nodes for external unit & real loads
        U_unit(keep_DOF,keep_DOF)=allU(:,1:numel(keep_DOF)); % Nodal displacement when the unit loads are applied
        U(keep_DOF,1:N_loadcase)=allU(:,(1:N_loadcase)+numel(keep_DOF)); % Nodal displacement when the actual loads are applied 
        f_int_unit=zeros(N_member,N_node*D);
        for r=1:N_member   
            if AAP(r)==1
                inds=GSCP(r,:);
                m=inds(1);
                n=inds(2);
                indii=[D*m+1-(D:-1:1) D*n+1-(D:-1:1) ];
                L_mat=keep_L_mat(2*r-1:2*r,:);
                Ttool=L_mat*U_unit(indii,:); 
                Sig_unit=(Ttool(2,:)-Ttool(1,:))*M_elasticity/L_member(r);  
                f_int_unit(r,:)=Sig_unit*A(r);
            end
        end
        f_int_unit=f_int_unit'; %axiaql forces when unit loads are applied
    else %no displacmeent constriant, do not calculate unit load response
        U(keep_DOF,1:N_loadcase)=Kg_red\Fext_red;
        f_int_unit=[];
    end
    displacement=U;
    f_int_ext=zeros(N_member,size(Fext,2));
    for r=1:N_member
        if AAP(r)==1
            inds=GSCP(r,:);
            m=inds(1);
            n=inds(2);
            indii=[D*m+1-(D:-1:1) D*n+1-(D:-1:1) ];
            L_mat=keep_L_mat(2*r-1:2*r,:);
            UE=L_mat*U(indii,:);
            Sig=(UE(2,:)-UE(1,:))*M_elasticity/L_member(r);  
            f_int_ext(r,:)=Sig*A(r);
            Sig_ratio(r,:)=(Sig/SigT_all).*(Sig>0) +  abs(Sig/SigC_all).*(Sig<0);
            if N_member == 180
                temp = L_member(r)/r_buckling(sections(60) == A(r));
                if temp <= 120
                    sigma_b = 1300 - (temp^2)/24;
                else
                    sigma_b = 10^7 / (temp^2);
                end
                buck_ratio(r,:)=abs(f_int_ext(r,:)/sigma_b)   .* (Sig<0);
            else
                buck_ratio(r,:)=abs(f_int_ext(r,:)/(pi^2*I(r)/L_member(r)^2))   .* (Sig<0);
            end
        end
    end
    % return the critical case
Def_ratio_all=(max(abs(U),[],2)/Dis_all)';
buck_ratio(buck_ratio>0)=0; % no buck ratio in my problems
buck_ratio_all=max(buck_ratio,[],2)';
Sig_ratio_all=max(Sig_ratio,[],2)';
end


 


 
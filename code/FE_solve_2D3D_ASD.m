%This function performs finite element analysis of the truss and returns
%the output data including the response to the unit load loads, when AISC-ASD design specifications govern. 
function [L_member,vol, Def_ratio_all, Slend_ratio_all,  Sig_ratio_all,f_int_ext,displacement,f_int_unit]=FE_solve_2D3D_ASD(NAP,coordinates,DOF_constrained,A,Ri,AAP,GSCP,Fext,M_elasticity,Fy,Dis_all,Max_Kcond,D)
X=coordinates(1:D:end);
Y=coordinates(2:D:end);
N_loadcase=size(Fext,2);
N_member=numel(AAP);
N_node=numel(NAP);
if D==3
    Z=coordinates(3:D:end);
end
L_member=A*0;%initial value of member lengths
Kg=zeros(N_node*D,N_node*D); % initial values of global stiffness matrix
if D==3 %spatial truss analysis
    for r=1:N_member
        if AAP(r)==1 % if the member exist, apply its contribution to the global stiffness matrix
            inds=GSCP(r,:);
            m=inds(1);
       	    n=inds(2);
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
        keep_L_mat(2*r-1:2*r,:)=L_mat; % matrices converting nodal coordinates to global coordinates
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
keep_DOF=setdiff(1:D*N_node,omit_DOF);
Kg_red2=Kg(keep_DOF,:);
Kg_red=Kg_red2(:,keep_DOF);
Fext_red=Fext(keep_DOF,:);
stable=(rcond(Kg_red)>(1/Max_Kcond));
if ~stable
    Def_ratio_all=zeros(1,D*N_node);
    Slend_ratio_all=zeros(1,numel(AAP));
    Sig_ratio_all=Slend_ratio_all;
    displacement=Fext*0;
    f_int_ext=zeros(N_member,N_loadcase);
    req_AAP=D*sum(NAP)-numel(DOF_constrained);
    nulity=numel(keep_DOF)-rank(Kg_red);
    %penalt=nulity+(sum(AAP)-req_AAP)/req_AAP; % this one is better in 
    %directing towards stable topologies (224 bar with no symmetry, for
    %example)
    penalt=(nulity+sum(AAP)-req_AAP)/req_AAP;
    penalt=nulity+(sum(AAP)-req_AAP)/req_AAP;
    vol=1e100*(1+penalt);
    f_int_unit=[];
else
    vol=(A.*AAP)*L_member';
    U=Fext*0; % initial value of node displacement when the external force is applied 
    cc=sqrt(2*pi^2*M_elasticity/Fy);
    lambdam=L_member./Ri;
    numerator=(1-lambdam.^2/2/cc^2)*Fy;
    denominator=5/3+3/8/cc*lambdam-lambdam.^3/8/cc^3;
    slender_ratio=zeros(N_member,size(Fext,2));
    Sig_ratio=zeros(N_member,size(Fext,2));
    Fext_unit=eye(numel(keep_DOF)); %unit load matrix
    if Dis_all<1e10 %There are disoplcement constriant, find the unit load response
        U_unit=zeros(numel(coordinates));%initial value of node displacement when the unit force is applied 
        allU=Kg_red\[Fext_unit Fext_red]; %Dispalcememnts of all nodes for external unit & real loads
        U_unit(keep_DOF,keep_DOF)=allU(:,1:numel(keep_DOF)); % Nodal displacement when the unit loads are applied
        U(keep_DOF,1:N_loadcase)=allU(:,(1:N_loadcase)+numel(keep_DOF)); % Nodal displacement when the actual loads are applied 
        % Kg_red_inv=inv(Kg_red);
        % U_unit(keep_DOF,keep_DOF)=Kg_red_inv;
        % U(keep_DOF,1:N_loadcase)=Kg_red_inv*Fext_red;
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
        f_int_unit=f_int_unit'; % axiaql forces when unit loads are applied
    else %no displacement constriant, do not calculate unit load response
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
            slender_ratio(r,:)=lambdam(r)/300.*(f_int_ext(r,:)>=0) +  lambdam(r)/200.*(f_int_ext(r,:)<0);
            if lambdam(r)<cc % inelastic buckling
                 tanesh=abs(Sig/(numerator(r)/denominator(r)));
            else  % M_elasticity buckling
                 tanesh=abs(Sig/(12*pi^2*M_elasticity/23/lambdam(r)^2));
            end
            Sig_ratio(r,:)=(Sig/Fy/.6).*(Sig>=0) +  tanesh.*(Sig<0);
        end
    end
Def_ratio_all=(max(abs(U),[],2)/Dis_all)';
Slend_ratio_all=max(slender_ratio,[],2)';
Sig_ratio_all=max(Sig_ratio,[],2)';
end




 
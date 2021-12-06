% 3/19/2016: Critical error in computation of the critical member force (exclusing coupled
% member when subject to one load case) was corrected
function [A,Ri]=resize_FSDII_ASD(L_member,dis_all,f_int_ext,displacement,AAP,A,Ri,Fy,M_elasticity,Pzarib,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,c_sec_red)
% L_link: Length of the members
%Vol: Volume of truss
%Def_ratio: absolute value of the ratio of deflection to the allowable deflection
%dis_all: Allowable deflection
%buck_ratio: Buckling ratio, square root
%stress_ratio: ratio fo the maximum stress to the allowable stress (absolute)
%Kg_red: inverse of the reduceds stiffness matrix of the truss
%f_int_ext: force of members
%Deflection: Deflection of DOFs
%keep_DOF: DOF that are not constrained
%keep_L_mat: Convert global coordinates to element coordinates, for all members
%AAP,A,TA,elastic,Max_A,Min_A: Like before
%Pzarib
%Max_Avar: maximum variation (ratio) in section area
%bias_round: bias in rounding the sections towards the upoper value
%section_No: section list
%sym_AAP: symmetric members
%zaribpenalt
initialA=A;
time0=toc;
if   max(var(A(sym_A).*AAP(sym_A)))>.000001 % Checkpoint for coupled member sections
   disp('Error: The coupled member sections found to be unequal. Aborting the optimization')
   A(sym_A).*AAP(sym_A)
   var(A(sym_A).*AAP(sym_A))
   Errorr        
end     
% Finding independent sections
indep_members=setdiff(Avar_ind.*AAP(Avar_ind),0);
dep_members=setdiff(1:numel(AAP),indep_members);
maxA=A*0;
maxA(indep_members)=round_A_W(initialA(indep_members)*Max_Avar,AAP(indep_members),section_no,ones(1,numel(indep_members))); %The upper limit of A (vector) wrt move limit
FSDred=exp(1-1*Pzarib); %reliance on FSD for reducing a crosss section
%avail_sec=sections(section_no);
%if size(f_int_ext,2)==1
%   f_int_ext=[f_int_ext f_int_ext];
%end
%if size(f_int_ext,2)>1 % The structure is subject to multiple load cases. Find the maximum and minimum member force for member-based resizing 
    max_f_int_ext=max(f_int_ext',[],1);
    min_f_int_ext=min(f_int_ext',[],1);
    %Assign the most critical axial force to all coupled members
    max_f_int_ext(sym_A(1:end,:))=repmat(max(max_f_int_ext(sym_A(1:end,:))),size(sym_A,1),1); 
    min_f_int_ext(sym_A(1:end,:))=repmat(min(min_f_int_ext(sym_A(1:end,:))),size(sym_A,1),1);
    critical_f_int_ext=[max_f_int_ext;min_f_int_ext]';
%else %The structure is subject to a single load case
 %   critical_f_int_ext=f_int_ext;
%end
activeAvar_ind= intersect(Avar_ind, find(AAP)); % active independent members
for r=activeAvar_ind %now resize active independent members one after another
    Max_Avar_dec=(Max_Avar-1)*FSDred(r)+1; 
    [A(r),Ri(r),~,~]=find_proper_section_ASD(section_no,A(r),L_member(r),critical_f_int_ext(r,:)',Fy,M_elasticity,Max_Avar,Max_Avar_dec);
end
%Assign the value of dependent members
A(sym_A(2:end,:))=repmat(A(sym_A(1,:)),size(sym_A,1)-1,1);
Ri(sym_A(2:end,:))=repmat(Ri(sym_A(1,:)),size(sym_A,1)-1,1);



if dis_all<1e10 % Now resize for deflections
    active_member=find(AAP);
    %Def=fe_kol(:,active_member)*    (    f_int_ext(active_member,:)  .* repmat( (L_member(active_member).*(1./A(active_member)))',1,size(f_int_ext,2))     )/M_elasticity;
    deltaU=f_int_unit(:,active_member)*    (    f_int_ext(active_member,:)  .* repmat( (L_member(active_member).*(1./A(active_member)-1./initialA(active_member)))',1,size(f_int_ext,2))     )/M_elasticity;
    displacement=deltaU+displacement;
    % calculating average CE of coupled members
    W_sym_AAP=ones(1,size(sym_A,2));
    W_sym_AAP=[W_sym_AAP;(abs(diff(sym_A))>0)]; % This is used to remove repetion of a member in sym_AAP matrix
    iter=0;
    while max(max(abs(displacement)))>dis_all & iter<40000 %gradually increase section of the most effective member till all displacement-constraints are satisfied
        iter=iter+1;
        ind=find_largest_element_ind(abs (displacement));% most critical [DOF, load case];
        Ucr=displacement(ind(1),ind(2));
        SEcoeff=(f_int_unit(ind(1),:)'.*f_int_ext(:,ind(2))/M_elasticity)'; % coefficient of save effectivenss, (SEcoeff=fF/(E))
        SEcoeff(sym_A(1,:))=sum(W_sym_AAP.*SEcoeff(sym_A).*L_member(sym_A) ) ./ sum(W_sym_AAP.*L_member(sym_A)+1e-100); %The average value of cost effectiveness is asigned for coupled members
        SEcoeff(dep_members)=0; %never select a dependent member for resizing
        L_link_total=L_member;
        L_link_total(sym_A(1,:))=sum(L_member(sym_A(:,:)).*W_sym_AAP);
        goalUcr=sign(Ucr)*   max(dis_all,abs((1-c_sec_red)*Ucr));% taregt displacememnt
        [recA,changeind]=   resize_for_critical_displacement2(AAP(indep_members),A(indep_members),maxA(indep_members),SEcoeff(indep_members),Ucr,goalUcr  , L_link_total(indep_members),section_no  );
        recRi=find_Ri(recA,section_no);
        if numel(changeind)==0  %effective areas has reached the limit, cannot reduce U
            break 
        else
            oldA=A;
            A(indep_members)=recA;
            Ri(indep_members)=recRi;
            A(sym_A(2:end,:))=repmat(A(sym_A(1,:)),size(sym_A,1)-1,1);
            Ri(sym_A(2:end,:))=repmat(Ri(sym_A(1,:)),size(sym_A,1)-1,1);
            changed_members=find(~(A==oldA));% all members that weree changed
            delta_U=f_int_unit(:,changed_members)*(f_int_ext(changed_members,:).* repmat( [L_member(changed_members).* (1./A(changed_members)-1./oldA(changed_members))]',1,size(f_int_ext,2)) )/M_elasticity; % change in all displacements because of chenge in the sections
            displacement=delta_U+displacement; %update all displacements
        end
    end
end
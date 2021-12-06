% Revised by Ali Ahrari on 11-18-2015
%This function resizes a member such that constraints are satisfied.
function A=resize_FSDII_simp(L_member,Def_all,buck_ratio,stress_ratio,f_int_ext,displacement,AAP,A,M_elasticity,Pzarib,Max_Avar,section_no,sym_A,Avar_ind,f_int_unit,c_sec_red)
% L_member: Length of the members
%Vol: Volume of truss
%Def_ratio: absolute value of the ratio of deflection to the allowable deflection
%Def_all: Allowable deflection
%buck_ratio: Buckling ratio, square root
%stress_ratio: ratio fo the maximum stress to the allowable stress (absolute)
%Kg_red: inverse of the reduceds stiffness matrix of the truss
%Fint: internal force, which is axial force of members subject to external load cases
%U: Displacement of DOFs
%keep_DOF: DOF that are not constrained
%keep_L_mat: Convert global coordinates to element coordinates, for all members
%AAP,A,TA,elastic,Max_A,Min_A: Like before
%Pzarib
%Max_Avar: maximum variation (ratio) in section area
%bias_round: bias in rounding the sections towards the upoper value
%section_No: section list
%sym_AAP: symmetric members
%f_unit [DxN_n,Max_Link): member axial force when truss is subject to unit load. 
%the k-th line is the axial force in the k-th unit unit load 
%zaribpenalt
%c_sec_red=.05;
initialA=A;
minA=initialA*0;maxA=initialA*0; %lower and upper resizing limit considering move limits
%Ri0=Ri;
if section_no>0
    avail_sections=sections(section_no);
    availA=avail_sections(:,1)';
    %availRi=avail_sections(:,2)';
    N_sections=numel(availA);
end
if   max(var(A(sym_A).*AAP(sym_A)))>.000001 % Checkpoint for coupled member sections 
     disp('Error: The coupled member sections found to be unequal. Aborting the optimization')
     A(sym_A).*AAP(sym_A)
     Errorr
end     %checkpoint completed
% Finding independent sections
indep_members=setdiff(Avar_ind.*(AAP(Avar_ind)),0); %independent and active members
dep_members=setdiff( 1:numel(AAP),indep_members); %dependent or pasive memebrs
FSDred=exp(1-1*Pzarib); % reliance on FSD for reducing a crosss section
Max_Avar_dec=(Max_Avar-1).*FSDred+1; %the ratio (>1) that a section can be decreased
critical=(max([stress_ratio;buck_ratio])); %the critical constriant ratio between stress and buckling-change the name of buck ratio since it includes the alpha effect and squre root
minA(indep_members)=(initialA(indep_members)./Max_Avar_dec(indep_members)); %The lower limit of A (vector) wrt move limit
maxA(indep_members)=(initialA(indep_members)*Max_Avar); %The upper limit of A (vector) wrt move limit
[minA(indep_members),~,~]=round_A_W(minA(indep_members),ones(1,numel(indep_members)),section_no,ones(1,numel(indep_members))); %round it up (only independent active sections)
[maxA(indep_members),~,~]=round_A_W(maxA(indep_members),ones(1,numel(indep_members)),section_no,ones(1,numel(indep_members))); %round it up
A(indep_members)=max(minA(indep_members),initialA(indep_members).*critical(indep_members)); 
A(indep_members)=min(maxA(indep_members),A(indep_members));%  continuous resized sections, move limits applied
[A(indep_members),~]=round_A_W(A(indep_members),AAP(indep_members),section_no,AAP(indep_members));  % discrete values of sections by rounding up the continuous values
A(sym_A(2:end,:))=repmat(A(sym_A(1,:)),size(sym_A,1)-1,1); %the dependent sectionss are updated
%Arank(sym_AAP(2:end,:))=repmat(Arank(sym_AAP(1,:)),size(sym_AAP,1)-1,1); % the rank of sections in the given section list
A4=A;
A=A.*AAP+(1-AAP).*initialA;
if norm(A4-A)>0
    'eshkal dar resize'
    [A4;A]'
    pause
end
iter=0;
if Def_all<1e10 % Now resize for displacement constriants
    activeAAP=find(AAP);
   % Def=f_unit(:,activeAAP)*    (    Fint(activeAAP,:)  .* repmat( (L_member(activeAAP).*(1./A(activeAAP)))',1,size(Fint,2))     )/elastic %displacement after member-based resizing
    deltaU=f_int_unit(:,activeAAP)*    (    f_int_ext(activeAAP,:)  .* repmat( (L_member(activeAAP).*(1./A(activeAAP)-1./initialA(activeAAP)))',1,size(f_int_ext,2))     )/M_elasticity;
    displacement=deltaU+displacement;
    % caculating average CE of coupled members
    W_sym_AAP=ones(1,size(sym_A,2));
    W_sym_AAP=[W_sym_AAP;(abs(diff(sym_A))>0)]; % This is used to remove repetion of a member in sym_AAP matrix
    iter=0;
  %  meghdar=(max(max( abs(U)))/Def_all)^.05;
    while  max(max( abs(displacement)))>Def_all & iter<40000 %gradually increase section of the most effective member till all U-constraints are satisfied
        iter=iter+1;
        %W_sym_AAP=sym_AAP*0+1; %fagaht baraye moghayese ba old vesion
        %  oldA=A;oldDef=U;Def5=U;A5=A;
        iter=iter+1;
        ind=find_largest_element_ind(abs(displacement));% index of the most critical [DOF, load case];
        Ucr=displacement(ind(1),ind(2));
        SEcoeff=(f_int_unit(ind(1),:)'.*f_int_ext(:,ind(2))/M_elasticity)'; % coefficient of save effectivenss, (SEcoeff=fF/(E))
        SEcoeff(sym_A(1,:))=sum(W_sym_AAP.*SEcoeff(sym_A).*L_member(sym_A) ) ./ sum(W_sym_AAP.*L_member(sym_A)+1e-100); %The average value of cost effectiveness is asigned for coupled members
        SEcoeff(dep_members)=0; %never select a dependent member for resizing
        L_link_total=L_member;
        L_link_total(sym_A(1,:))=sum(L_member(sym_A(:,:)).*W_sym_AAP); %total length of coupled members is considered for indep members. 
        goalUcr=sign(Ucr)*   max(Def_all,abs((1-c_sec_red)*Ucr));% taregt displacememnt
        %now increase size wer critical U considering save effectiveness
        [recA,changeind]=   resize_for_critical_displacement(AAP(indep_members),A(indep_members),maxA(indep_members),SEcoeff(indep_members),Ucr,goalUcr  , L_link_total(indep_members),section_no  );
        if numel(changeind)==0  %effective areas has reached the limit, cannot reduce U
            break 
        else
            %recRi=Ri(indep_members);
            % [recA(changeind),recRi(changeind)]=increase_section_ASD(section_no,recA(changeind),minlambdam); %eslah shavad, 200 or 300
            oldA=A;
            A(indep_members)=recA; %update A
            %  Ri(indep_members)=recRi;
            A(sym_A(2:end,:))=repmat(A(sym_A(1,:)),size(sym_A,1)-1,1); %update the coupled members
            %   Ri(sym_AAP(2:end,:))=repmat(Ri(sym_AAP(1,:)),size(sym_AAP,1)-1,1);
            %  keepU=U;
            %  U= f_unit(:,activeAAP)*    (    Fint(activeAAP,:)  .* repmat( (L_member(activeAAP).*(1./A(activeAAP)-0))',1,size(Fint,2))     )/elastic;
            % delta_Def=f_unit(:,activeAAP)*    (    Fint(activeAAP,:)  .* repmat( (L_member(activeAAP).*(1./A(activeAAP)-1./initialA(activeAAP)))',1,size(Fint,2))     )/elastic
            changed_members=find(~(A==oldA));% all members that weree changed
            delta_U=f_int_unit(:,changed_members)*(f_int_ext(changed_members,:).* repmat( [L_member(changed_members).* (1./A(changed_members)-1./oldA(changed_members))]',1,size(f_int_ext,2)) )/M_elasticity; % change in all displacements because of chenge in the sections
            displacement=delta_U+displacement; %update all displacements
        end
    end
end




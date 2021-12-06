function [Agoal]=find_Agoal_ASD_224(L_link, Felement, AAP,A,Ri,Fy, elastic, Pzarib, Max_Avar, section_no,sym_AAP, Avar_ind)

% L_link: Length of the members
%Vol: Volume of truss
%Def_ratio: absolute value of the ratio of deflection to the allowable deflection
%Def_all: Allowable deflection
%buck_ratio: Buckling ratio, square root
%stress_ratio: ratio fo the maximum stress to the allowable stress (absolute)
%Kg_red: inverse of the reduceds stiffness matrix of the truss
%Felement: force of members
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
% Def0=Deflection;
time0=toc;
Max_Link=numel(AAP);
if   max(var(A(sym_AAP).*AAP(sym_AAP)))>.000001
             max(var(A(sym_AAP)))
            'A voroodi be resize sym nist'
             vvvv=A(sym_AAP).*AAP(sym_AAP)
              sym_AAP
            
end     %checkpoint completed

% Finding independent sections
indepA=setdiff(Avar_ind.*(AAP(Avar_ind)==1),0);
depA_ind=setdiff( 1:numel(AAP),indepA);
%indepA_ind=setdiff(1:Max_Link,depA_ind);
A0=A;
Ri0=Ri;
FSDA=exp(1-1*Pzarib);
%avail_sec=sections(section_no);
time05=toc;
keep_Felement=Felement;
if size(Felement,2)==1
   Felement=[Felement Felement];
end

if size(Felement,2)>1
maxcriticFelement=max(Felement');
mincriticFelement=min(Felement');
maxcriticFelement(sym_AAP)=repmat(max(maxcriticFelement(sym_AAP)),size(sym_AAP,1),1);
mincriticFelement(sym_AAP)=repmat(min(mincriticFelement(sym_AAP)),size(sym_AAP,1),1);
criticFelement=[maxcriticFelement;mincriticFelement]';
end

Felement=keep_Felement;
minlambdam=zeros(1,numel(AAP));
ind_activeAvar_ind=find(AAP(Avar_ind)==1);
activeAvar_ind= Avar_ind(ind_activeAvar_ind);

for r=activeAvar_ind
 
        Max_Avar_dec=(Max_Avar-1)*FSDA(r)+1; 
         
        [A(r),Ri(r)]=find_proper_section_eff4(section_no,A(r),Ri(r),L_link(r),criticFelement(r,:)',Fy,elastic,Max_Avar,Max_Avar_dec);
        
       % alaki=find_proper_section(section_no,A0(r),L_link(r),criticFelement(r,:)',Fy,elastic,Max_Avar,Max_Avar_dec);

        if 0%~(A(r)==alaki(1))
               KK=rng;
               KK=double(KK.Seed);
               filename=['error' num2str(KK) '.xls']
               xlswrite(filename,0)

              'eshkal dar resize_eff7'
               [A0(r) A(r) alaki(1) ]
                criticFelement(r,:)'/(.6*Fy)
              pause
        end
       % alaki=find_proper_section(section_no,A0(r),L_link(r),criticFelement(r,:)',Fy,elastic,Max_Avar,Max_Avar_dec);
  
end
%[alaki1,alaki2]=max(A(sym_AAP));
A(sym_AAP(2:end,:))=repmat(A(sym_AAP(1,:)),size(sym_AAP,1)-1,1);


minlambdam(sym_AAP(2:end,:))=repmat(minlambdam(sym_AAP(1,:)),size(sym_AAP,1)-1,1);
Ri(sym_AAP(2:end,:))=repmat(Ri(sym_AAP(1,:)),size(sym_AAP,1)-1,1);

Agoal=A;
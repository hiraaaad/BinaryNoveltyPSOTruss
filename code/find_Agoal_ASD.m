%starting 11-19-2015, incomplete
function A=find_Agoal_ASD(L_member,f_int_ext,AAP,A,Ri,Fy,M_elasticity,section_no,sym_A,Avar_ind,opt)
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

% Finding independent sections
if size(f_int_ext,2)>1 % The structure is subject to multiple load cases. Find the maximum and minimum member force for member-based resizing
 % attention: coupled members must have identical                        % length
    max_f_int_ext=max(f_int_ext');
    min_f_int_ext=min(f_int_ext');
else
    max_f_int_ext=f_int_ext';
    min_f_int_ext=f_int_ext';
end
%Assign the most critical axial force to all coupled members
max_f_int_ext(sym_A(1:end,:))=repmat(max(max_f_int_ext(sym_A(1:end,:))),size(sym_A,1),1); 
min_f_int_ext(sym_A(1:end,:))=repmat(min(min_f_int_ext(sym_A(1:end,:))),size(sym_A,1),1);
critical_f_int_ext=[max_f_int_ext;min_f_int_ext]';
activeAvar_ind= intersect(Avar_ind, find(AAP)); % active independent members
Max_Avar=1e100;
Max_Avar_dec=1+1e-12;
for r=activeAvar_ind %now resize active independent members one after another
    [A(r),Ri(r),~,ratios]=find_proper_section_ASD2(section_no,A(r),L_member(r),critical_f_int_ext(r,:)',Fy,M_elasticity,Max_Avar,Max_Avar_dec);
     ratios=ratios.^[opt.ratiosUnsucc(1) opt.ratiosUnsucc(2)]+opt.ratiosUnsucc(3)*(ratios>1);
    A(r)=A(r)*max([1 ratios]);
end
%Assign the value of dependent members
A(sym_A(2:end,:))=repmat(A(sym_A(1,:)),size(sym_A,1)-1,1);

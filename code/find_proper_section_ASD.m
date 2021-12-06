% revised 11-16-2015 by Ali Ahrari
% This function resize a member for the given axial forces (f_int_ext)
% regarding the stress and slenderness ratio constraints according to
% AISC-ASD specifications (9th-edition)
function [A,Ri,rnk_A,acceptable]=find_proper_section_ASD(section_no,A,L,f_int_ext,Fy,M_elasticity,Max_Avar,Max_Avar_dec)
%section_no (scalar): section number, which determines the available sections, including the radii of gyrations 
%A (scalar): cross section area of the current section
%Ri (scalar): radius of gyration of the current section
%minRi: The minimum required radius of gyration (
%L (scalar): Length of the current member
%f_int_ext (a vector of size two): The highest and lowest axial forces in the member 
%Fy (scalar): Yield strength of the material according to AISC-ASD
%M_elasticity (scalar): Modulus of elasticity of the material
%Max_Avar (scalar): maximum variation during resizing, which controls the move limit, equivalent to parameter c_A in the FSD-ESII paper)
%Max_Avar_dec (scalar): This parameter specifies the maximum reduction of the section area during resizing, which is equal or smaller than Max_Avar;
avail_sec=sections(section_no); %list of available sections (areas and radii of gyration)
minA=(A/Max_Avar_dec); % the lower bound on section area for member-based resizing
maxA=min(A*Max_Avar,max(avail_sec(:,1))); % the upper bound on section area for member-based resizing
[minA,~,rnk_minA]=round_A_W(minA,1,section_no,1); %round_up and find the rank of the minA in the given set
[maxA,~,rnk_maxA]=round_A_W(maxA,1,section_no,1); %round_up and find the rank of the maxA in the given set

% member-based resizing is performed depending on which one is the case:
%i) member undergoes only tension (both numbers in f_int_ext are positive)
%ii)member undergoes only compression (both numbers in f_int_ext are negative)
%iii) member undergoes tension and compression (the first number in f_int_ext is positive and the second one is negative)
if numel(f_int_ext)==2 %more than one load case is applied to the structure
    if f_int_ext(2)>=0 % member is in tension in all load cases
        f_int_ext=f_int_ext(1); % consider the maximum force only
    elseif f_int_ext(1)<=0 %member force is tension in all load cases
        f_int_ext=f_int_ext(2);  %consider the minimum force only (with largest absolute value)
    end %otherwise consider both axial forces in f_int_ext (case III is encountered)
end
% now f_int_ext is a scalar (case I or case II), or a vector of size two (case III)
if numel(f_int_ext)==1 % there is only one axial force
    if f_int_ext>=0 %the axial force is tensile 
        
        A=min(maxA,max(f_int_ext/(.6*Fy),minA));
        acceptable_stress=(f_int_ext/(.6*Fy))<=maxA; %acceptability of the member wrt stress constraint
        [A,Ri,rnk_A]=round_A_W(A,1,section_no,1); % round up the calculated area, find the location of the new area in the given set and the corresponding radius of gyration
        acceptable_slenderness=((L/Ri)<=300); %acceptability of the member wrt slenderness constraint
        while (~acceptable_slenderness) & (rnk_A<rnk_maxA) % go to the next section until the slenderness constraint is satisfied or the largest section is reached
            rnk_A=rnk_A+1;
            A=avail_sec(rnk_A,1);
            Ri=avail_sec(rnk_A,2);
            acceptable_slenderness=((L/Ri)<=300);
        end %either this member satisfies all member-based constraint or the cross section area cannot be increased anymore
    else %the axial force is compressive 
        tryNo=0;
        rnk_A=rnk_minA-1; %start with the smallest section (we  add one to this value in the next two lines)
        while rnk_A<rnk_maxA %search for the smallest cross section within the range [minA, maxA] that satisfies stress constraint
            rnk_A=rnk_A+1; %try the next cross section
            A=avail_sec(rnk_A,1);
            Ri=avail_sec(rnk_A,2);
            lambdam=L/Ri;
            cc=sqrt(2*pi^2*M_elasticity/Fy); %value of C_c according to AISC-ASD specifications 
            acceptable_stress=0;
            acceptable_slenderness=0;
            if lambdam<=200 % slenderness is satisfied
                acceptable_slenderness=1;
                tryNo=tryNo+1;
                if lambdam<=cc %inelastic buckling happens
                    numerator=(1-lambdam.^2/2/cc^2)*Fy;
                    denominator=5/3+3/8/cc*lambdam-lambdam.^3/8/cc^3;
                    sig_all=numerator/denominator; %allowable stress
                else %elastic buckling happens
                    sig_all=12*pi^2*M_elasticity/23/lambdam^2; %allowable stress
                end
                force_all=sig_all*A;% magnitude of the allowable axial force 
                if force_all>=abs(f_int_ext) %stress constraint is also satisfied
                    acceptable_stress=1;
                    break; % exit the while loop because the proper section 
                end
            end 
        end %end while: the proper rank has been found
    end
else %first load is tensile, the second one is compressive (case III)
    A=min(maxA,max(f_int_ext(1)/(.6*Fy),minA));
    acceptable_stress(1)=((f_int_ext(1)/(.6*Fy))<=maxA);
    [~,~,rnk_A]=round_A_W(A,1,section_no,1); % find the location of the new area in the given set 
    rnk_A=rnk_A-1;
    while rnk_A<rnk_maxA %search for the smallest cross section within the range [minA, maxA] that satisfies stress constraint
        rnk_A=rnk_A+1; %try the next cross section
        A=avail_sec(rnk_A,1);
        Ri=avail_sec(rnk_A,2);
        lambdam=L/Ri;
        cc=sqrt(2*pi^2*M_elasticity/Fy); %value of C_c according to AISC-ASD specifications 
        acceptable_slenderness=0;
        acceptable_stress(2)=0;
        if lambdam<200 % slenderness is satisfied
            acceptable_slenderness=1;
            if lambdam<=cc %inelastic buckling happens
                numerator=(1-lambdam.^2/2/cc^2)*Fy;
                denominator=5/3+3/8/cc*lambdam-lambdam.^3/8/cc^3;
                sig_all=numerator/denominator; %allowable stress
            else %elastic buckling happens
                sig_all=12*pi^2*M_elasticity/23/lambdam^2; %allowable stress
            end
            force_all=sig_all*A;% magnitude of the allowable axial force 
            if force_all>=abs(f_int_ext(2)) %stress constraint is satisfied
                acceptable_stress(2)=1;
                break; % exit the while loop because the proper section 
            end
        end 
    end %end while: the proper rank has been found
end
acceptable=min([acceptable_stress acceptable_slenderness]);
%[minA maxA minA0 maxA0 rnk_minA rnk_maxA  rnk1_ini Max_Avar Max_Avar Max_Avar_dec];


    
    


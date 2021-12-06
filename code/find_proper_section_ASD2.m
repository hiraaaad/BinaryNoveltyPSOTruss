% This function resizes a member for the given axial forces (f_int_ext)
% regarding the stress and slenderness ratio constraints according to
% AISC-ASD specifications (9th-edition)
function [A,Ri,rnk_A,ratios]=find_proper_section_ASD2(section_no,A,L,f_int_ext,Fy,M_elasticity,Max_Avar,Max_Avar_dec)
% member-based resizing is performed depending on which one is the case:
%i) member undergoes only tension (both numbers in f_int_ext are positive)
%ii)member undergoes only compression (both numbers in f_int_ext are negative)
%iii) member undergoes tension and compression (the first number in f_int_ext is positive and the second one is negative)
avail_sec=sections(section_no); %list of available sections (areas and radii of gyration)
minA=(A/Max_Avar_dec); % the lower bound on section area for member-based resizing
maxA=min(A*Max_Avar,max(avail_sec(:,1))); % the upper bound on section area for member-based resizing
[minA,~,rnk_minA]=round_A_W(minA,1,section_no,1); %round_up and find the rank of the minA in the given set
[maxA,~,rnk_maxA]=round_A_W(maxA,1,section_no,1); %round_up and find the rank of the maxA in the given set
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
        [A,Ri,rnk_A]=round_A_W(A,1,section_no,1); % round up the calculated area, find the location of the new area in the given set and the corresponding radius of gyration
        stress_ratio=(f_int_ext/A)/(.6*Fy); %acceptability of the member wrt stress constraint
        slender_ratio=(L/Ri)/300; %acceptability of the member wrt slenderness constraint
        while (slender_ratio>1) & (rnk_A<rnk_maxA) % go to the next section until the slenderness constraint is satisfied or the largest section is reached
            rnk_A=rnk_A+1;
            A=avail_sec(rnk_A,1);
            Ri=avail_sec(rnk_A,2);
            slender_ratio=(L/Ri)/300;
            stress_ratio=(f_int_ext/A)/(.6*Fy); %acceptability of the member wrt stress constraint
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
            if lambdam<=200 | (rnk_A>=rnk_maxA)% slenderness is satisfied or the upper limit of A is reached
                tryNo=tryNo+1;
                if lambdam<=cc %inelastic buckling happens
                    numerator=(1-lambdam.^2/2/cc^2)*Fy;
                    denominator=5/3+3/8/cc*lambdam-lambdam.^3/8/cc^3;
                    sig_all=numerator/denominator; %allowable stress
                else %elastic buckling happens
                    sig_all=12*pi^2*M_elasticity/23/lambdam^2; %allowable stress
                end
                force_all=sig_all*A;% magnitude of the allowable axial force 
                stress_ratio=abs(f_int_ext)/force_all;
                slender_ratio=lambdam/200;
                if stress_ratio<=1 %stress constraint is also satisfied
                     break; % exit the while loop because the proper section was found
                end
            end 
        end %end while: the proper section has been found
    end
else %first load is tensile, the second one is compressive (case III)
    A=min(maxA,max(f_int_ext(1)/(.6*Fy),minA));
    [~,~,rnk_A]=round_A_W(A,1,section_no,1); % find the location of the new area in the given set 
    stress_ratio(1)=(f_int_ext(1)/A)/(.6*Fy);
    rnk_A=rnk_A-1;
    while rnk_A<rnk_maxA %search for the smallest cross section within the range [minA, maxA] that satisfies stress constraint
        rnk_A=rnk_A+1; %try the next cross section
        A=avail_sec(rnk_A,1);
        Ri=avail_sec(rnk_A,2);
        lambdam=L/Ri;
        cc=sqrt(2*pi^2*M_elasticity/Fy); %value of C_c according to AISC-ASD specifications 
        if (lambdam<200) | (rnk_A>=rnk_maxA)% slenderness is satisfied ot the upper limit of A is reached
            slender_ratio=lambdam/200;
            if lambdam<=cc %inelastic buckling happens
                numerator=(1-lambdam.^2/2/cc^2)*Fy;
                denominator=5/3+3/8/cc*lambdam-lambdam.^3/8/cc^3;
                sig_all=numerator/denominator; %allowable stress
            else %elastic buckling happens
                sig_all=12*pi^2*M_elasticity/23/lambdam^2; %allowable stress
            end
            force_all=sig_all*A;% magnitude of the allowable axial force 
            stress_ratio(2)=abs(f_int_ext(2))/force_all;
            if stress_ratio(2)<=1 %stress constraint is satisfied
                break; % exit the while loop because the proper section has been found
            end
        end 
    end %end while: the proper rank has been found
end
ratios=[max(stress_ratio) max(slender_ratio)];


    
    


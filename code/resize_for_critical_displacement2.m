function [Anew,changeAind]=resize_for_critical_displacement(AAP,A,maxA,SEcoeff,Ucr,targetUcr,L_link,section_no) 
%This function solves the problem of finding whoch memebrs and how much to increase to reduce the critical displacement (U) to targetdelta
%A: current section area
%maxA: The upper limit on A determined by move limits and the avaialble sections
%SEcoeff: coefficient in save effectiveness, (Ff/E)
%Ucr: critical displacememnt to be reduced 
%targetUcr: the target value of the critical displacement
%L_link: length of the members
%section_no: section number list
initialCEcoeff=SEcoeff;
initialA=A;
newAcon=A;
changeAind=[]; %the members which are resized
if section_no>0 % discrete sections
    avail_sec=sections(section_no)';
    availA=avail_sec(1,:); % The available section list
    setsize=numel(availA); % the number of members in the list
end
maxiter=100;% The maximum numbe of iterations in bisection search to find the target save effectriveness
targetCE_base=.95; %<1, multiplier  when searching for CEtarget
SE=sign(Ucr).*(SEcoeff./initialA.^2).* (initialA<maxA); % save effectivenss of all members
maxSE=max(SE);
max_power=100;
min_power=0;


if  max(SE)>0 %it is possible to reduce the critical dispalcment 
    iteration=0;
    estimated_deltaU=0;
    while  iteration<maxiter %Now gradually decrease CEtarget unless the displacment constriants aren t satisfied  
        newAcon=initialA;
        target_power=(max_power+min_power)/2;
        SEtarget=maxSE*targetCE_base^target_power; % The upper limit for searching for CE_target
      %  SEtarget=targetCE_mult*SEtarget;
        iteration=iteration+1;
        changeAind=find(SE>SEtarget) ; %These members have smaller greated CS than CEtarget, therefore, their section should increase
        newAcon(changeAind)=sqrt(sign(Ucr)*SEcoeff(changeAind)./SEtarget); %new continuous values for the members with CE>CE_target
        newAcon(changeAind)=min([newAcon(changeAind);maxA(changeAind)]); %Apply the move limit to the increased section  
        estimated_deltaU=sum(SEcoeff(changeAind).* (1./newAcon(changeAind)-1./A(changeAind)).*L_link(changeAind)); %The estimated change in the current dispalcmement 
        if (max_power-min_power)<1
              break;
         end
         if (sign(Ucr)*(Ucr+estimated_deltaU-targetUcr)>0) % increase exponent to increase target SE
                    min_power=target_power;
        else
                    max_power=target_power;
        end
     %   [min_power target_power max_power];
        

    end
    %now SEtarget is good enough to reduced Ucr to targetUcr
    newAdis=newAcon;% find the proper discrete sections, since newAcon consists of continuous values
    rnk_newAdis= newAdis*0+1; %initial values (unimporatnt) of the location of the section values in the available section list
    [newAdis(changeAind),~,rnk_newAdis(changeAind)]=round_A_W(newAcon(changeAind),AAP(changeAind),section_no,zeros(1,numel((changeAind))));% round down the sections to the discrete values, nothing happens if sections can be continuous;
    % rnk_newAdis is the rank of member section in the given set
    % newAdis is the rounded section to the lower value
    if norm(newAdis(changeAind)-initialA(changeAind))>0 %if some change has occured, calculate deltaU caused by this change
        pridicted_deltaU=sum(SEcoeff(changeAind).* (1./newAdis(changeAind)-1./initialA(changeAind)).*L_link(changeAind));
    else
        pridicted_deltaU=0;
    end
    if section_no>0 & (sign(Ucr)*(pridicted_deltaU+Ucr-targetUcr))>0   % the sections are discrete and the target Ucr has noot been reached yet they are discrete, some of them should be rounded to th upper value (select the next larger section)
        newAdis0=newAdis; %some sections should chnage to the next available section in the list 
        next_rnk=min(setsize,rnk_newAdis+1); % next rank in the available section list
        nextAdis=availA(next_rnk); % the next section after newAdis
      %  priority_changeAind=(newAcon(changeAind)-newAdis(changeAind))./(nextAdis(changeAind)-newAdis(changeAind)).*(newAdis(changeAind)< maxA(changeAind));    % priority to rounding to the upper value
    
        priority_changeAind=-sign(Ucr)*SEcoeff(changeAind)./(nextAdis(changeAind)-newAdis(changeAind)+1e-100).*(1./nextAdis(changeAind)-1./newAdis(changeAind))     .*  (newAdis(changeAind)< maxA(changeAind))   ;
       
 [~,ind5]=sort(priority_changeAind,'descend');
        ind6=changeAind(ind5); % priority of members for increasing to the next avaialoble value
        clear ind5
        no_of_candid_max=numel(priority_changeAind); % maximum number of candiate sections to increase 

        for  candid_no=1:no_of_candid_max
        
            pridicted_deltaU=pridicted_deltaU+(SEcoeff(ind6(candid_no))   .* (1./nextAdis(ind6(candid_no))-1./newAdis0(ind6(candid_no)))  .*   L_link(ind6(candid_no))); %apply the effecr of increase of section to the criticval dispalcememtn conmstriant
         %  innew= [candid_no  pridicted_deltaU pridicted_deltaU+Ucr targetUcr ]
           newAdis(ind6(candid_no))=nextAdis(ind6(candid_no)); %change no_of_candid sections with highest priority to the next available section 
       %     pridicted_deltaU4=sum(SEcoeff   .* (1./newAdis) .*   L_link)-Ucr;

        %    if abs((pridicted_deltaU4-pridicted_deltaU4)/Ucr)>1e-6  % This is a checkpoint only: error in calculation of predicted_deltaU
         %       'eshkal'
         %       ind6
        %        [pridicted_deltaU0+Ucr sum(SEcoeff   .* (1./newAdis0) .*   L_link) pridicted_deltaU4 pridicted_deltaU Ucr  sum(SEcoeff   .* (1./initialA) .*   L_link)]
        %        format longe
        %        [pridicted_deltaU4 pridicted_deltaU Ucr]
        %        'warning in resizse subproblem'
      %          erorrrrrr
     %       end
            if (sign(Ucr)*(pridicted_deltaU+Ucr-targetUcr))<=0 % enough
                   break;
             end

        end % The required number of sections were rounded to the upper value, based on their priority
    end %sections are discrete now
end % Ucr should be targetUcr unless the maximum move limit is reached or targetSE limit is reached 
Anew=initialA;
if numel(changeAind)>0 % some member sections have changed so far, update them
    Anew(changeAind)=newAdis(changeAind);
end
if 0 % these are only checkpoints for the algorithm
    maxpossibleDdelta4= sum(-sign(Ucr)*((SEcoeff./maxA) - (SEcoeff./A)) .*  (Anew<maxA)); %if zero, def cannot get better
    computed_Def=sum(SEcoeff./(Anew).*L_link);
    if abs(computed_Def-(pridicted_deltaU+Ucr))/Ucr>1e-6 % if condition, then an error has occured in computation of dispalcmenterror 
        'tafavot dar subproblem15'
        format longe
        [computed_Def pridicted_deltaU  Ucr  pridicted_deltaU+Ucr pridicted_deltaU0]
        N_change
        format shorte
        initialCEcoeff-SEcoeff
        'warning, these walues should  be equal'
        %erororrrrr
    end
    SE=(  sign(Ucr).*  (SEcoeff./Anew.^2)   .* (Anew<maxA));
    if maxpossibleDdelta4>0  &  (sign(Ucr)*(pridicted_deltaU+Ucr-targetUcr))>0 & max(SE>SEtarget) % checkpoint: must not happen
        'chera nashod bishatresh konim?'
        maxpossibleDdelta4
        [pridicted_deltaU Ucr  targetUcr]
        SE=(  sign(Ucr).*  (SEcoeff./Anew.^2)   .* (Anew<maxA));
        sort(SE)
        'eshkal dar resize subproblem'
        [SEtarget max(SE)]   
        [~,ind22]=max(SE)
        [ind22 A(ind22) AAP(ind22) initialA(ind22) maxA(ind22) targetUcr pridicted_deltaU+Ucr]
        ind6(1:no_of_candid)
        errorrrrr  
    end
end
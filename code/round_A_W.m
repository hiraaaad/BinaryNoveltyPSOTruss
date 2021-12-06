function [A_round,Ri_round,keeprnk]=round_A_W(A,AAP,list_no,W)
%w is a vector f size A= 0:down, 1:up, .5: uniform stochastic
W = double(W); % Hirad addition

A-AAP-W; %ignore this line
if max(abs(W-.5))>.5 %checkpoint
    'error in round_A_W'
    'rounding option is not supported'
    'Aborting the process ...'
    error_occured
end

A_round=A;
Ri_round=A;
keeprnk=A;
if list_no>0
    
    list=sections(list_no);
    availA=list(:,1)';
    Asize=numel(A);
    listsize=numel(availA);
    AAP_active_ind=find(AAP);
    for k=AAP_active_ind
        if A(k)<=availA(1)
            rnk=1;
        elseif A(k)>=availA(end)
            rnk=listsize;
        else
            tavan=abs(log(.5)./log(W(k)));
            rnk=sum(availA<=A(k));
            p=(A(k)-availA(rnk))./(availA(rnk+1)-availA(rnk));
            rnk=rnk+ (p>(rand^tavan));
        end
        A_round(k)=availA(rnk);
        if list_no>100
            Ri_round(k)=list(rnk,2);
        end
        keeprnk(k)=rnk;
    end    
end


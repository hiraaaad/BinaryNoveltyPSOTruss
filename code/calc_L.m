function L_link=calc_L(AAP,nodes,TA,D)

X=nodes(1:D:end);
Y=nodes(2:D:end);
if D==3
Z=nodes(3:D:end);
end
L_link=AAP*0;
Kg=zeros(numel(X)*D,numel(X)*D); % initial values of global stiffness matrix
if D==3
    for r=1:numel(AAP)
        if AAP(r)>.99 % if exists
            inds=TA(r,:);
            m=inds(1);
       	    n=inds(2);
            L_link(r)=norm([X(n)-X(m) Y(n)-Y(m) Z(n)-Z(m)]);
        end
    end
elseif D==2
    for r=1:numel(AAP)
        if AAP(r)>.99
            inds=TA(r,:);
            m=inds(1);
            n=inds(2);
            L_link(r)=norm([X(m)-X(n) Y(m)-Y(n)]);
          
         end
    end
end
 
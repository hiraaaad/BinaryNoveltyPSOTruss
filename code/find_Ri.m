function Ri=find_Ri(A,section_no)
sec=sections(section_no);
for k=1:numel(A)
   [dis,ind]=min(abs(A(k)-sec(:,1)));
   Ri(k)=sec(ind,2);
end


function [NAP,NAPconnector]=find_NAP(TA,AAP)

NAP=(1:max(max(TA)))*0;
NAPconnector=NAP*0;
  for k=1:numel(AAP)
       if AAP(k)==1
            NAP(TA(k,:))=1;
            NAPconnector(TA(k,:))=NAPconnector(TA(k,:))+1;


      end
   end

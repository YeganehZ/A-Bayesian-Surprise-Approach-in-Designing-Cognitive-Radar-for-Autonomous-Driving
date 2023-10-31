function [R_Localized, R_Localized_indx]=nearest_neighbour(landa,beta,NRR,R_Library,R_Loc)

NRR_b=ceil(NRR/length(landa));
NRR_l=NRR-(length(landa)*(NRR_b-1));

X=[NRR_b ,NRR_l];
for i=1:length(beta)
   for j=1:length(landa)
     Y=[i j];  
     NY(j,:)=Y;
   end
   NYY(:,:,i)=NY;
end

cc=1;
for k=1:length(beta)
NYYY(1+(length(landa)*(cc-1)):(length(landa)*cc),:)=NYY(:,:,k);
cc=cc+1;
end

 indices=knnsearch(NYYY,X,'K',R_Loc);
 %indices(:,1)=[];
 nearestpoints=NYYY(indices,:);
 
R_Localized_indx=indices;
R_Localized=R_Library(:,:,indices);
%Plot_localized_actions(landa,beta,nearestpoints,X)




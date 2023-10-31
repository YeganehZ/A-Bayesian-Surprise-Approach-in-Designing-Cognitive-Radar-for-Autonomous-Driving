function ExpSurprise=Planning(m,cycle,Pest,H,Q,F,R_Localized)

R_Loc=size(R_Localized,3);
Spred_new=zeros(m,m,R_Loc); 
ExpSurprise=zeros(1,R_Loc);

Ppred = F*Pest*F' + Q;

 for i=1:R_Loc
  Spred_new(:,:,i)=R_Localized(:,:,i) + H*Ppred*H';
  ExpSurprise(i)= trace(pinv(R_Localized(:,:,i)*pinv(Spred_new(:,:,i))))+ 0.5*log(det(R_Localized(:,:,i)*pinv(Spred_new(:,:,i))))-m;
 end


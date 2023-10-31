function [NORM1est_ind,NORM1_ind,NORM2est_ind,NORM2_ind] = Performance_Test(K,n,x, xest)

NORM1est_ind=zeros(n,K);
NORM1_ind=zeros(n,K);

    
for j=1:K
  for i=1:n    
    NORM1est_ind(i,j)=sqrt((x(i,j) - xest(i,j))'*(x(i,j) - xest(i,j)));
    NORM1_ind(i,j)=sqrt(x(i,j)'*x(i,j));
  end
end
    
NORM2est_ind=zeros(n,K);
NORM2_ind=zeros(n,K);

for j=1:K
   for i=1:n 
    NORM2est_ind(i,j)=(x(i,j) - xest(i,j))'*(x(i,j) - xest(i,j));
    NORM2_ind(i,j)=x(i,j)'*x(i,j);
   end
end

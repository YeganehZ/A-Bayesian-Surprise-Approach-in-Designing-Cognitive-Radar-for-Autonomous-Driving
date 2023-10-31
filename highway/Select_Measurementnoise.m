function R_library=Select_Measurementnoise(landa,beta,dx,B)

%landa=0.4*10^(-5);
%beta=10^(10);

c=3*10^8;
fc=77*10^9;
R_library=zeros(2,2,length(landa)*length(beta));
d0=2000;
dy=0; d=sqrt(dy^2+dx^2);
SNR=(d0/d)^4;

R= [(c)^2/(2*SNR),  -(c^2 * B)/(2*pi*fc*SNR);...
    -(c^2* B)/(2*pi*fc*SNR), (c^2/((2*pi*fc)^2*SNR))]; 

for i=1:length(landa)
  for j=1:length(beta)
    R_library(:,:,length(beta)*(i-1)+j) = [R(1,1)*(landa(i))^2,  -(beta(j)*(c*landa(i))^2)/(2*pi*fc*SNR);...
                         -(beta(j)*(c*landa(i))^2)/(2*pi*fc*SNR),   R(2,2)*((1/(2*landa(i)^2))+(2*beta(j)^2*landa(i)^2) ) ]; 
  end
end





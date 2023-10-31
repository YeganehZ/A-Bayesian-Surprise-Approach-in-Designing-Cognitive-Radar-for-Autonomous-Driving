clc
clear all

% This program implements the system performance of our proposed cognitive 
% radar algorithm based on the Expectation of Bayesian Surprise for the highway 
% driving scenario with constant acceleration 


K=100;  % time index 
MC=100;  % number of Monte Carlo runs

delta=0.1; % sample time
q=0.01;    % state noise variance

% Transition Matrix
F = [1 delta 0 0 0; ...
    0 1 0 0 0;...
    -delta -delta^2/2 1 delta delta^2/2;...
    0 0 0 1 delta;...
    0 0 0 0 1]; 

% Measurement Matrix
H=[0 0 1 0 0; 0 0 0 1 0]; 

% State Noise Covariance
Q = [(delta^4)/4    (delta^3)/2   -(delta^5)/12  0   0; ...
     (delta^3)/2     delta^2     -(delta^4)/6   0   0; ...
    -(delta^5)/12  -(delta^4)/6   (delta^6)/18  (delta^5)/12   (delta^4)/6;...
        0            0              (delta^5)/12   (delta^4)/4   (delta^3)/2;...
        0            0              (delta^4)/6   (delta^3)/2     delta^2]*q; 

    
n=size(H,2);   % Dimension of state vector, state estimate, .... 
m=size(H,1);   % Dimension of the measurement vector, innovation sequence

% Pre-setting the vectors and matrices:

x=zeros(n,K);     % state vector
z=zeros(m,K);     % measurement vector

xpred=zeros(n,K); % predicted state vector k|k-1
xest=zeros(n,K);  % estimate state vector k|k


% Error Covariance
Ppred=zeros(n,n,K);   % predicted error covariance P(k|k-1)
Pest=zeros(n,n,K);    % estimate error covariance at P(k|k) 
Spred=zeros(m,m,K);   % predicted innovation covariance matrix 
Sest=zeros(m,m,K);    % estimated innovation covariance matrix
nu=zeros(m,K);        % innovation vector
KG=zeros(n,m,K);      % Kalman Gain matrix


NORM1est_ind_MC=zeros(n,K,MC);  % defines the norm 2 vector between xest and x for all monte-carlo simulation 
NORM1_ind_MC=zeros(n,K,MC);     % defines the norm 2 vector of x for all monte-carlo simulation 
NORM2est_ind_MC=zeros(n,K,MC);  % defines the norm 2 vector between xest and x for all monte-carlo simulation 
NORM2_ind_MC=zeros(n,K,MC);     % defines the norm 2 vector of x for all monte-carlo simulation 
RMSRE_ind=zeros(n,K);           % Root Mean Square Relative Error (RMSRE) (without considering consistancy test)
AEE_ind=zeros(n,K);             % Average Euclidean relative error (ARE)
GAE_ind=zeros(n,K);             % Harmonic average relative error (HRE)
HAE_ind=zeros(n,K);             % Geometric average relative error (GRE)

x0 =[25; 3; 100; 23; 2];      % true initial state vector
xpred0=[24; 3; 80; 23; 2];    % initial state estimation mean vector 
Ppred0=diag([100,1,100, 100, 1]);    % initial state estimation covariance matrix
distance=100;    % initial distance between two vehicles
B=10^(8);        % Bandwidth 

landa=[10^(-6):0.1*10^(-6):10^(-5), 1.1*10^(-5):0.1*10^(-5):10^(-4)] ; % pulse duration values
beta= [-10^(12) :0.2*10^(12): -0.2*10^(12),   0.2*10^(12):0.2*10^(12):10^(12)]; % chirp rate values
R_library=Select_Measurementnoise(landa,beta,distance,B);       % This function builds the whole action library (measurement 
                                                                % noise covariance library)
NR = 1: (length(landa)*length(beta));      % counts the action library 


SEED=200;   % setting the seed of the noise 
rng(SEED);

for mc=1:MC    

% Select initial action from library
NRR = NR(randperm(length(NR),1));  % Select the initial action index randomly based on uniform distirbution 
R=R_library(:,:,NRR);              % set the action

% Set initial conditions of the filter
x(:,1)= x0;
xpred(:,1)= xpred0;
xest(:,1)= xpred0;
Ppred(:,:,1) =Ppred0;
Pest(:,:,1) =Ppred0;
Spred(:,:,1)=R + H*Pest(:,:,1)*H';
KG(:,:,1)=Pest(:,:,1)*H'*pinv(R);
z(:,1)=H*x(:,1);
nu(:,1)=z(:,1)-H*xpred(:,1);

Rindx=zeros(1,K);  % record the index of the selected actions at each time instant 
Rindx(1,1)=NRR;    % the index of the intial action 
R_Loc=25; % Size of the localized action library (number of actions in the localized action library)


% build the state vector
  for j=1:K-1
    x(:,j+1)=F*x(:,j)+ mvnrnd(zeros(1,n),Q)';  
  end
  
 for j=2:K

  z(:,j)=H*x(:,j)+ mvnrnd(zeros(1,m),R)'; % Measurement noise vector 


  %% Kalman Filter Estimation and Prediction

  [xpred(:,j),Ppred(:,:,j),Spred(:,:,j)]=timeupdate(xest(:,j-1),Pest(:,:,j-1),R,H,Q,F); 
  [xest(:,j), Pest(:,:,j), KG(:,:,j),nu(:,j)] = measurementupdate(xpred(:,j),Ppred(:,:,j),Spred(:,:,j),z(:,j),H);  

  % creates a localized action library using k-nearest neighbor algorithm
  [R_Localized,  R_Localized_indx]=nearest_neighbour(landa,beta, Rindx(1,j-1),R_library,R_Loc); % localize action

   %% Planning 
   ExpSurprise=Planning(m,j,Pest(:,:,j),H,Q,F,R_Localized);
 
   %% Action Selection
  
    Rnewindx=find(ExpSurprise == max(ExpSurprise(:)));
    R=R_library(:,:,Rnewindx);
    Rindx(1,j)= R_Localized_indx(1,Rnewindx);

  end
 
 %% Computing the second and first norm to calculate RMSRE, HRE, ARE, GRA
  
[NORM1est_ind,NORM1_ind,NORM2est_ind,NORM2_ind] = Performance_Test(K,n,x, xest); 
NORM1est_ind_MC(:,:,mc)=NORM1est_ind;
NORM1_ind_MC(:,:,mc)=NORM1_ind; 
NORM2est_ind_MC(:,:,mc)=NORM2est_ind;
NORM2_ind_MC(:,:,mc)=NORM2_ind;  

mc  
end

%% Compute the RMSRE and Absolute RMSRE
NOM1=sum(NORM2est_ind_MC,3); 
DIM1=sum(NORM2_ind_MC,3);
for j=1:K
    for i=1:n
    RMSRE_ind(i,j)=sqrt(NOM1(i,j))/sqrt(DIM1(i,j));
    end
end

ABS_RMSRE_ind=(1/(K-1))*sum(RMSRE_ind(:,2:K),2);


%% Compute the ARE and Absolute ARE
NOM2=sum(NORM1est_ind_MC,3); 
DIM2=sum(NORM1_ind_MC,3);
for j=1:K
    for i=1:n
    AEE_ind(i,j)=NOM2(i,j)/DIM2(i,j);
    end
end

ABS_AEE_ind=(1/(K-1))*sum(AEE_ind,2);


%% Compute the GRA and Absolute GRA

NOM3=sum(log(NORM1est_ind_MC),3); 
DIM3=sum(log(NORM1_ind_MC),3); 

for j=1:K
    for i=1:n
    GAE_ind(i,j)=(1/MC).*(NOM3(i,j)-DIM3(i,j));
    end
end

ABS_GAE_ind=(1/(K-1))*sum(GAE_ind(:,2:K),2);

%% Compute the HRA and Absolute HRA

NOM4=sum((NORM1est_ind_MC).^(-1),3); 
DIM4=sum((NORM1_ind_MC).^(-1),3); 

for j=1:K
    for i=1:n
      HAE_ind(i,j)=(NOM4(i,j))^(-1)/(DIM4(i,j))^(-1);
    end
end

ABS_HAE_ind=(1/(K-1))*sum(HAE_ind(:,2:K),2);


f=2:K;  
figure, 
p= plot(f, (RMSRE_ind(1,f)), 'b -', f, AEE_ind(1,f),'b --',f, HAE_ind(1,f),'b :');
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
leg2 = legend('\textbf{RMSRE}', '\textbf{ARE}', '\textbf{HRE}');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',22);
xlabel({'\textbf{Time Samples($k$)}'},'Interpreter','latex','FontSize',22,'Color','k'); 
ylabel('$\mathbf{E}[\mathcal{S}^{B}_{k}] (v^0_x)$','Interpreter','latex','FontSize',22,'Color','k');
set(gca,'FontSize',22)
set(gca,'FontName','Times New Roman')
set(gcf,'color','white')
hold on      

figure, 
p= plot(f, (RMSRE_ind(3,f)), 'b -', f, AEE_ind(3,f),'b --',f, HAE_ind(3,f),'b :');
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
leg2 = legend('\textbf{RMSRE}', '\textbf{ARE}', '\textbf{HRE}');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',22);
xlabel({'\textbf{Time Samples($k$)}'},'Interpreter','latex','FontSize',22,'Color','k'); 
ylabel('$\mathbf{E}[\mathcal{S}^{B}_{k}] (d_x)$','Interpreter','latex','FontSize',22,'Color','k');
set(gca,'FontSize',22)
set(gca,'FontName','Times New Roman')
set(gcf,'color','white')
hold on  

figure, 
p= plot(f, (RMSRE_ind(4,f)), 'b -', f, AEE_ind(4,f),'b --',f, HAE_ind(4,f),'b :');
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
leg2 = legend('\textbf{RMSRE}', '\textbf{ARE}', '\textbf{HRE}');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',22);
xlabel({'\textbf{Time Samples($k$)}'},'Interpreter','latex','FontSize',22,'Color','k'); 
ylabel('$\mathbf{E}[\mathcal{S}^{B}_{k}] (v^1_x)$','Interpreter','latex','FontSize',22,'Color','k');
set(gca,'FontSize',22)
set(gca,'FontName','Times New Roman')
set(gcf,'color','white')
hold on      







function [Xest, Pest, KG,nu] = measurementupdate(Xpred,Ppred,Spred,z,H)
KG = Ppred*H'*pinv(Spred); %% Kalman gain 
nu= z - H * Xpred;
Xest= Xpred + KG * nu;
Pest= Ppred - KG * Spred * KG';

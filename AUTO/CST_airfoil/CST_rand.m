function [wlr,wur,coord] = CST_rand(wl,wu,Sigma,N)
% CST Airfoil Parameterization Randomization
% Inputs: 
%        N: Number of randomized correlated coefficients
%       wl: Mean lower surface CST coefficients
%       wu: Mean upper surface CST coefficients
%    sigma: Covariance matrix
%
% Outputs:
%      wlr: iid Lower surface CST coefficients
%      wur: iid Upper surface CST coefficients
%    coord: Random coordinate collection from CST parameterization

mu = [wl,wu];
R = chol(Sigma);
z = repmat(mu,N,1) + randn(N,length(mu))*R;

wlr = z(:,1:length(wl));
wur = z(:,length(wl)+1:end);

for i = 1:N 
    coord(:,:,i) = CST_airfoil(wlr(i,:),wur(i,:),0,200);
end
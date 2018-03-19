function [betaK] = saen(y,X,K,AL)
%-< Sequential adaptive elastic net (SAEN) approach >---- 
% Finds the SAEN estimate using c-PW-WEN algorithm.
%-< Outputs >-----
% betaK = Final K-sparse SAEN estimate
% Code by Muhammad Naveed Tabassum and Esa Ollila, Aalto University. <2018>

if nargin < 4, AL = 1; end % Default case is Lasso (\alpha=1)
[~,p] = size(X);
betaK = zeros(p,1);
w = ones(p,1);      % initial weights
supp = 1:p;         % initial support set
for itr=3:-1:1
    [beta0,~,A] = cpwwen(y,X(:,supp), itr*K, AL,w, itr==1);
    nzi = A(1:end-1);
    w = 1 ./ abs(beta0(nzi));	% current weights
    supp = supp(nzi);           % current support set
end
betaK(supp) = beta0(nzi);    % Final K-sparse solution

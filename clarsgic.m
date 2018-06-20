function [estK, A,beta] = clarsgic(y,X, A,B,gama,maxp)
%< Complex-valued LARS with generalized information criterion (GIC) for
%  finding the sparsity level and respective solution.
% Code by Muhammad Naveed Tabassum and Esa Ollila, Aalto University. <2018>

[n,p] = size(X);
if nargin < 6, maxp = min(n,p); end
if nargin < 5, gama = 2; end

scale_est = zeros(maxp,1);    % Unbiased scale estimate
scale_est(1) = (y'*y)/n;      % For null, i.e., all zeros case
nzB = B(:,2:end)~=0;
beta = zeros(p,1);

for k = 1:maxp-1
    Xa = X(:, nzB(:,k));    % active-set of predictors
    b_ls = B(nzB(:,k),k+1);	% LSE for current active-set | Xa \ y
    r = y - Xa*b_ls;        % Residual for current active-set
    RSS = r'*r;
    scale_est(k+1) = RSS / (n-k);
end
switch gama
    case 0
        c = log(n);                 % c-LARS-GIC0
    case 1
        c = log(n) * log(log(p));	% c-LARS-GIC1
    otherwise
        c = log(p) * log(log(n));	% c-LARS-GIC2
end
GIC = n*log(scale_est) + (0:maxp-1)'*c;	
indx = find(GIC<=min(GIC), 1);
estK = indx - 1;    % since we start from 0
A = A(1:estK);
beta(A) = X(:, A) \ y;
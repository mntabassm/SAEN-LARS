function [betaK,Lam,A] = cpwwen(y,X,K,AL,w,debias)
%-< Complex-valued Path Wise Weighted Elastic Net (c-PW-WEN) >----
% Finds the knots (and respective solution) for WEN using the c-LARS-WLasso.
%-< Outputs >-----
% betaK = Final K-sparse WEN estimate
% Lam = {lam0,...,lamK}, i.e., set of knot-values from lam0 to lamK
% A = support set in the predictor-entering order
% Code by Muhammad Naveed Tabassum and Esa Ollila, Aalto University. <2018>

[n, p] = size(X);
if nargin < 6, debias = 'true'; end	% Default case is to debias the solution
if nargin < 5, w = ones(p,1); end	% Default case is uniformly weighted
if nargin < 4, AL = 1; end          % Default case is Lasso (\alpha=1)

maxp = min(K+1, n);
m = length(AL);
RSS = zeros(1, m);
Lam = zeros(maxp, m);    % EN knots for m-alphas
A = zeros(maxp, m);      % Predictor-entering orders for m-alphas
B = zeros(p, m);

[lam,A(:,1),B1] = clarswlasso(y,X,0,0,maxp,w);
Lam(:,1) = lam;     B(:,1) = B1(:,maxp);
lamen0 = lam;       A(1,:) = A(1,1); % 1st entering var is same for all-alphas
Xa = X(:, abs(B(:,1))~=0);	% active-set of predictors
b_ls = Xa \ y;              % LSE for current active-set
r = y - Xa*b_ls;            % Residual for current active-set
RSS(1) = r'*r;
ya = [y; zeros(p,1)]; % augmented form of y

for ii= 2:m
    al = AL(ii);
    Aen = zeros(maxp-1,1);	Ben = zeros(size(B1));
    lamen = zeros(maxp,1);	lamen(1) = lam(1)/al; % lam(alpha=1)/alpha_i
    for jj=2:maxp
        Xa = [X; sqrt(lamen0(jj)*(1-al))*eye(p)]; % augmented form of X
        [tmpL,tmpA,tmpB] = clarswlasso(ya,Xa,0,0,jj,w);
        lamen(jj) = tmpL(jj)/al;        Aen(jj-1) = tmpA(jj);
        eta_k = lamen(jj)*(1-al);
        Ben(:,jj) = (1+eta_k)*tmpB(:,jj); % correction for double shrinkage
    end
    A(2:end,ii) = Aen;      B(:,ii) = Ben(:,maxp);
    Lam(:,ii) = lamen;      lamen0 = lamen;
    
    Xa = X(:, abs(B(:,ii))~=0);	% active-set of predictors
    b_ls = Xa \ y;              % LSE for current active-set
    r = y - Xa*b_ls;            % Residual for current active-set
    RSS(ii) = r'*r;
end
ali = find(RSS<=min(RSS),1); % best \alpha's index
betaK = B(:,ali);	% Final K-sparse solution
Lam = Lam(:,ali);	% Knots (lambda values)
A = A(:,ali);       % Order of predictors in which they entered in solution

if debias % debiased solution using LSE of active predictors
    betaK(A(1:end-1)) = X(:,A(1:end-1)) \ y;
end

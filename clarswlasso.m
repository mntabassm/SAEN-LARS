function [lam,A,B1] = clarswlasso(y,X, intcpt,scaleX,maxp,w)
%-< Complex-valued LARS for Weighted Lasso (c-LARS-WLasso) >----
% Finds the knots (and respective solutions) for WLasso
% using the complex-valued extention of LARS method.
%-< Outputs >-----
% lam = {lam0,...,lamK}, i.e., set of knot-values from lam0 to lamK
% A = support set in the predictor-entering order
% B1 = set of solutions using {lam0,...,lamK} for Lasso(alpha=1 in EN)
% Code by Muhammad Naveed Tabassum and Esa Ollila, Aalto University. <2018>

[n, p] = size(X);

if nargin < 6, w = ones(p,1); end
if nargin < 5, maxp = min(n,p); end
if nargin < 4, scaleX=1; end
if nargin < 3, intcpt=1; end

if intcpt % if intercept is in the model, center the data
    meanX= mean(X); meany = mean(y);
    X = bsxfun(@minus, X, meanX);
    y = y-meany;
end
if scaleX
    % Default: scale X to have column norms eq to sqrt(n) (n = nr of rows)
    sdX = sqrt(sum(X.*conj(X)));
    X   = bsxfun(@rdivide, X, sdX);
end

X = X./repmat(w',n,1);  % for weighted Lasso from Adaptive Lasso paper (section 3.5)
lam = zeros(1,maxp);
B1 = zeros(p,maxp);
beta = zeros(p,1);
r  = y - X*beta; % residual vector
c = X'*r;        % correlation of predictors with residual vector
[lam0, A] = max(abs(c));  % compute lambda0
lam(1)=lam0;
Delta = zeros(p,1);
for k=2:maxp
    delta = (1/lam(k-1))*(X(:,A) \ r) ;
    Delta(A) = delta;
    notA = setdiff(1:p,A);
    al = zeros(1,p);
    for ell=notA
        cl = X(:,ell)'*r;
        bl = X(:,ell)'*(X(:,A)*delta);
        ax2 = abs(bl)^2 - 1;
        bx1 = 2*(lam(k-1)-real(cl*conj(bl)));
        cx0 = abs(cl)^2 - lam(k-1)^2;
        qp = [ax2, bx1, cx0];
        ALz = roots(qp);
        
        if all(ALz>0)
            al(ell)  = min(ALz);
        else
            al(ell)  = subplus(max(ALz));
        end
        
        if ~isreal(ALz)
            fprintf('complex root\n');
            al(ell) = 0;
        end
        c(ell) = cl - al(ell)*bl;
    end
    [almin, tmp] = min(al(notA));
    j = notA(tmp);
    lam(k) = lam(k-1) - almin;
    A = [A, j];
    beta = beta + almin*Delta;
    B1(:,k) = beta;
    r = y - X*beta;
end
B1 = B1./repmat(w,1,maxp); % from Adaptive Lasso paper (section 3.5)
if scaleX, B1 = bsxfun(@rdivide,B1, sdX'); end
if intcpt, B1 = [meany - meanX*B1;  B1]; end

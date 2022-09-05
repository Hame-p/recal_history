function [b,Stats] = ck_stat_regress(y,X)


% function  [beta,Stats] = ck_stat_regress(Y,X)
%
% Stats.model = [r2 F prob s2 SSE RSS TSS LL];
% Stats.B_T t-values for each beta
% Stats.B_P p-values 

[n,ncolX] = size(X);
if ~isvector(y) || numel(y) ~= n
    error(message('stats:regress:InvalidData'));
end

% Use the rank-revealing QR to remove dependent columns of X.
[Q,R,perm] = qr(X,0);
p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
if p < ncolX
    warning(message('stats:regress:RankDefDesignMat'));
    R = R(1:p,1:p);
    Q = Q(:,1:p);
    perm = perm(1:p);
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(ncolX,1);
b(perm) = R \ (Q'*y);


nu = max(0,n-p);                % Residual degrees of freedom
yhat = X*b;                     % Predicted responses at each data point.
r = y-yhat;                     % Residuals.
normr = norm(r);
if nu ~= 0
  rmse = normr/sqrt(nu);      % Root mean square error.
else
  rmse = NaN;
end
s2 = rmse^2;                    % Estimator of error variance.

% There are several ways to compute R^2, all equivalent for a
% linear model where X includes a constant term, but not equivalent
% otherwise.  R^2 can be negative for models without an intercept.
% This indicates that the model is inappropriate.
SSE = normr.^2;              % Error sum of squares.
RSS = norm(yhat-mean(y))^2;  % Regression sum of squares.
TSS = norm(y-mean(y))^2;     % Total sum of squares.
r2 = 1 - SSE/TSS;            % R-square statistic.
if p > 1
  F = (RSS/(p-1))/s2;      % F statistic for regression
else
  F = NaN;
end
prob = fpval(F,p-1,nu); % Significance probability for regression

% log likelihood
sigma = std(r);
LL = -n*log(2*pi)/2-n*log(sigma^2)/2-(1/(2*sigma^2))*(yhat-y)'*(yhat-y);


Stats.model = [r2 F prob s2 SSE RSS TSS LL];

% compute SE of each parameter
RI = R\eye(p);
B_SE(perm,:) = rmse*sqrt(sum(abs(RI).^2,2));
B_T = b ./ B_SE;
B_P = 2 * normcdf(-abs(B_T));
Stats.B_T = B_T;
Stats.B_P = B_P;


return;
function p = fpval(x,df1,df2)
%FPVAL F distribution p-value function.
%   P = FPVAL(X,V1,V2) returns the upper tail of the F cumulative distribution
%   function with V1 and V2 degrees of freedom at the values in X.  If X is
%   the observed value of an F test statistic, then P is its p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FINV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 2010 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2010/11/08 02:37:46 $

if nargin < 3, 
    error(message('stats:fpval:TooFewInputs')); 
end

xunder = 1./max(0,x);
xunder(isnan(x)) = NaN;
p = fcdf(xunder,df2,df1);

    
function[Lo,Up,pval]=ck_stat_bootCI_percentile(x,statfun,alpha,B1,xnull)
%           
%      [Lo,Up]=confintp(x,statfun,alpha,B1,PAR1,...)
%
%      Confidence interval of the estimator of a parameter
%      based on the bootstrap percentile method  
%
%     Inputs:
%           x - input vector data 
%     statfun - the estimator of the parameter given as a Matlab function   
%      alpha  - level of significance (default alpha=0.05)  
%          B1 - number of bootstrap resamplings (default B1=199)   
%    PAR1,... - other parameters than x to be passed to statfun
%
%     Outputs:
%         Lo - The lower bound 
%         Up - The upper bound
%
%     Example:
%
%     [Lo,Up] = confintp(randn(100,1),'mean');


%  Created by A. M. Zoubir and  D. R. Iskander
%  May 1998
%
%  References:
% 
%  Efron, B.and Tibshirani, R.  An Introduction to the Bootstrap.
%               Chapman and Hall, 1993.
%
%  Hall, P. Theoretical Comparison of Bootstrap Confidence
%               Intervals. The Annals of Statistics, Vol  16, 
%               No. 3, pp. 927-953, 1988.
%
%  Zoubir, A.M. Bootstrap: Theory and Applications. Proceedings 
%               of the SPIE 1993 Conference on Advanced  Signal 
%               Processing Algorithms, Architectures and Imple-
%               mentations. pp. 216-235, San Diego, July  1993.
%
%  Zoubir, A.M. and Boashash, B. The Bootstrap and Its Application
%               in Signal Processing. IEEE Signal Processing Magazine, 
%               Vol. 15, No. 1, pp. 55-76, 1998.


if (exist('B1')~=1), B1=199; end;
if (exist('alpha')~=1), alpha=0.05; end;

vhat=feval(statfun,x);
vhatstar=bootstrp(B1,statfun,x);

q1=floor(B1*alpha*0.5);
q2=B1-q1+1;
st=sort(vhatstar);
Lo=st(q1);
Up=st(q2);


pval = min(sum(vhatstar<xnull),sum(vhatstar>xnull))/B1;
pval = 2*pval;
function [B,stats,stats_each,contrast] = ck_stat_glm(X,Y,C)
% Glm Interface with contrasts
%
% function [B,stats,stats_each,contrast] = ck_stat_glm(X,Y,C)
%
% compute GLM outputs for data Y and predictors X
% Y = [nsamples]
% X = [nsamples npredictors] with 1st predictor = constant
%
% C Contrast vector between variables or cell array of contrasts
%
% Y = X*B
% uses pinv to solve problem
%
% stats_each contains the partial and semi-partial correlations and
% asociated stats.
% Partial correlation:  relative contribution of each regressor to the data 
% (how much each explains of y when we control for the influence of all other xs on both y and xi)
%
%
% Semipartial correlation: relative contribution of each regressor to the total variance 
% (how much x1 explains of y when we control for the influence of all other xs on x1 only) 
% -> measure of effect size, or unique variance explained
%  stats_each.semiP(var,:) = [semiR2,F,p];


s = size(X);
if rank(X)~=s(2)
%  fprintf('Model is rank deficient\n');
end

B = pinv(X)*Y; % inv(X'*X)*X'*Y


Yhat = X*B;
Res = Y - Yhat;
stats.SStotal = norm(Y-mean(Y)).^2;
stats.SSeffect = norm(Yhat-mean(Yhat)).^2;
stats.SSerror  = norm(Res).^2;
stats.df = rank(X)-1;
stats.dferror = length(Y) - rank(X);
stats.R2 = 1-(stats.SSerror / stats.SStotal);
stats.F = (stats.SSeffect / stats.df) / (stats.SSerror / stats.dferror);
stats.p = 1 - fcdf(stats.F,stats.df,stats.dferror);
stats.Var_e=(Res'*Res)/(stats.dferror);  % -> Random effects covariance parameters:

stats_each=[];
% compute partial and semipartial correlations
if nargout > 2
  % ------------------------------------
  % semi-partial correlation coefficient
  % ------------------------------------
  % let's think of the model and what it means it terms of geometry.
  % the data Y can be described as a vector and a point in R_20
  % a space with 20 dimensions. We then establish a model X with 5
  % regressors; that is we look for a combination of these 5 vectors which
  % will get as close as possible to Y. To find the actual contribution of x1
  % to the data for this model one needs to look at how much x1 explains
  % to the total variance, ie we want to compare the R2 between the full and
  % a reduced model without x1 - the difference will be how much x1 explains in Y.
  
  for j=1:size(X,2)
    puse = ones(1,size(X,2));
    puse(j) = 0; puse = find(puse);
    
    Xreduced    = X(:,puse); % reduced model all minus 1st regressor
    Breduced    = inv(Xreduced'*Xreduced)*Xreduced'*Y;
    Yhatreduced = Xreduced*Breduced;
    Resreduced  = Y - Yhatreduced;
    
    dfreduced       = rank(Xreduced) -1 ;
    dferrorreduced = length(Y) - dfreduced - 1;
    SSeffectreduced = norm(Yhatreduced-mean(Yhatreduced)).^2;
    SSerrorreduced  = norm(Resreduced-mean(Resreduced)).^2;
    R2reduced       = SSeffectreduced / stats.SStotal;
    Freduced       = (SSeffectreduced / dfreduced) / (SSerrorreduced / dferrorreduced);
    preduced       = 1 - fcdf(Freduced,dfreduced,dferrorreduced);

    Semi_Partial_corr_coef = stats.R2 - R2reduced;
    dfe_semi_partial_coef  = stats.df - dfreduced;
    F_semi_partail_coef    = (Semi_Partial_corr_coef*stats.dferror) / ...  % variance explained by x1
      ((1-stats.R2)*dfe_semi_partial_coef); % unexplained variance overall
    p_semi_partial_coef    = 1 - fcdf(Semi_Partial_corr_coef, stats.df, dfe_semi_partial_coef); % note df is from the full model
    
    stats_each.semiP(j,:) = [Semi_Partial_corr_coef,F_semi_partail_coef,p_semi_partial_coef];
    
    % --------------------------------
    % partial correlation coefficient
    % --------------------------------
    % As we shall see below, this is easily obtained using projections - but
    % for now let just think about what we want to measure. We are interested
    % in knowing the correlation between y and x1 controlling for the effect of
    % the other xs, that is removing the effect of other xs. Compred to
    % semi-parital coef we also want to remove the effect of xs on y or if you
    % prefer we want to compute how much of x1 we need to get to the point
    % defined by the xs which is the closest to Y
    
    
    x = X(:,j);
    Btmp = inv(Xreduced'*Xreduced)*Xreduced'*x;
    xhat = Xreduced*Btmp;
    Resx = x - xhat;
    
    % the correlation between Resreduced and Resx is the partial coef
    stats_each.P(j,:) = corr(Resx,Resreduced);
  end
end

contrast =[];
% compute contrast
if nargin ==3
  if ~iscell(C)
    C = {C};
  end
  for i=1:length(C)
    Cuse = C{i}(:);
    % C*B1 is distributed around the true C*B with known variance
    contrast.t_stat(i)=Cuse'*B/sqrt(stats.Var_e*Cuse'*pinv(X'*X)*Cuse);
  end
  contrast.pval=2*tcdf(-abs(contrast.t_stat),stats.dferror); %CHECK DF HERE!!!!!!!!
end






return;

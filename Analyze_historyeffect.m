
% Analysis of VE and VAE biases across multiple datasets (exps 1-10 in the manuscript)
% Experiments consisted of sequences of AV-A trials, sometimes interrupred by V trials (AV-A-V)
%
% Focus here is the dependence of VE/VAE on the multisensory discrepancy (V-A position) 
% across multiple AV trials.

clear; close all

% load data
load('DataExp1to10.mat')

% Dataset{}  contains the data for each experiment 
% Dataset{1}(trial,:) contains as columns 1-8 the following for each AV-A trial pair:
% c1: AV trial V stimulus location
% c2: AV trial A stimulus location
% c3: A trial A stimulus location
% c4: dVA: V - A position in the AV trial (a.k.a, audio-visual discrepancy)
% c5: AV trial response
% c6: VE bias (defined as A pos minus response in the AV trial)
% c7: A trial response
% c8: VAE (defined as A pos minus the mean response for this A location across all A trials)

%------------------------------------------------------------------------

VAE_all=[];
VE_all=[];
c=1; % counter across exps and s
for d=1:length(Dataset) % experiments
  n = length(Dataset{d});
  for s=1:n % participants
    tmp = Dataset{d}{s}; % local copy
    
    % ---------------------------------------------------------------------------
    % analysis 1: Select sequences of 3 AV-A pairs (a,b,c) not interrupted with a trial with discrepancy (DVA) of zero.
    % compute proportional bias contingent on relative sign of previous discrepancies
    % ---------------------------------------------------------------------------

    % create the following conditions
    % {'c','b=c','a=b=c','a=b~=c'};
    for cond=1:4
      data_local{cond}=[];
    end
    
    for t = 3:size(tmp,1) % triplets of trial pairs
      dvaseq = tmp(t+[-2:0],4); % DVA
      vae = tmp(t,8); % biases
      ve = tmp(t,6); 
      
      tmpdata = [ve,vae,dvaseq(3)];
      if ~sum(dvaseq==0)   % only current trial
        conduse = 1;
        data_local{conduse} = cat(1,data_local{conduse},tmpdata);
        if sign(dvaseq(3))==sign(dvaseq(2))          % n-1 same as n
          conduse = 2;
          data_local{conduse} = cat(1,data_local{conduse},tmpdata);
        end
        if abs(sum(sign(dvaseq)))==3          % all the same
          conduse = 3;
          data_local{conduse} = cat(1,data_local{conduse},tmpdata);
        elseif  sign(dvaseq(3))~=sign(dvaseq(2)) && sign(dvaseq(3))~=sign(dvaseq(1))
          % n-1 and n-2 the same and different than n
          conduse = 4;
          data_local{conduse} = cat(1,data_local{conduse},tmpdata);
        end
      end
    end % t
    
    % derive proportionnal biases relative to the DVA on trial C
    for cond=1:4
      m = data_local{cond}(:,3);
      VE_all(c,cond) = mean(data_local{cond}(:,1)./m);
      VAE_all(c,cond) = mean(data_local{cond}(:,2)./m);
    end
    
    % -------------------------------------------------------------------------
    % collect single trial data incl. all trials (Also zero DVA) for regression modelling
    X=[];
    for t = 3:size(tmp,1)
      dvaseq = tmp(t+[-2:0],4);
      vae = tmp(t,8);
      ve = tmp(t,6);
      X = cat(1,X,[ve,vae,dvaseq']);
    end
    
    
    % -------------------------------------------------------------------------
    % analysis 2: partial regression with n-trials into past and comparison of model prediction power
    % -------------------------------------------------------------------------
    % setup model
    model = X(:,[3:5]);
    model(:,end+1) = 1;
    % partial regression beta's
    [beta] = ck_stat_glm(model,X(:,2));
    Partial_beta_vae(c,:) = beta(1:end-1);
    
    [beta]  = ck_stat_glm(model,X(:,1));
    Partial_beta_ve(c,:) = beta(1:end-1);
    
    % -------------------------------------------------------------------------
    % analysis 3: model comparison of models with and without past 
    % -------------------------------------------------------------------------
    M={};
    M{1} =[3,4]; % trial c only
    M{2} =[2,3,4]; % trial b,c only
    M{3} =[1:4]; % trial a,b,c only
    
    for p=1:3 % model
      % call regression model and get LL
      [b,Stats] = ck_stat_regress(X(:,1),model(:,M{p}));
      LL = Stats.model(end);
      AIC = 2*length(M{p})-2*(LL);
      VEmodel(c,p,:) = [AIC];
      [b,Stats] = ck_stat_regress(X(:,2),model(:,M{p}));
      LL = Stats.model(end);
      AIC = 2*length(M{p})-2*(LL);
      VAEmodel(c,p,:) = [AIC];
    end
    c=c+1;
  end %sub
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10);clf;
Colorvector(1,:) = [250,150,0]/250; %
Colorvector(2,:) = [250,150,0]/250; %
Colorvector(3,:) = [250,150,0]/250; %
Colorvector(4,:) = [0,150,250]/250; %

%-----------------------------------------------------------------------
% Display proportional biases - Analysis 1
%-----------------------------------------------------------------------
Label = {'c','b=c','a=b=c','a=b~=c'};

% VE
subplot(2,2,1); hold on
line([0.5 4.5],[0 0],'LineWidth',2,'color',[0.5 0.5 0.5]);
ckboxplotcompact(VE_all,[1:size(VE_all,2)],1,Colorvector);
axis([0.5 4.5 -0.6 1.2]);
set(gca,'XTick',[1:4],'XTickLabel',Label,'YTick',[-0.6:0.3:1.2]);
ylabel('Proportional bias');
text(2.5,1.4,'VE','FontSize',15,'FontWeight','Bold')
xtickangle(45)
grid on
set(gca,'FontSize',13);

% VAE
subplot(2,2,2); hold on
line([0.5 4.5],[0 0],'LineWidth',2,'color',[0.5 0.5 0.5]);
ckboxplotcompact(VAE_all,[1:size(VAE_all,2)],1,Colorvector);
axis([0.5 4.5 -0.2 0.6]);
set(gca,'XTick',[1:4],'XTickLabel',Label,'YTick',[-0.2:0.2:0.6]);
ylabel('Proportional bias');
text(2.4,0.69,'VAE','FontSize',15,'FontWeight','Bold')
xtickangle(45)
grid on
set(gca,'FontSize',13);


% Median and bootstrap p-values 
for t=1:4
  x = VE_all(:,t);
  [Lo,Up,p]=ck_stat_bootCI_percentile(x,'median',0.01,8000,0);
  CIs_VE(t,:) =[median(x) Lo Up];
  x = VAE_all(:,t);
  [Lo,Up,p]=ck_stat_bootCI_percentile(x,'median',0.01,8000,0);
  CIs_VAE(t,:) =[median(x) Lo Up];
end
fprintf('VE prop bias (Median and CI): \n');
for t=1:4
  fprintf('%1.3f [%1.3f, %1.3f]  \n',CIs_VE(t,1),CIs_VE(t,2),CIs_VE(t,3));
end
fprintf('VAE prop bias (Median and CI): \n');
for t=1:4
  fprintf('%1.3f [%1.3f, %1.3f]  \n',CIs_VAE(t,1),CIs_VAE(t,2),CIs_VAE(t,3));
end



%-----------------------------------------------------------------------
% Display partial betas - Analysis 2
%-----------------------------------------------------------------------

lgd = {'a','b','c'};
% VE
subplot(2,2,3); hold on
line([0.5 3.5],[0 0],'LineWidth',2,'color',[0.5 0.5 0.5]);
ckboxplotcompact(Partial_beta_ve(:,[1:3]),[1:3],1);
ylabel('Partial beta'); xlabel('trial');
axis([0.5 3.5 -0.4 1]);
set(gca,'XTick',[1:3],'XTickLabel',lgd,'YTick',[-0.4:0.2:1]);
grid on
set(gca,'FontSize',13);
% VAE
subplot(2,2,4); hold on
line([0.5 3.5],[0 0],'LineWidth',2,'color',[0.5 0.5 0.5]);
ckboxplotcompact(Partial_beta_vae(:,[1:3]),[1:3],1);
ylabel('Partial beta');xlabel('trial');
axis([0.5 3.5 -0.11 0.25]);
set(gca,'XTick',[1:3],'XTickLabel',lgd,'YTick',[-0.1:0.1:0.2]);
grid on
set(gca,'FontSize',13);


% Median and bootstrap p-values 
for t=1:3
  x = Partial_beta_ve(:,t);
  [Lo,Up,p]=ck_stat_bootCI_percentile(x,'median',0.01,8000,0);
  CIs_VEB(t,:) =[median(x) Lo Up];
  x = Partial_beta_vae(:,t);
  [Lo,Up,p]=ck_stat_bootCI_percentile(x,'median',0.01,8000,0);
  CIs_VAEB(t,:) =[median(x) Lo Up];
end
fprintf('VE partial betas (Median and CI):  \n');
for t=1:3
  fprintf('%1.3f [%1.3f, %1.3f]  \n',CIs_VEB(t,1),CIs_VEB(t,2),CIs_VEB(t,3));
end
fprintf('VAE partial betas (Median and CI): \n');
for t=1:3
  fprintf('%1.3f [%1.3f, %1.3f]  \n',CIs_VAEB(t,1),CIs_VAEB(t,2),CIs_VAEB(t,3));
end


% ----------------------------------
% we compute the  ratio of beta's
% to be more robust, we compute the ratio of the medians, and use boostrap CIs

x = Partial_beta_ve(:,[2,3]);
ratio = median(x(:,1))./median(x(:,2));
for b=1:8000
  rr = ceil(rand(1,205)*205);
  x = Partial_beta_ve(rr,[2,3]);
  ratioboot(b) = median(x(:,1))./median(x(:,2));
end
[Ci] = prctile(ratioboot,[0.5,99.5]);

fprintf('VE ratio of median partial betas (CI):  ');
fprintf('%1.3f [%1.3f, %1.3f]  \n',ratio,Ci(1),Ci(2))
 
x = Partial_beta_vae(:,[2,3]);
ratio = median(x(:,1))./median(x(:,2)); 
for b=1:8000
  rr = ceil(rand(1,205)*205);
  x = Partial_beta_vae(rr,[2,3]);
  ratioboot(b) = median(x(:,1))./median(x(:,2));
end
[Ci] = prctile(ratioboot,[0.5,99.5]);

fprintf('VAE ratio of median partial betas (CI):  ');
fprintf('%1.3f [%1.3f, %1.3f]  \n',ratio,Ci(1),Ci(2))
 

% direct comparison of relative growth 
x = Partial_beta_ve(:,[2,3]);
ratiove = median(x(:,1))./median(x(:,2));
x = Partial_beta_vae(:,[2,3]);
ratiovae = median(x(:,1))./median(x(:,2));
delta = ratiovae./ratiove;
for b=1:8000
  rr = ceil(rand(1,205)*205);
  x = Partial_beta_ve(rr,[2,3]);
  ratiove = median(x(:,1))./median(x(:,2));
  x = Partial_beta_vae(rr,[2,3]);
  ratiovae = median(x(:,1))./median(x(:,2));
  deltaboot(b) = ratiovae./ratiove;
end
[Ci] = prctile(deltaboot,[0.5,99.5]);



%-----------------------------------------------------------------------
% Comparison of model fit - Analysis 3
%-----------------------------------------------------------------------
%VE
fprintf('VE Model comparison:  \n');
aic =  sum(sq(VEmodel(:,:,1)));
aic = aic-min(aic)
WAIC_VE = exp(-0.5*aic)./(sum(exp(-0.5*aic)))
%VAE
fprintf('VAE Model comparison:  \n');
aic =  sum(sq(VAEmodel(:,:,1)));
aic = aic-min(aic)
WAIC_VAE = exp(-0.5*aic)./(sum(exp(-0.5*aic)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure with average biases for each dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);clf;
cmap = jet(10);
cmap(8,:) = [0.8 0.8 0];

for d=1:length(Dataset)
  n = length(Dataset{d});
  
  Mean_biases =[];
  for s=1:n %participant
    tmp = Dataset{d}{s};
    dva = unique(tmp(:,4));
    for dv=1:length(dva)
      j = find(tmp(:,4)==dva(dv));
      Mean_biases(s,dv,:) = mean(tmp(j,[6,8]));
    end
  end
  
  subplot(1,2,1); hold on
  h = errorbar(dva,mean(Mean_biases(:,:,1)),sem(Mean_biases(:,:,1)));
  set(h,'Color',cmap(d,:));
  set(h,'LineWidth',1.5);
  xlabel('DVA [deg]');
  ylabel('Bias [deg]')
  axis([-55 55 -20 20]);
  set(gca,'XTick',[-50:25:50],'YTick',[-20:5:20]);
  text(0,22,'VE','FontSize',14,'FontWeight','Bold')
  
  
  subplot(1,2,2); hold on
  h = errorbar(dva,mean(Mean_biases(:,:,2)),sem(Mean_biases(:,:,2)));
  set( h,'LineWidth',1.5);
  set(h,'Color',cmap(d,:));
  xlabel('DVA [deg]');
  ylabel('Bias [deg]')
  axis([-55 55 -4 4]);
  set(gca,'XTick',[-50:25:50],'YTick',[-4:2:4]);
  text(-2,4.4,'VAE','FontSize',14,'FontWeight','Bold')
  text(-40,3.9-d*0.3,sprintf('Exp %d',d),'Color',cmap(d,:),'FontSize',8);
  
end





%% --------------------------------------------------------------------------------------------------
%% Additional analysis on the cumulative growth after removing the expected cumulative influence of a VAE
% cumulative growth of VE
d1 = VE_all(:,3)-VE_all(:,2); % 3-2
[Lo,Up,p]= ck_stat_bootCI_percentile(d1,'median',0.01,8000,0);
CI_growth(1,:) =[median(d1) Lo Up];

fprintf('VE cumulative growth uncorrected (Median and CI):  \n');
fprintf('%1.3f [%1.3f, %1.3f]  \n',CI_growth(1,1),CI_growth(1,2),CI_growth(1,3));


% cumulative growth corrected for VAE
a = VE_all(:,3)-VAE_all(:,2);
b = VE_all(:,2)-VAE_all(:,1);
d2 = a-b; 
[Lo,Up,p]=ck_stat_bootCI_percentile(d2,'median',0.01,8000,0);
CI_growth(2,:) =[median(d2) Lo Up];
fprintf('VE cumulative growth corrected (Median and CI):  \n');
fprintf('%1.3f [%1.3f, %1.3f]  \n',CI_growth(2,1),CI_growth(2,2),CI_growth(2,3));



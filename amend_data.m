% 2022 November Hame Park
% Upon review - 
% VE and VAE computed as R - mean(A) and R - Apos respectively, inserted in
% the beginning of each subject dataset. Now all dataset has two types of
% VE and VAE. I also added the mean A trial response in the third column. 
% The below dataset description has also been amended accordingly.


% Analysis of VE and VAE biases across multiple datasets (exps 1-10 in the manuscript)
% Experiments consisted of sequences of AV-A trials, sometimes interrupred by V trials (AV-A-V)
%
% Focus here is the dependence of VE/VAE on the multisensory discrepancy (V-A position) 
% across multiple AV trials.



clear; close all

% Dataset{}  contains the data for each experiment 
% Dataset{1}(trial,:) contains as columns 1-10 the following for each AV-A trial pair:
% -amended 11.2022

% new column indices:
% c1: VE bias (defined as AV trial response minus mean of all A trial response for each A stimulus location)
% c2: VAE bias (defined as A trial response minus A trial stimulus
% location)
% c3: AV trial V stimulus location
% c4: AV trial A stimulus location
% c5: A trial A stimulus location
% c6: dVA: V - A position in the AV trial (a.k.a, audio-visual discrepancy)
% c7: AV trial response
% c8: VE bias (defined as AV trial response minus AV trial A stimulus location)
% c9: A trial response
% c10: VAE bias (defined as A trial response minus mean of all A trial response for each A stimulus location)


% load data
load('DataExp1to10.mat')
%------------------------------------------------------------------------
% here still use the old column indices
% c1: AV trial V stimulus location
% c2: AV trial A stimulus location
% c3: A trial A stimulus location
% c4: dVA: V - A position in the AV trial (a.k.a, audio-visual discrepancy)
% c5: AV trial response
% c6: VE bias (defined as AV trial response minus AV trial A stimulus location)
% c7: A trial response
% c8: VAE bias (defined as A trial response minus mean of all A trial response for each A stimulus location)

for d = 1:length(Dataset)
    % get new VAE: A trial response minus A stimulus location. 
    newVAE = cellfun(@(x) x(:, 7) - x(:, 3), Dataset{d}, 'UniformOutput', false);
    % get new VE: AV trial response minus mean(A stimulus location) per
    % location.
    % A stimulus positions
    aloc = unique(Dataset{d}{1}(:, 3)); % this is always in ascending order
    % get mean A trial response for each A stimulus position
    mua = cellfun(@(x) (permute(repmat(x(:, 7), 1, 5), [2 1])*(x(:, 3)==aloc'))./repmat(sum(x(:, 3)==aloc', 1), 5, 1), ... 
        Dataset{d}, 'UniformOutput', false); % this may not be monotonously increasing
    mua = cellfun(@(x) x(1, :), mua, 'UniformOutput', false);
    % make it into a matrix the size of 5 x trials
    mua = cellfun(@(x, y) repmat(x', 1, size(y, 1)), mua, Dataset{d}, 'UniformOutput', false);
    % assign the mean A trial response into corresponding A stimulus in AV
    % trial
    muApos = cellfun(@(x, y) x(y(:, 2)'==aloc), mua, Dataset{d}, 'UniformOutput', false);
    % get new VE: AV trial response minus mean of all A trial response per
    % position
    newVE = cellfun(@(x, y) x(:, 5)-y, Dataset{d}, muApos, 'UniformOutput', false);
    % put new arrays in the beginning of datasets
    Dataset{d} = cellfun(@(x, y, z) cat(2, y, z, x), Dataset{d}, newVE, newVAE, 'UniformOutput', false);
end

save('DataExp1to10_new.mat', 'Dataset')

%% plot new and old VE and VAE
load('DataExp1to10_new.mat', 'Dataset')
% from here use the new column indices described above.
vis = figure;
ii = 1;
for d = 1:length(Dataset)
    dva = unique(Dataset{d}{1}(:, 6));
    % plot VE, w.r.t stim location (original)
    vemu = permute(cellfun(@(x) bsxfun(@(x, y) (y'*(x==dva'))./sum(x==dva'), x(:, 6), x(:, 8)), Dataset{d}, 'UniformOutput', false), [2 1]);
    subplot(10, 4, ii), plot(dva, cell2mat(vemu)), hold on; plot(dva, mean(cell2mat(vemu), 1), '-ok', 'LineWidth', 2)
    if ii == 1, title('VE (R - stim.loc)'), end
    ii = ii+1;
    % plot VE, w.r.t mean A response
    vemu = permute(cellfun(@(x) bsxfun(@(x, y) (y'*(x==dva'))./sum(x==dva'), x(:, 6), x(:, 1)), Dataset{d}, 'UniformOutput', false), [2 1]);
    subplot(10, 4, ii), plot(dva, cell2mat(vemu)), hold on; plot(dva, mean(cell2mat(vemu), 1), '-ok', 'LineWidth', 2)
    if ii == 2, title('VE (R - mean(Astim.loc))'), end
    ii = ii+1;
    % plot VAE, w.r.t stim location
    vaemu = permute(cellfun(@(x) bsxfun(@(x, y) (y'*(x==dva'))./sum(x==dva'), x(:, 6), x(:, 2)), Dataset{d}, 'UniformOutput', false), [2 1]);
    subplot(10, 4, ii), plot(dva, cell2mat(vaemu)), hold on; plot(dva, mean(cell2mat(vaemu), 1), '-ok', 'LineWidth', 2)
    if ii == 3, title('VAE (R - stim.loc)'), end
    ii = ii+1;
    % plot VAE, w.r.t mean A response (original)
    vaemu = permute(cellfun(@(x) bsxfun(@(x, y) (y'*(x==dva'))./sum(x==dva'), x(:, 6), x(:, 10)), Dataset{d}, 'UniformOutput', false), [2 1]);
    subplot(10, 4, ii), plot(dva, cell2mat(vaemu)), hold on; plot(dva, mean(cell2mat(vaemu), 1), '-ok', 'LineWidth', 2)
    if ii == 4, title('VAE (R - mean(Astim.loc))'), end
    ii = ii+1;
end


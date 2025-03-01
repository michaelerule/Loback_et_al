function [Output] = StatLin_Concat_match(varargin)
%__________________________________________________________________________
% Description: This version matches the number of trials used for the
%              train and test set for each day in the generated 
%              concatenated "dataset" to match the average for each
%              individual day (generated by StatLin_AvgMult_kCV.m). 
% Written 18 December 2018 by AL
% Note: I verified on 18/12/2018 that the random permutations are indeed 
%       different between the different V.nperms iterations.
%__________________________________________________________________________

%% -- Initializations: -- 
V.mouse        = 4;
V.BaseSessions = [6];                  
V.behavvars    = [3 6 4];     %3=forward position (units: m); 6=forward speed
V.trial_type   = 'all';       %2 = black right-turn trial;
                              %3 = white left-turn trial;
                              %'all' specifies both trial types
V.preproc        = 1;         %1=preprocess data (linear model); 0=don't preprocess (for affine model) 
V.save           = 1;         %1=save output; 0=don't save 
V.viz            = 0;         %1=visualize results; 0=don't plot
V.N              = 'all';     %'all'=use all; else use float # specified
V.nperms         = 10;        %# of permutations to run for concat
V                = parseargs(V, varargin{:});   

LoadDir = '/home/mrule/Workspace2/PPC_data/';
fprintf(1,'Loading datasets from %s,\n',LoadDir);
fprintf(1,'please re-define LoadDir to match the path to the data files on your system.\n');
fprintf(1,'Datasets are available from the Harvey lab.\n');

rng(4); %set rng
MAE_cache = NaN(length(V.behavvars),V.nperms);

%% -- Load Data Regarding # of Trials for Train & Test Set: -- 

loadfile = ['m0' num2str(V.mouse) '_bs6_AvgBestPerfkCV_behavs23564.mat'];
O        = load([pwd '/Results/' loadfile],'O'); O=O.O; 
V.Diffs  = O.Diffs; %ensure use correct matching  
Output.BaseSessions = V.BaseSessions; %cache 
Output.Diffs        = V.Diffs;        %cache 

%% -- Compute Sessions: --
TestSessions = NaN(length(V.BaseSessions),length(V.Diffs)); 
for d_idx = 1:length(V.BaseSessions)
    TestSessions(d_idx,:) = V.BaseSessions(d_idx)+V.Diffs; 
end 
allsessions   = unique([V.BaseSessions TestSessions(:)']);
V.allsessions = allsessions; %base session values

V.avg_ntrain_trials = round(O.avg_ntrain_trials/length(allsessions)); %avg. trials to include in train set *per session*
V.avg_ntest_trials  = round(O.avg_ntest_trials /length(allsessions)); %avg. trials to include in test set *per session*

%% -- Data Preparation for Each Random Permutation: -- 
%goodtrials_perm_cell = cell(1,V.nperms); 
for rp=1:V.nperms
    
    [Z_train,Z_test,X_train,X_test,~] = PrepareData_Concat_match(V); 
    %goodtrials_perm_cell{rp} = gtrials; %to verify different permutations
    if V.preproc
       [Z_train,Z_test,X_train,X_test,...
        ~,X_train_mean] = PreProcess(Z_train,Z_test,X_train,X_test); 
    else 
       X_train_mean   = mean(X_train,1);
    end 

    % -- Fit Linear Model via LS Regression:
    [M_star,~] = FitM_UC(Z_train',X_train'); 

    % -- Evaluate Model Performance on Held-Out Test Data: -- 
    MAE_cache(:,rp) = getDecodePerf(Z_test, X_test, M_star); 
    
    if V.viz
       visualise_perf(Z_test,X_test,M_star,X_train_mean); 
    end 
end %over random permutations

Output.MAE_cache = MAE_cache; %save

%% -- Perform Mann-Whitney U Test: --
%forward position 
[p,h,stats]       = ranksum(O.avg_bestdecodeperf_CV(1,:),MAE_cache(1,:)); 
Output.p_fpos     = p;
Output.h_fpos     = h;
Output.stats_fpos = stats; 
%forward velocity
[p,h,stats]       = ranksum(O.avg_bestdecodeperf_CV(2,:),MAE_cache(2,:)); 
Output.p_fvel     = p;
Output.h_fvel     = h;
Output.stats_fvel = stats; 
%view angle
[p,h,stats]       = ranksum(O.avg_bestdecodeperf_CV(3,:),MAE_cache(3,:)); 
Output.p_va       = p;
Output.h_va       = h;
Output.stats_va   = stats;

%% -- Save: -- 
if V.save
    savename = ['m0' num2str(V.mouse) 'bs' num2str(V.BaseSessions) '_ConcatMatch_behavs364.mat'];
    save([pwd '/Results/' savename], 'Output');
end

end %main

%% -- Internal Functions: -- 
% (1)______________________________________________________________________
function [Z_train,Z_test,X_train,X_test,goodtrials_cell] = PrepareData_Concat_match(S)
    %{
    MR20190722

    Prepares data for concateated analyses. 

    Parameters
    ----------
    S.allsessions: vector
        List of sessions for concatenated spans of days. 
    S.N:
        Whether or not to shuffle neurons
    S.avg_ntrain_trials: int
        Number of trials use for training
    S.avg_ntest_trials: int
        Number of trials used for testing
    %}
    LoadDir = ['../Data/m0' num2str(S.mouse) '/'];
    bd_vec  = S.allsessions; %base session values

    % -- Find consistent set of good neurons over days:
    goodneurons_pre = []; 
    for b_idx = 1:length(bd_vec)
        bd    = bd_vec(b_idx); 
        Datafile = sprintf('%sm%02d_s%02d.mat',LoadDir,s.mouse,bd); %MR20190722
        session_obj = load(Datafile, 'session_obj'); 
        session_obj = session_obj.session_obj; 
        gn_pre      = find(~isnan(session_obj.confidenceLabel)); 
        if isempty(goodneurons_pre)
            goodneurons_pre = gn_pre; 
        else 
            goodneurons_pre = intersect(goodneurons_pre, gn_pre); 
        end 
    end 
    if ischar(S.N)
        goodneurons = goodneurons_pre;
    else
        goodneurons = goodneurons_pre(randperm(S.N)); 
    end 
    
    % -- Data Prep for Concatenated Days (into Single Dataset):
    Z_train = [];
    Z_test  = [];
    X_train = [];
    X_test  = [];
    total_ntrials_cell = NaN(1,length(bd_vec));  %init cache
    goodtrials_cell    = cell(1,length(bd_vec)); %init cache
    
    for b_idx = 1:length(bd_vec) %over sessions to concatenate
        bd = bd_vec(b_idx); 

        % -- Find Time Bins Comprising Correct Trials of Specified Type: 
        O_tbininfo = ComputeActivity_SpecifyTrials('mouse',S.mouse,'session',bd); 
                 
        if ischar(S.trial_type)
           trials_to_include = O_tbininfo.trialtimes_correct_all;
        elseif S.trial_type==2
           trials_to_include = O_tbininfo.trialtimes_correct_BR;
        elseif S.trial_type==3
           trials_to_include = O_tbininfo.trialtimes_correct_WL; 
        end 

        % -- Partition total trials into training & test set:
        total_ntrials_cell(b_idx) = size(trials_to_include,1); 
        goodtrials_cell{b_idx}    = trials_to_include(randperm(size(trials_to_include,1)),:);
    end 

    % ** Ensure proportion of each day is consistent for train & test: **
    unif_triallen    = min(total_ntrials_cell); 

    % MR20190722 goodtrials_cell2 stores the start and ending
    % timepoints of each good trial from each session.
    goodtrials_cell2 = cellfun(@(y) y(1:unif_triallen,:), goodtrials_cell, 'UniformOutput', 0);
    ntrials_train    = S.avg_ntrain_trials; %matches avg # of trials per day for train set
    ntrials_test     = S.avg_ntest_trials;  %matches avg # of trials per day for test set

    for b_idx = 1:length(bd_vec) %over sessions to concatenate
        bd = bd_vec(b_idx); 
        % -- Prepare Data: --
        if bd<10
           Datafile = ['m0' num2str(S.mouse) '_s0' num2str(bd) '.mat']; 
        else
           Datafile = ['m0' num2str(S.mouse) '_s' num2str(bd) '.mat']; 
        end 
        session_obj = load([LoadDir Datafile], 'session_obj'); 
        session_obj = session_obj.session_obj; 

        % Compute time bin indices to include:
        bins_to_include_train = [];
        bins_to_include_test  = [];
        for bi=1:ntrials_train
            addon = goodtrials_cell2{b_idx}(bi,1):goodtrials_cell2{b_idx}(bi,2); 
            bins_to_include_train = [bins_to_include_train; addon(:)]; %#ok
        end 
        for bi=1:ntrials_test
            addon = goodtrials_cell2{b_idx}(bi+ntrials_train,1):goodtrials_cell2{b_idx}(bi+ntrials_train,2);
            bins_to_include_test  = [bins_to_include_test; addon(:)];  %#ok
        end 

        % Prepare training set (all data) of Ca2+ dF/F values:
        Z = session_obj.timeSeries.calcium.data{1,1}; %{2,1} = use deconvolved
        Z_train = [Z_train; Z(goodneurons,bins_to_include_train)']; %#ok
        Z_test  = [Z_test; Z(goodneurons,bins_to_include_test)'];   %#ok 

        % Prepare train set of kinematic values: 
        Y = session_obj.timeSeries.virmen.data(S.behavvars,:);
        if sum(isnan(Y))>0, error('NaN entries/n'); end 
        X_train = [X_train; Y(:,bins_to_include_train)'];           %#ok
        X_test  = [X_test; Y(:,bins_to_include_test)'];             %#ok
    end 
end 
%__________________________________________________________________________
% (2)______________________________________________________________________
function [Z_train,Z_test,X_train,X_test,...
          Z_train_mean,X_train_mean] = PreProcess(Z_train,Z_test,...
                                                  X_train,X_test)
    % Pre-processing: 
    %Z-score:
    Z_train_mean = mean(Z_train,1); 
    %Z_train_std = std(Z_train,0,1); 
    Z_train      = (Z_train - Z_train_mean); %./Z_train_std; 
    Z_test       = (Z_test - Z_train_mean);  %./Z_train_std;

    %0-center:
    X_train_mean = mean(X_train,1);
    X_train      = (X_train - X_train_mean);
    X_test       = (X_test - X_train_mean); 
end 
%__________________________________________________________________________
% (3)______________________________________________________________________
function [M_star,xi] = FitM_UC(Z,X)
    M_star = X*Z'*inv(Z*Z'); %#ok
    xi     = (norm(X - (M_star*Z),'fro'))^2; %RSS 
end 
%__________________________________________________________________________
% (4)______________________________________________________________________
function [mae_kin] = getDecodePerf(Z_test,X_test,M_b)
    [ntest_samples,p] = size(X_test);
    % Initialize Caches: 
    xk_preds = zeros(p,ntest_samples);
    evec_kin = zeros(p,ntest_samples);
    for k=1:ntest_samples
        zk              = Z_test(k,:);
        xk              = X_test(k,:);    %actual kinematic value
        xk_hat          = M_b*zk(:);      %model pred for kinematics
        evec_kin(:,k)   = xk(:) - xk_hat(:);
        xk_preds(:,k)   = xk_hat;
    end 

    mae_kin = mean(abs(evec_kin),2);
end
%__________________________________________________________________________
% (5)______________________________________________________________________ 
function [] = visualise_perf(Z_test,X_test,M_star,X_train_mean)
    
    % -- Initializations: --
    Labels = {'Lateral Position (m)', 'Forward Position (m)', ...
              'Lateral Speed (m/s)', 'Forward Speed (m/s)', ...
              'View Angle ($^{\circ}$)'};
    figure;
    [ntest_samples,p] = size(X_test);
    % Initialize Caches: 
    xk_preds = zeros(p,ntest_samples);
    evec_kin = zeros(p,ntest_samples);
    for k=1:ntest_samples
        zk              = Z_test(k,:);
        xk              = X_test(k,:);       %actual kinematic value
        xk_hat          = M_star*zk(:);      %model pred for kinematics
        evec_kin(:,k)   = xk(:) - xk_hat(:);
        xk_preds(:,k)   = xk_hat;
    end
    
    % -- Main computations: -- 
    for idx = 1:p
        subplot(p,1,idx);
        plot1 = plot(X_test(:,idx)+X_train_mean(idx),'k','linewidth',2); 
        plot1.Color(4) = 0.3;
        hold on;
        plot2 = plot(xk_preds(idx,:)+X_train_mean(idx),'b','linewidth',2);
        plot2.Color(4) = 0.5;
        xlim([0 1000]); ax=gca; ax.XTick = [0 1000]; %axis square; 
        ylabel(Labels{idx},'Interpreter','Latex');
        if idx==p
            xlabel('Time Bin'); 
        end 
    end 
    % -- Format font etc. -- 
    a=findobj(gcf); alltext=findall(a,'Type','text');
    allaxes=findall(a,'Type','axes');
    set(allaxes,'FontName','Helvetica','FontWeight','Normal','FontSize',16);
    set(alltext,'FontName','Helvetica','FontWeight','Normal','FontSize',20);
end 
%__________________________________________________________________________

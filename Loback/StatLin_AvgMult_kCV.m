function [O] = StatLin_AvgMult_kCV(varargin)
%__________________________________________________________________________
% Description: For part (1ii,b) of the paper analyses, i.e. aims to
%              address whether the decoding performance for the 
%              concatenated dataset is slightly worse vs. the 
%              average of the best linear decoder performance on each
%              individual day. 
%              This version performs k-fold CV. 
% Written 18 Dec 2018 by AL
%__________________________________________________________________________

%% -- Initializations: -- 
V.mouse        = 3;     
V.BaseSessions = [7];           %#ok
V.Diffs        = [-5 -1 0 1 4];             
V.behavvars    = [3 6 4];       %3=forward position (units: m); 6=forward speed
V.trial_type   = 'all';         %2=black right-turn trial;
                                %3=white left-turn trial;
                                %'all' specifies both trial types
V.kfold          = 10;          %# of folds (k) for cross-validation 
V.preproc        = 1;           %1=preprocess data (linear model); 0=don't preprocess (for affine model) 
V.save           = 1;           %1=save output; 0=don't save 
V.N              = 'all';       %'all'=use all; else use float # specified
V.deconv         = 0;           %1=use deconvolved signals; 0=use raw deltaF/F
V                = parseargs(V, varargin{:});  


LoadDir = '/home/mrule/Workspace2/PPC_data/';
fprintf(1,'Loading datasets from %s,\n',LoadDir);
fprintf(1,'please re-define LoadDir to match the path to the data files on your system.\n');
fprintf(1,'Datasets are available from the Harvey lab.\n');

O.BaseSessions = V.BaseSessions;
O.Diffs        = V.Diffs; 
                
%% -- Compute Best Linear Decoder Perf. on Each Day & Avg.: -- 
[allsessions,goodneurons] = getsessions(V); 
V.goodneurons             = goodneurons; 
p                         = length(V.behavvars); 
MAE_cache                 = NaN(p,V.kfold,length(allsessions)); %init
M_cache                   = cell(1,length(allsessions)); %init cache
Ztrain_cache              = M_cache; %init 
Xtrain_cache              = M_cache; %init
Ztest_cache               = M_cache; %init
Xtest_cache               = M_cache; %init
xi_cache                  = NaN(1,length(allsessions));  %init cache

for s_idx = 1:length(allsessions)
    sess = allsessions(s_idx);
    S    = V;
    S.session = sess;
    % -- Prepare data for this session:
    [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold, ...
     avg_ntrain_trials,avg_ntest_trials] = PrepareData_kfoldCV(S); 
    if V.preproc
    [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = ...
         PreProcess_kCV(Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold);
    end 
    
    % -- Computations for each CV partition for this session:
    for pt=1:V.kfold
        % Fit forward static linear model via LS regression:
        [M_star,xi] = FitM_UC(Ztrain_kfold{pt}',Xtrain_kfold{pt}');
        % Assess decoding performance: 
        MAE_cache(:,pt,s_idx) = getDecodePerf(Ztest_kfold{pt}, ...
                                              Xtest_kfold{pt}, M_star);
    end %over CV folds 
    M_cache{s_idx}      = M_star; 
    xi_cache(s_idx)     = xi; 
    Ztrain_cache{s_idx} = Ztrain_kfold{pt}';
    Xtrain_cache{s_idx} = Xtrain_kfold{pt}';
    Ztest_cache{s_idx}  = Ztest_kfold{pt}'; 
    Xtest_cache{s_idx}  = Xtest_kfold{pt}'; 
end %over all sessions 

O.M_cache      = M_cache; 
O.xi_cache     = xi_cache; 
O.Ztrain_cache = Ztrain_cache;
O.Xtrain_cache = Xtrain_cache;
O.Ztest_cache  = Ztest_cache;
O.Xtest_cache  = Xtest_cache; 

% Now compute average:
avg_bestdecodeperf_CV = mean(MAE_cache,3);  %average of best perf over sessions

% Cache results:
O.avg_bestdecodeperf_CV = avg_bestdecodeperf_CV;
O.avg_ntrain_trials     = avg_ntrain_trials;
O.avg_ntest_trials      = avg_ntest_trials; 

%% -- Save: -- 
if V.save
    savename = ['m0' num2str(V.mouse) '_bs' num2str(V.BaseSessions) '_AvgBestPerfkCV_behavs364.mat'];
    save([pwd '/Results/' savename],'O');
end

end

%% ** Internal Functions: ** 
% (1)______________________________________________________________________
function [allsessions,goodneurons] = getsessions(S)
    global LoadDir;

    % Compute test days: 
    TestSessions = NaN(length(S.BaseSessions),length(S.Diffs)); 
    for d_idx = 1:length(S.BaseSessions)
        TestSessions(d_idx,:) = S.BaseSessions(d_idx)+S.Diffs; 
    end 
    allsessions = unique([S.BaseSessions TestSessions(:)']);
    bd_vec      = allsessions; %base session values 

    % Find consistent set of good neurons over days:
    goodneurons_pre = []; 
    for b_idx = 1:length(bd_vec)
        bd    = bd_vec(b_idx); 
        if bd<10
           Datafile = ['m0' num2str(S.mouse) '_s0' num2str(bd) '.mat']; 
        else
           Datafile = ['m0' num2str(S.mouse) '_s' num2str(bd) '.mat']; 
        end 
        session_obj = load([LoadDir Datafile], 'session_obj'); 
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
end 
%__________________________________________________________________________
% (2)______________________________________________________________________
function [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold, ...
          avg_ntrain_trials,avg_ntest_trials] = PrepareData_kfoldCV(S)
    global LoadDir;

    % -- Specify k: --
    k            = S.kfold;
    Ztrain_kfold = cell(1,k); 
    Ztest_kfold  = cell(1,k);
    Xtrain_kfold = cell(1,k);
    Xtest_kfold  = cell(1,k); 
    
    if S.session<10
       Datafile = ['m0' num2str(S.mouse) '_s0' ...
                   num2str(S.session) '.mat']; 
    else
       Datafile = ['m0' num2str(S.mouse) '_s' ...
                   num2str(S.session) '.mat']; 
    end 

    session_obj = load([LoadDir Datafile], 'session_obj'); 
    session_obj = session_obj.session_obj; 
    if isnan(S.goodneurons(1))
       goodneurons = find(~isnan(session_obj.confidenceLabel));
    else
       goodneurons = S.goodneurons; 
    end 
    
    % -- Find Time Bins Comprising Correct Trials of Specified Type: -- 
    % Note that this removes inter-trial intervals (ITIs)
    O_tbininfo = ComputeActivity_SpecifyTrials('mouse',S.mouse,...
                                               'session',S.session); 

    if ischar(S.trial_type)
       trials_to_include = O_tbininfo.trialtimes_correct_all;
    elseif S.trial_type==2
       trials_to_include = O_tbininfo.trialtimes_correct_BR;
    elseif S.trial_type==3
       trials_to_include = O_tbininfo.trialtimes_correct_WL;
    end 

    ntotal_trials = length(trials_to_include); 
    
    ntest_trials_cache  = NaN(1,k);
    ntrain_trials_cache = NaN(1,k); 
    
    % Determine partitions for k-fold CV:
    for pt=1:k %define each test & train set
        test_idx1 = floor(((pt-1)/k)*ntotal_trials + 1); %trial indices 
        test_idx2 = floor((pt/k)*ntotal_trials);
        test_idxs = test_idx1:test_idx2; 
        train_idxs = setdiff(1:ntotal_trials,test_idxs); 
        ntest_trials_cache(pt)  = length(test_idxs);
        ntrain_trials_cache(pt) = length(train_idxs); 
         
        bins_to_include_test  = [];
        bins_to_include_train = [];
        for bi=test_idxs
            addon = trials_to_include(bi,1):trials_to_include(bi,2);
            bins_to_include_test = [bins_to_include_test; addon(:)]; %#ok<AGROW>
        end 
        for bi=train_idxs
            addon = trials_to_include(bi,1):trials_to_include(bi,2);
            bins_to_include_train = [bins_to_include_train; addon(:)]; %#ok<AGROW>
        end 
        
        % Prepare training set (all data) of Ca2+ dF/F values:
        Z = session_obj.timeSeries.calcium.data{1,1};
        Ztrain_kfold{pt} = Z(goodneurons,bins_to_include_train)';
        Ztest_kfold{pt}  = Z(goodneurons,bins_to_include_test)'; 
        
        % Prepare training set of kinematic values: 
        X = session_obj.timeSeries.virmen.data(S.behavvars,:);
        Xtrain_kfold{pt} = X(:,bins_to_include_train)';
        if sum(isnan(Xtrain_kfold{pt}))>0, error('NaN entries/n'); end
        Xtest_kfold{pt}  = X(:,bins_to_include_test)';
    end %over CV partitions 
    % Compute average # of trials for train & test sets:
    avg_ntrain_trials = mean(ntrain_trials_cache);
    avg_ntest_trials  = mean(ntest_trials_cache); 
end 
%__________________________________________________________________________
% (3)______________________________________________________________________
function [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = ...
         PreProcess_kCV(Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold)
    
    % -- Get # of CV partitions:
    k = length(Ztrain_kfold); 
    
    for pt = 1:k
        % Pre-processing: 
        %0-center:
        Z_train_mean     = mean(Ztrain_kfold{pt},1); 
        Ztrain_kfold{pt} = (Ztrain_kfold{pt} - Z_train_mean);
        Ztest_kfold{pt}  = (Ztest_kfold{pt} - Z_train_mean); 

        X_train_mean     = mean(Xtrain_kfold{pt},1);
        Xtrain_kfold{pt} = (Xtrain_kfold{pt} - X_train_mean);
        Xtest_kfold{pt}  = (Xtest_kfold{pt} - X_train_mean); 
    end 
end 
%__________________________________________________________________________
% (4)______________________________________________________________________
function [M_star,xi] = FitM_UC(Z,X)
        M_star = X*Z'*inv(Z*Z'); %#ok
        xi     = (norm(X - (M_star*Z),'fro'))^2; %RSS 
end
%__________________________________________________________________________
% (5)______________________________________________________________________
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

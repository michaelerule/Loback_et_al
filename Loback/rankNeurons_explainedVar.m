function [rankedNeurons_ev,best_ev,best_mae] = rankNeurons_explainedVar(varargin)
%__________________________________________________________________________
% Description: Rank orders the PPC neurons i \in [N] based on 
%              the fraction of variance explained using this subset
%              of neurons for prediction of the kinematic variable
%              values (i.e., R^2 for OLS). 
%              This uses Michael's great idea of storing the covariance 
%              & cross-covariance to speed up greedy search of the
%              "best" K-subset for each value K (1 <= K <= N). 
% Note:        This is a prerequisite analysis for main Fig. 2e
%
% Written 21 March 2019 by AL
% Updated 22 March 2019 by AL : Now also returns the mean absolute error
%                               (MAE) for each best subset
%__________________________________________________________________________

%% -- Initializations: --
V.mouse      = 4;          % mouse ID (can specify multiple IDs)
V.Nthresh    = 200;        % choose sessions that have >= this 
                           % threshold of high-confidence neurons           
V.behavvar   = 3;          % 3=forward position (units: m); 6=forward speed; 4=view angle
V.trial_type = 'all';      % 2=black right-turn trial;
                           % 3=white left-turn trial;
                           % 'all' specifies both trial types
V.nxval      = 10;         % # of folds (k) for cross-validation 
V            = parseargs(V, varargin{:});

%% -- Main Computations: -- 
S = V;
% Find all sessions with >= V.Nthresh high-confidence neurons:
% sess_list = FindSess_Nthresh(S);

S.session = 1; %sess_list(1); %Just choose one session for now

% Initial data preparation:
[Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold, ...
 avg_ntrain_trials,avg_ntest_trials] = PrepareData_kfoldCV(S); %#ok

% Initialize covariance structure:
Caa = cell(1,V.nxval);
Cab = cell(1,V.nxval); 

% Pre-process data (normalize): 
for k = 1:V.nxval 
       [Ztrain_kfold{k}, Ztest_kfold{k}, Xtrain_kfold{k}, Xtest_kfold{k}] = ...
        PreProcess(Ztrain_kfold{k}, Ztest_kfold{k}, Xtrain_kfold{k}, ...
                   Xtest_kfold{k}); 
       
       % Store covariances, cross-covariances
       Caa{k} = Ztrain_kfold{k}' * Ztrain_kfold{k}; 
       Cab{k} = Ztrain_kfold{k}' * Xtrain_kfold{k}; 
end 

% Track the set of units we've added / not added
unused = 1:size(Ztrain_kfold{1},2); % full index set of all neurons
used = []; 
best_ev  = []; % init cache for best R^2
best_mae = []; % init cache for best MAE

while ~isempty(unused)
    % Consider adding every possible unused neuron to current set
    ev = [];
    for unit = unused
        ev = [ev explained_variance([used unit], Caa, Cab, Ztest_kfold,...
                                    Xtest_kfold, V.nxval)]; %#ok 
    end 
    % Find the neuron that adds the most information
    ib = find(ev == max(ev), 1, 'first'); 
    % Add this best neuron to the set and iterate
    u = unused(ib);
    best_ev  = [best_ev ev(ib)]; %#ok
    % fprintf('%02d ', u);
    % fprintf('%02f\n', ev(ib)); 
    unused   = setdiff(unused,u); 
    used     = [used u]; %#ok
    best_mae = [best_mae compute_mae(used, Caa, Cab, Ztest_kfold,...
                                     Xtest_kfold, V.nxval)]; %#ok
end 

rankedNeurons_ev = used;

end

%% -- Internal Functions: -- 
% (0)______________________________________________________________________
function [sess_list] = FindSess_Nthresh(S)
    %Init: 
    LoadDir   = ['../Data/m0' num2str(S.mouse) '/'];
    sess_list = []; 
    
    for sess = 1:20
        if sess<10
            Datafile = ['m0' num2str(S.mouse) '_s0' ...
            num2str(sess) '.mat']; 
        else
            Datafile = ['m0' num2str(S.mouse) '_s' ...
            num2str(sess) '.mat']; 
        end 
        session_obj = load([LoadDir Datafile], 'session_obj'); 
        session_obj = session_obj.session_obj; 
        goodneurons = find(session_obj.confidenceLabel==1 | ...
                           session_obj.confidenceLabel==2); 
        if (length(goodneurons)>= S.Nthresh) && ...
           session_obj.numConditions(sess)==2
            sess_list = [sess_list sess]; %#ok<AGROW>
        end 
    end 
end 
%__________________________________________________________________________
% (1)______________________________________________________________________
function [Ztrain_kfold, Ztest_kfold, Xtrain_kfold, Xtest_kfold, ...
          avg_ntrain_trials, avg_ntest_trials] = PrepareData_kfoldCV(S)
      
    % -- Specify k: --
    k            = S.nxval;
    Ztrain_kfold = cell(1,k); 
    Ztest_kfold  = cell(1,k);
    Xtrain_kfold = cell(1,k);
    Xtest_kfold  = cell(1,k); 
    
    LoadDir = ['../Data/m0' num2str(S.mouse) '/'];
    if S.session<10
       Datafile = ['m0' num2str(S.mouse) '_s0' ...
                   num2str(S.session) '.mat']; 
    else
       Datafile = ['m0' num2str(S.mouse) '_s' ...
                   num2str(S.session) '.mat']; 
    end 

    session_obj = load([LoadDir Datafile], 'session_obj'); 
    session_obj = session_obj.session_obj; 
    goodneurons = find(~isnan(session_obj.confidenceLabel));
    
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
    for pt=1:k % define each test & train set
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
        X = session_obj.timeSeries.virmen.data(S.behavvar,:);
        Xtrain_kfold{pt} = X(:,bins_to_include_train)';
        if sum(isnan(Xtrain_kfold{pt}))>0, error('NaN entries/n'); end
        Xtest_kfold{pt}  = X(:,bins_to_include_test)';
    end % over CV partitions 
    % Compute average # of trials for train & test sets:
    avg_ntrain_trials = mean(ntrain_trials_cache);
    avg_ntest_trials  = mean(ntest_trials_cache); 
end 
%__________________________________________________________________________
% (2)______________________________________________________________________
function [Z_train, Z_test, X_train, X_test,...
          Z_train_mean, X_train_mean] = PreProcess(Z_train, Z_test,...
                                                  X_train, X_test)
    % Implement pre-processing 
    % zero-center
    Z_train_mean = mean(Z_train,1); 
    Z_train      = (Z_train - Z_train_mean); 
    Z_test       = (Z_test - Z_train_mean);  

    % zero-center
    X_train_mean = mean(X_train, 1);
    X_train      = (X_train - X_train_mean);
    X_test       = (X_test - X_train_mean); 
end 
%__________________________________________________________________________
% (3)______________________________________________________________________
function [avg_ev] = explained_variance(us, Caa, Cab, testA, testB, nxval)

    ev_vec = NaN(1, nxval); 

    for k = 1:nxval
        % Compute the OLD weights from the training data
        % Use the cached covariances for speed!
        w_star = Caa{k}(us,us)\Cab{k}(us);
        
        % Compute the decoded estimate
        yhat = testA{k}(:,us) * w_star; 
        
        % Compute R^2 (variance explained)
        SSres = compute_RSS(testB{k}, yhat); 
        SStot = sum(testB{k}.^2); 
        ev_vec(k) = 1 - (SSres/SStot); % R^2
    end
    
    avg_ev = mean(ev_vec); % mean explained variance, averaged over CV folds
end 
%__________________________________________________________________________
% (4)______________________________________________________________________
function [SSres] = compute_RSS(testy, yhat)
    res = testy(:) - yhat(:);
    SSres = sum(res.^2); 
end
%__________________________________________________________________________
% (5)______________________________________________________________________
function [avg_MAE] = compute_mae(us, Caa, Cab, testA, testB, nxval)

    mae_vec = NaN(1, nxval);
    for k = 1:nxval
        % Compute the OLD weights from the training data
        w_star = Caa{k}(us,us)\Cab{k}(us);
        
        % Compute the decoded estimate
        yhat = testA{k}(:,us) * w_star; 
        
        % Compute MAE
        mae_vec(k) = mean(abs(testB{k} - yhat)); 
    end 
    
    avg_MAE = mean(mae_vec); 
end 
%__________________________________________________________________________
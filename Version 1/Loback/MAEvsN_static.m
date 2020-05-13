function [] = MAEvsN_static(varargin)
%__________________________________________________________________________
% Description: Prepares data for extrapolation to assess performance  
%              of the *static decoder* to address the question: 
%              How precisely could we decode with a larger population?
%              (For main Fig. 2) 
%
% Written 31 January 2019 by AL
%__________________________________________________________________________

%% -- Initializations: -- 
V.mice         = 4;         %mouse ID (can specify multiple IDs) 
V.Nthresh      = 150;       %choose sessions that have >= this threshold 
                            %of high-confidence neurons           
V.behavvars    = [3 6 4];   %3=forward position (units: m); 6=forward speed
V.trial_type   = 'all';     %2 = black right-turn trial;
                            %3 = white left-turn trial;
                            %'all' specifies both trial types
V.training_range = [0,0.8];
V.valid_range    = [0.8,1];   
V.preproc        = 1;       %1=preprocess data; 0=don't preprocess (for affine model)  
V.niters         = 50;      %# of iterations (since use permutations) 
V                = parseargs(V, varargin{:});  

% Specify N values:
step_size = 5;                 
Nvec      = 50:step_size:V.Nthresh;  
mcell     = cell(1,length(V.mice));

%% -- Main Computations: -- 
for m_idx = 1:length(V.mice) 
    mouse = V.mice(m_idx); 
    S     = V;
    S.mouse = mouse; 
    
    %Find all sessions with >= V.Nthresh high-confidence neurons: 
    sess_list = FindSess_Nthresh(S);
    
    %Init cache:
    MAE_cell = cell(length(sess_list),S.niters);
    
    for s_idx = 1:length(sess_list) 
        S.session = sess_list(s_idx);
        %Initial data preparation:
        [Z_train,Z_test,X_train,X_test] = PrepareData(S);
        
        for iter = 1:S.niters
            %Init cache: 
            MAE_kin_cache = NaN(length(S.behavvars),length(Nvec));

            for N_ix = 1:length(Nvec)
                O = StatLinModel(Z_train,Z_test,X_train,X_test, ...
                                 'mouse',S.mouse, ...
                                 'preproc',S.preproc, ...
                                 'N', Nvec(N_ix), 'viz',0, 'save',0); 
                MAE_kin_cache(:,N_ix) = O.mean_abserror_kin; 
            end 

            MAE_cell{s_idx,iter} = MAE_kin_cache; %cache 
        end %over permutation iterations
    end %over sessions
    
    %Cache: 
    mcell{m_idx} = MAE_cell(:)'; 
    
    % -- Process Results for Exporting to Python: -- 
    R_cache = NaN(p,length(Nvec)); %init 
    for idx = 1:p
        mean_kmvar = cell2mat(cellfun(@(x) x(idx,:), mcell{m_idx}, ...
                              'UniformOutput',0)'); 
        R_cache(idx,:) = mean(mean_kmvar); 
    end 
end 

%% -- Save: -- 
savename = ['MAEvsN_static_m0' num2str(V.mouse) '_Nthresh' ...
            num2str(V.Nthresh) '_' num2str(step_size) '.csv'];
save([pwd '/Results/' savename], 'Output');

% Put into .csv format to be read into Python: 
csvwrite(savename,[Nvec(:) R_cache(:)]);

end %main 

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
end %FindSess_Nthresh 
%__________________________________________________________________________
% (1)______________________________________________________________________
function [Z_train,Z_test,X_train,X_test] = PrepareData(S)
    LoadDir     = ['../Data/m0' num2str(S.mouse) '/'];
    train_ratio = S.training_range(end); 
    if S.session<10
        Datafile = ['m0' num2str(S.mouse) '_s0' ...
                    num2str(S.session) '.mat']; 
    else 
        Datafile = ['m0' num2str(S.mouse) '_s' ...
                    num2str(S.session) '.mat']; 
    end 

    session_obj = load([LoadDir Datafile], 'session_obj'); 
    session_obj = session_obj.session_obj; 
    goodneurons = find(session_obj.confidenceLabel==1 | ...
                       session_obj.confidenceLabel==2); 

    % -- Find Time Bins Comprising Correct Trials of Specified Type: -- 
    O_tbininfo = ComputeActivity_SpecifyTrials('mouse',S.mouse,...
                 'session',S.session); 
                 
    if ischar(S.trial_type)
       trials_to_include = O_tbininfo.trialtimes_correct_all;
    elseif S.trial_type==2
       trials_to_include = O_tbininfo.trialtimes_correct_BR;
    elseif S.trial_type==3
       trials_to_include = O_tbininfo.trialtimes_correct_WL; 
    end 

    % Partition total trials into training & test set:
    total_ntrials = size(trials_to_include,1); 
    goodtrials = trials_to_include(randperm(size(trials_to_include,1)),:);
 
    % ** Compute number of train & test trials: **
    ntrials_train = floor(total_ntrials*train_ratio);
    ntrials_test  = total_ntrials-ntrials_train; 

    % Compute time bin indices to include:
    bins_to_include_train = [];
    bins_to_include_test  = [];
    for bi=1:ntrials_train
        addon = goodtrials(bi,1):goodtrials(bi,2); 
        bins_to_include_train = [bins_to_include_train; addon(:)]; %#ok
    end 
    for bi=1:ntrials_test
        addon = goodtrials(bi+ntrials_train,1):goodtrials(bi+ntrials_train,2);
        bins_to_include_test  = [bins_to_include_test; addon(:)];  %#ok
    end 

    % Prepare training set (all data) of Ca2+ dF/F values:
    Z = session_obj.timeSeries.calcium.data{1,1}; %{2,1} = use deconvolved
    Z_train = Z(goodneurons,bins_to_include_train)'; 
    Z_test  = Z(goodneurons,bins_to_include_test)';   

    % Prepare train set of kinematic values: 
    Y = session_obj.timeSeries.virmen.data(S.behavvars,:);
    if sum(isnan(Y))>0, error('NaN entries/n'); end 
    X_train = Y(:,bins_to_include_train)';          
    X_test  = Y(:,bins_to_include_test)';             
end
%__________________________________________________________________________

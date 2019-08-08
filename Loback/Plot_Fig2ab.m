function [] = Plot_Fig2abc(varargin)
%__________________________________________________________________________
% Description: Generates main Fig. 2 for paper (first data analysis
%              result - i.e., static linear decoding results).
%
% Updated 4 Feb. 2019 by AL : added 10-fold CV to assess performance
%__________________________________________________________________________

%% -- Initializations: -- 
V.mice         = [1 3 4 5];             
V.behavvars    = [3 6 4];    %3=forward position (units: m); 6=forward speed; 4=view angle
V.trial_type   = 'all';      %2 = black right-turn trial;
                             %3 = white left-turn trial;
                             %'all' specifies both trial types
V.preproc      = 1;          %1=preprocess data; 0=don't                  
V.kfold        = 10;         %# of folds (k) for cross-validation 
V.deconv       = 0;          %1=use deconvolved; 0=use raw deltaF/F
V.Nthresh      = 200;        %only include sessions w/ >= Nthresh neurons
V.panelB       = 0;
V.panelC       = 1; 
V.kfold        = 10;         %# of folds (k) for cross-validation 
V              = parseargs(V, varargin{:});

p          = length(V.behavvars);
Labels1    = {'MAE Position (m)', 'MAE Speed (m/s)', ...
              'MAE View Angle ($^{\circ}$)'}; 
Labels2    = {'Forward Position (m)', 'Forward Speed (m/s)', ...
              'View Angle ($^{\circ}$)'};

%savedir    = [pwd '/PreparedData/'];
savedir    = ['/home/mrule/Desktop/'];
savename   = 'Fig1_panelA_results.mat';

%% -- Define colors for plots: -- 
tableau20 = {[ 31, 119, 180], [174, 199, 232], [255, 127,  14], [255, 187, 120], ...    
             [ 44, 160,  44], [152, 223, 138], [214,  39,  40], [255, 152, 150], ... 
             [148, 103, 189], [197, 176, 213], [140,  86,  75], [196, 156, 148], ...    
             [227, 119, 194], [247, 182, 210], [127, 127, 127], [199, 199, 199], ...   
             [188, 189,  34], [219, 219, 141], [ 23, 190, 207], [158, 218, 229]}; 
% Scale the RGB values to the [0, 1] range, 
% which is the format matplotlib accepts  
for i=1:length(tableau20)  
    tableau20{i} = tableau20{i}./255.0; 
end 

kin_colors = {tableau20{1} tableau20{3} tableau20{5}}; 

LoadDir = '/home/mrule/Workspace2/PPC_data/';
fprintf(1,'Loading datasets from %s,\n',LoadDir);
fprintf(1,'please re-define LoadDir to match the path to the data files on your system.\n');
fprintf(1,'Datasets are available from the Harvey lab.\n');

%% -- Panel c (Summary Fig): -- 
if ~(exist([savedir savename],'file')==2) %load if exists
    MAE_cell  = cell(length(V.mice),1);  
    sess_cell = cell(length(V.mice),1); 

    for m_idx = 1:length(V.mice)
        mouse = V.mice(m_idx); 
        S     = V;
        S.mouse = mouse;

        %Find all sessions with >= V.Nthresh high-confidence neurons: 
        sess_list = FindSess_Nthresh(S);
        sess_cell{m_idx} = sess_list;

        %Init caches:
        MAE_midx = NaN(p,V.kfold,length(sess_list)); %init

        %Fit affine model for each session in sess_list: 
        for s_idx = 1:length(sess_list)
            S.session = sess_list(s_idx);
            % -- Prepare data for this session:
            [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold, ...
             avg_ntrain_trials,avg_ntest_trials] = PrepareData_kfoldCV(S); %#ok
            
            if V.preproc
            [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = ...
                 PreProcess_kCV(Ztrain_kfold,Ztest_kfold,...
                 Xtrain_kfold,Xtest_kfold); 
            end 
            
            % -- Computations for each CV partition for this session:
            for pt=1:V.kfold
                % Fit forward static linear model via LS regression:
                [M_star,~] = FitM_UC(Ztrain_kfold{pt}',Xtrain_kfold{pt}');
                % Assess decoding performance: 
                MAE_midx(:,pt,s_idx) = getDecodePerf(Ztest_kfold{pt}, ...
                                                      Xtest_kfold{pt}, ...
                                                      M_star);
            end %over CV folds
    
        end %over sessions in sess_list

        MAE_cell{m_idx} = reshape(mean(MAE_midx,2),...
                                  [p,length(sess_list)]); %cache
    end %over mice

    % -- Save: --
    save([savedir savename], 'MAE_cell', 'sess_cell'); 
else
    %Load MAE_cell & sess_cell 
    load([savedir savename]); %#ok
end %if 

% -- Visualize plot for panel a: -- 
figure(2); clf;
%Panel c: 
sp_idxs = [3 6 9]; 
ylims_c = {[0.2 0.6], [0.02 0.15], [0 25]}; 
avg_MAE_cache = NaN(p,length(V.mice)); 
std_MAE_cache = NaN(p,length(V.mice)); 
for kin_idx = 1:p
    for m_idx = 1:length(V.mice)
        avg_MAE_cache(kin_idx,m_idx) = mean(MAE_cell{m_idx}(kin_idx,:)); 
        std_MAE_cache(kin_idx,m_idx) = std(MAE_cell{m_idx}(kin_idx,:)); 
    end 
    subplot(3,3,sp_idxs(kin_idx)); 
    e = errorbar(1:length(V.mice),avg_MAE_cache(kin_idx,:),...
                 std_MAE_cache(kin_idx,:),'Color',[kin_colors{kin_idx}], ...
                 'LineStyle','none','Marker','o','LineWidth',3); 
    e.Color(4) = 0.2;
    axis square; box off; 
    ax = gca; 
    ax.XTick = 1:1:length(V.mice);
    ax.XTickLabel = {'m1' 'm3' 'm4' 'm5'}; 
    ylabel(Labels1{kin_idx},'Interpreter','Latex');
    xlim([0.5 length(V.mice)+0.5]); 
    ylim(ylims_c{kin_idx}); 
    ax.YTick      = ylims_c{kin_idx}; 
end 

% -- Format font etc. -- 
a=findobj(gcf); alltext=findall(a,'Type','text');
allaxes=findall(a,'Type','axes');
set(allaxes,'FontName','Helvetica','FontWeight','Normal','FontSize',18);
set(alltext,'FontName','Helvetica','FontWeight','Normal','FontSize',20);

%% -- Panels a & b (examples): -- 
%Load MAE_cell & sess_cell 
load([savedir savename]); %#ok 
m_indices = [3 4]; 
sp_idxs   = [1 2; 4 5; 7 8]; 
ylims_a   = {[-2 6], [-0.5 1], [-150 100]}; 

for i = 1:length(m_indices)
    m_idx = m_indices(i); 
    for kin_idx = 1:p
        best_sidx = find(MAE_cell{m_idx}(kin_idx,:) == ...
                    min(MAE_cell{m_idx}(kin_idx,:))); 
        best_sess = sess_cell{m_idx}(best_sidx);      %#ok 

        % -- Prepare data: --
        S         = V; 
        S.mouse   = V.mice(m_idx);
        S.session = best_sess; 
        S.training_range = [0,1-(1/V.kfold)];
        S.valid_range    = [(1/V.kfold),1];  
        [Z_train,Z_test,X_train,X_test] = PrepareData(S);
        ones_add_train = ones(1,size(Z_train,1)); 
        Z_tilde_train  = [Z_train ones_add_train(:)];
        %Fit affine model: 
        O = fit_affine_model(Z_tilde_train', Z_test', ...
                             X_train', X_test');
        subplot(3,3,sp_idxs(kin_idx,i));
        plot1 = plot(X_test(:,kin_idx),'k','linewidth',3); 
        plot1.Color(4) = 0.3;
        hold on;
        plot2 = plot(O.xk_preds(kin_idx,:),'color',...
                     kin_colors{kin_idx},'linewidth',3);
        plot2.Color(4) = 0.5;
        axis square; xlim([0 300]); 
        ax=gca; ax.XTick = [0 300];
        ylim(ylims_a{kin_idx}); 
        ax.YTick = ylims_a{kin_idx}; 
        if i==1
            ylabel(Labels2{kin_idx},'Interpreter','Latex');
        end 
        if kin_idx==p
            xlabel('Time Bin','Interpreter','Latex'); 
        end 
        box off; 
    end %over kinematic vars
end %over mice

% -- Format font etc. -- 
a=findobj(gcf); alltext=findall(a,'Type','text');
allaxes=findall(a,'Type','axes');
set(allaxes,'FontName','Helvetica','FontWeight','Normal','FontSize',18);
set(alltext,'FontName','Helvetica','FontWeight','Normal','FontSize',20);

end %main

%% -- Internal Functions: --
% (0)______________________________________________________________________
function [sess_list] = FindSess_Nthresh(S)
    global LoadDir
    %Init: 
    sess_list = []; 
    
    for sess = 1:200
        Datafile = sprintf('m%02d_s%02d.mat',S.mouse,sess);
        try
            session_obj = load([LoadDir Datafile], 'session_obj'); 
            session_obj = session_obj.session_obj; 
            goodneurons = find(session_obj.confidenceLabel==1 | ...
                               session_obj.confidenceLabel==2); 
            if (length(goodneurons)>= S.Nthresh) && ...
               session_obj.numConditions(sess)==2
                sess_list = [sess_list sess]; %#ok<AGROW>
            end 
        catch
             fprintf(1,'Error opening %s, ignoring\n',[LoadDir Datafile]);
        end
    end 
end

%FindSess_Nthresh
%__________________________________________________________________________
% (1)______________________________________________________________________
function [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold, ...
          avg_ntrain_trials,avg_ntest_trials] = PrepareData_kfoldCV(S)
    global LoadDir
    
    % -- Specify k: --
    k            = S.kfold;
    Ztrain_kfold = cell(1,k); 
    Ztest_kfold  = cell(1,k);
    Xtrain_kfold = cell(1,k);
    Xtest_kfold  = cell(1,k); 
    
    Datafile = sprintf('m%02d_s%02d.mat',S.mouse,S.session);
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
% (2)______________________________________________________________________
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
function [Z_train,Z_test,X_train,X_test] = PrepareData(S)
    global LoadDir
    Datafile = sprintf('m%02d_s%02d.mat',S.mouse,S.session);
    session_obj = load([LoadDir Datafile], 'session_obj'); 
    session_obj = session_obj.session_obj; 
    goodneurons = find(session_obj.confidenceLabel==1 | ...
                   session_obj.confidenceLabel==2); 
    
    % -- Find Time Bins Comprising Correct Trials of Specified Type: -- 
    % Note that this removes inter-trial intervals (ITIs)
    O_tbininfo = ComputeActivity_SpecifyTrials('mouse',S.mouse,'session',S.session); 

    if ischar(S.trial_type)
       trials_to_include = O_tbininfo.trialtimes_correct_all;
    elseif S.trial_type==2
       trials_to_include = O_tbininfo.trialtimes_correct_BR;
    elseif S.trial_type==3
       trials_to_include = O_tbininfo.trialtimes_correct_WL;
    end 

    ntotal_trials = length(trials_to_include); 
    ntrials_train = floor(ntotal_trials*S.training_range(2));
    ntrials_test  = ntotal_trials-ntrials_train; 

    % Compute time bin indices to include:
    bins_to_include_train = [];
    bins_to_include_test  = []; 
    for bi=1:ntrials_train
        addon = trials_to_include(bi,1):trials_to_include(bi,2); 
        bins_to_include_train = [bins_to_include_train; addon(:)]; %#ok<AGROW>
    end 
    for bi=1:ntrials_test
        addon = trials_to_include(bi+ntrials_train,1):...
                trials_to_include(bi+ntrials_train,2);
        bins_to_include_test  = [bins_to_include_test; addon(:)];  %#ok
    end 

    % Prepare training set (all data) of Ca2+ dF/F values:
    if S.deconv
        Z   = session_obj.timeSeries.calcium.data{2,1};
    else
        Z   = session_obj.timeSeries.calcium.data{1,1};
    end 
    Z_train = Z(goodneurons,bins_to_include_train);
    Z_train = Z_train'; 
    Z_test  = Z(goodneurons,bins_to_include_test); 
    Z_test  = Z_test'; 

    % Prepare training set of kinematic values: 
    X = session_obj.timeSeries.virmen.data(S.behavvars,:);
    X_train = X(:,bins_to_include_train);
    if sum(isnan(X_train))>0, error('NaN entries/n'); end 
    X_train = X_train';
    X_test  = X(:,bins_to_include_test); 
    X_test  = X_test';
end
%__________________________________________________________________________
% (6)______________________________________________________________________
function [O] = fit_affine_model(Z_train,Z_test,X_train,X_test)
    [ntest_samples,p] = size(X_test');
    % -- Obtain analytical MSE parameters: -- 
    theta_star = X_train*Z_train'*inv(Z_train*Z_train'); %#ok
    W_star     = theta_star(:,1:end-1); 
    beta       = theta_star(:,end); 
    N          = size(Z_test,1);
    
    % Initialize Caches: 
    xk_preds   = zeros(p,ntest_samples); 
    evec_kin   = zeros(p,ntest_samples); 
    
    for k=1:ntest_samples
        zk            = Z_test(:,k);
        xk            = X_test(:,k);         %actual kinematic value
        xk_hat        = W_star*zk(:) + beta; %model prediction   
        evec_kin(:,k)   = xk(:) - xk_hat(:);
        xk_preds(:,k)   = xk_hat;
    end 
    
    MAE        = mean(abs(evec_kin),2);      %MAE on held-out test data
    MAE(1:2)   = MAE(1:2).*100;              %units now: cm & cm/s
    O.MAE      = MAE; 
    O.N        = N; 
    O.xk_preds = xk_preds;
    
end %fit_affine_model
%__________________________________________________________________________

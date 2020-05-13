function [Output] = DecayRate_StatLin_kfoldCV(varargin)
%__________________________________________________________________________
% Description: For part (1ii,a) of the paper analyses, i.e. aims to
%              assess rate of decay of the performance of the best 
%              linear decoder on session s for sessions s+1,s+2,...,s+D.
% Written 17 Dec 2018 by AL
%__________________________________________________________________________

%% -- Initializations: -- 
V.mouse          = 3;     
V.BaseSession    = 7;             %7      
V.Diffs          = [0 1 2 3];     %[-5 -1 0 1 4];             
V.behavvars      = [2 3 5 6 4];   %3=forward position (units: m); 6=forward speed
V.trial_type     = 'all';         %2=black right-turn trial;
                                  %3=white left-turn trial;
                                  %'all' specifies both trial types
V.kfold          = 10;            %# of folds (k) for cross-validation  
V.preproc        = 1;             %1=preprocess data; 0=don't preprocess
V.viz            = 1;             %1=visualize; 0=don't
V.save           = 1;             %1=save; 0=don't save 
V.N              = 'all';         %'all'=use all; else use # specified
V                = parseargs(V, varargin{:});  
                
Output.BaseSession = V.BaseSession;
Output.Diffs       = V.Diffs;

[allsessions,goodneurons,days] = getsessions(V); 
V.goodneurons                  = goodneurons; 
p                              = length(V.behavvars); 
MAE_abs                        = cell(1,length(allsessions)); %absolute MAE
MAE_percentchange              = MAE_abs; %units: percent change                                         

%% -- Fit Static Linear Model to Base Session: -- 
% Prepare data: 
S = V;
S.session = V.BaseSession; 
[Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = ...
 PrepareData_kfoldCV(S);
if V.preproc
   [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = ...
    PreProcess_kCV(Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold);
end

M_bd_kCV        = cell(1,V.kfold); 
MAE_baseday_kCV = NaN(p,V.kfold); 
for pt=1:V.kfold
    M_bd_kCV{pt} = FitM_UC(Ztrain_kfold{pt}',Xtrain_kfold{pt}'); 
    MAE_baseday_kCV(:,pt) = getDecodePerf(Ztest_kfold{pt},...
                                        Xtest_kfold{pt},M_bd_kCV{pt});
end 

MAE_abs{allsessions==V.BaseSession} = MAE_baseday_kCV;
MAE_percentchange{allsessions==V.BaseSession} = zeros(p,V.kfold); 

%% -- Test Performance on Other Days: -- 
test_idxs = find(allsessions~=V.BaseSession); 
S = V;
for tdi = 1:length(test_idxs)
    S.session       = allsessions(test_idxs(tdi));
    MAE_testday_kCV = NaN(p,V.kfold); %initialize
    MAE_pchange_kCV = NaN(p,V.kfold); %initialize
    [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = ...
                                           PrepareData_kfoldCV(S);
    if V.preproc
       [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = ...
        PreProcess_kCV(Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold);
    end
    for pt=1:V.kfold
        MAE_testday_kCV(:,pt) = getDecodePerf(Ztest_kfold{pt}, ...
                                              Xtest_kfold{pt}, ...
                                              M_bd_kCV{pt});
        
        %How much worse is this MAE vs. if used best decoder for this day?
        M_star = FitM_UC(Ztrain_kfold{pt}',Xtrain_kfold{pt}');  
        MAE_sameday = getDecodePerf(Ztest_kfold{pt},Xtest_kfold{pt}, ...
                                    M_star);
        MAE_diff = MAE_testday_kCV(:,pt) - MAE_sameday;
        MAE_pchange_kCV(:,pt) = ( MAE_diff ./ MAE_sameday)*100; 
    end 
    MAE_abs{test_idxs(tdi)} = MAE_testday_kCV;
    MAE_percentchange{test_idxs(tdi)} = MAE_pchange_kCV;
end

%% -- Compute Mean & Std: --
mean_MAE_abs = cell2mat(cellfun(@(x) mean(x,2), MAE_abs, ...
                        'UniformOutput', 0));
std_MAE_abs = cell2mat(cellfun(@(x) std(x,1,2), MAE_abs, ...
                        'UniformOutput', 0)); 
                   
mean_MAE_pdiff = cell2mat(cellfun(@(x) mean(x,2), MAE_percentchange, ...
                          'UniformOutput', 0));  
std_MAE_pdiff = cell2mat(cellfun(@(x) std(x,1,2), MAE_percentchange, ...
                          'UniformOutput', 0));                       

%% -- Save Results: -- 
Output.MAE_abs           = MAE_abs;
Output.MAE_percentchange = MAE_percentchange; 
Output.mean_MAE_abs      = mean_MAE_abs;
Output.std_MAE_abs       = std_MAE_abs;
Output.mean_MAE_pdiff    = mean_MAE_pdiff;
Output.std_MAE_pdiff     = std_MAE_pdiff; 
Output.goodneurons       = goodneurons; 
Output.days              = days; 
N                        = size(Ztrain_kfold{1},2); 
Output.N                 = N; 

if V.save
   savename = ['m0' num2str(V.mouse) '_DecayMAE_StatLin_kfoldCV_bs' ...
               num2str(V.BaseSession) '.mat'];
   save([pwd '/Results/' savename],'Output');
end 

%% -- Visualize Results: -- 
if V.viz
    DeltaDays = days - days(3); 
    figure(10); hold on;
    subplot(322); %forward position (y)
    yyaxis left; 
    plot1 = shadedErrorBar(DeltaDays,mean_MAE_abs(2,:),std_MAE_abs(2,:),...
                   'lineprops',{'b-o','markerfacecolor','b',...
                   'MarkerSize',4,'LineWidth',2}); 
    plot1.Color(4) = 0.1; 
    plot1.MarkerFaceAlpha = 0.1;
    xlim([DeltaDays(1) DeltaDays(end)]); axis square;
    ylabel('MAE Position (m)','Interpreter','Latex'); 
    ax = gca; ax.XTick = DeltaDays;
    yyaxis right; 
    plot2 = shadedErrorBar(DeltaDays,mean_MAE_pdiff(2,:),std_MAE_pdiff(2,:),...
                   'lineprops',{'r-o','markerfacecolor','r',...
                   'MarkerSize',4,'LineWidth',2});
    plot2.Color(4) = 0.1; 
    plot2.MarkerFaceAlpha = 0.1;
    ylabel('$\%$ Change vs. Same Day','Interpreter','Latex');
    title(['Mouse 3 ($N=' num2str(N) '$)'],'Interpreter','Latex');
    box off;

    subplot(324); %forward velocity (y velocity)
    yyaxis left; 
    plot1 = shadedErrorBar(DeltaDays,mean_MAE_abs(4,:),std_MAE_abs(4,:),...
                   'lineprops',{'b-o','markerfacecolor','b',...
                   'MarkerSize',4,'LineWidth',2}); 
    plot1.Color(4) = 0.1; 
    plot1.MarkerFaceAlpha = 0.1;
    xlim([DeltaDays(1) DeltaDays(end)]); axis square;
    ylabel('MAE Velocity (m/s)','Interpreter','Latex'); 
    ax = gca; ax.XTick = DeltaDays;
    yyaxis right; 
    plot2 = shadedErrorBar(DeltaDays,mean_MAE_pdiff(4,:),std_MAE_pdiff(4,:),...
                   'lineprops',{'r-o','markerfacecolor','r',...
                   'MarkerSize',4,'LineWidth',2});
    plot2.Color(4) = 0.1; 
    plot2.MarkerFaceAlpha = 0.1;
    ylabel('$\%$ Change vs. Same Day','Interpreter','Latex');
    box off;

    subplot(326); %view angle
    yyaxis left; 
    plot1 = shadedErrorBar(DeltaDays,mean_MAE_abs(5,:),std_MAE_abs(5,:),...
                   'lineprops',{'b-o','markerfacecolor','b',...
                   'MarkerSize',4,'LineWidth',2}); 
    plot1.Color(4) = 0.1; 
    plot1.MarkerFaceAlpha = 0.1;
    xlim([DeltaDays(1) DeltaDays(end)]); axis square;
    ylabel('MAE View Angle ($^\circ$)','Interpreter','Latex'); 
    ax = gca; ax.XTick = DeltaDays;
    yyaxis right; 
    plot2 = shadedErrorBar(DeltaDays,mean_MAE_pdiff(5,:),std_MAE_pdiff(5,:),...
                   'lineprops',{'r-o','markerfacecolor','r',...
                   'MarkerSize',4,'LineWidth',2});
    plot2.Color(4) = 0.1; 
    plot2.MarkerFaceAlpha = 0.1;
    ylabel('$\%$ Change vs. Same Day','Interpreter','Latex');
    xlabel('$\Delta$ Days','Interpreter','Latex');
    box off;

    % -- Format font etc. -- 
    a=findobj(gcf); alltext=findall(a,'Type','text');
    allaxes=findall(a,'Type','axes');
    set(allaxes,'FontName','Helvetica','FontWeight','Normal','FontSize',15);
    set(alltext,'FontName','Helvetica','FontWeight','Normal','FontSize',16);
end 

end %main

% (1)______________________________________________________________________
function [allsessions,goodneurons,days] = getsessions(S)
    global LoadDir;
    % Compute test days: 
    TestSessions = NaN(length(S.BaseSession),length(S.Diffs)); 
    for d_idx = 1:length(S.BaseSession)
        TestSessions(d_idx,:) = S.BaseSession(d_idx)+S.Diffs; 
    end 
    allsessions = unique([S.BaseSession TestSessions(:)']);
    bd_vec      = allsessions; %base session values 

    % Find consistent set of good neurons over days:
    goodneurons_pre = []; 
    days = NaN(size(allsessions)); 
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
        % Get the day that corresponds to each session ID:
        days(b_idx) = session_obj.deltaDays(bd); 
    end 
    
    if ischar(S.N)
        goodneurons = goodneurons_pre;
    else
        goodneurons = goodneurons_pre(randperm(S.N)); 
    end 
end 
% (2)______________________________________________________________________
function [Ztrain_kfold,Ztest_kfold,Xtrain_kfold,Xtest_kfold] = PrepareData_kfoldCV(S)
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
    
    % Determine partitions for k-fold CV:
    for pt=1:k %define each test & train set
        test_idx1 = floor(((pt-1)/k)*ntotal_trials + 1); %trial indices 
        test_idx2 = floor((pt/k)*ntotal_trials);
        test_idxs = test_idx1:test_idx2; 
        train_idxs = setdiff(1:ntotal_trials,test_idxs); 
         
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
end 
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
% (4)______________________________________________________________________
function [mae_kin] = getDecodePerf(Z_test,X_test,M_b)
    [ntest_samples,p] = size(X_test);
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
% (5)______________________________________________________________________
function M_star = FitM_UC(Z,X)
        M_star = X*Z'*inv(Z*Z'); %#ok
end 
%__________________________________________________________________________

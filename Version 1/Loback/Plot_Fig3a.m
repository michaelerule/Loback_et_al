function [] = Plot_Fig3abcef(varargin)
%__________________________________________________________________________
% Description: Creates panels a-c and e & f of main Fig. 3, 
%              i.e. everything except for the cartoon schematic
%              (which is panel d). 
%              Note that panel (d) will be genearted separately. 
% Written 04/01/2019 by AL 
% Updated 15/01/2019 by AL - added optimisation analysis 
% Updated 25/03/2019 by AL - now reports absolute magnitude for 
%                            bottom panels of Fig. 3e,f 
% Updated 26/03/2019 by AL - rearranged organization of panels
%__________________________________________________________________________

%% -- Initializations: -- 
V.m1          = 4;    % first mouse ID
V.m2          = 3;    % second mouse ID
V.bs1         = 6;    % base session #1
V.bs2         = 7;    % base session #2
V.Diffs1      = [-5 -1 0 1 5];
V.Diffs2      = [-5 -1 0 1 4];             
V             = parseargs(V, varargin{:});  

LoadDir = '/home/mrule/Workspace2/PPC_data/';
fprintf(1,'Loading datasets from %s,\n',LoadDir);
fprintf(1,'please re-define LoadDir to match the path to the data files on your system.\n');
fprintf(1,'Datasets are available from the Harvey lab.\n');

figure(300); clf;
k_idxs = [2 4 5]; %indices for kinetic vars
labels = {'MAE Position (m)', 'MAE Velocity (m/s)','MAE View Angle ($^\circ$)'};
ryaxes = {[0.2 1], [0.05 0.3], [5 25]};

%% -- Panel (a) & (b): --
% Left 

% MR20190722 decay in performance over time has been pre-computed;
% (or can be re-computed using DecayRate_StatLin_kfoldCV if those
% files are not available).

loadname1 = ['m0' num2str(V.m1) '_bs' num2str(V.bs1) '_DecayMAE_StatLin_kfoldCV.mat'];
if exist([pwd '/Results/' loadname1],'file')==2 % load results
   O1 = load([pwd '/Results/' loadname1],'Output'); 
   O1 = O1.Output; 
else % perform computations 
   O1 = DecayRate_StatLin_kfoldCV('mouse',V.m1, 'BaseSession',V.bs1, 'Diffs',V.Diffs1, 'viz',0); 
end 

DeltaDays1 = O1.days - O1.days(3);  
sps        = [1 6 11]; 

for idx = 1:3 
    k = k_idxs(idx); 
    subplot(3,5,sps(idx)); % first mouse, forward position (y)
    yyaxis left; 
    plot1 = shadedErrorBar(DeltaDays1, O1.mean_MAE_pdiff(k,:), O1.std_MAE_pdiff(k,:),...
        'lineprops', {'b-o','markerfacecolor','b','MarkerSize',4,'LineWidth',2});
    plot1.Color(4) = 0.1; 
    plot1.MarkerFaceAlpha = 0.1;
    if idx==2
        ylabel('$\%$ Change vs. Same Day','Interpreter','Latex');
    end 
    
    yyaxis right; 
    plot2 = shadedErrorBar(DeltaDays1, O1.mean_MAE_abs(k,:), O1.std_MAE_abs(k,:),...
        'lineprops', {'r-o','markerfacecolor','r','MarkerSize',4,'LineWidth',2}); 
    plot2.Color(4) = 0.1; 
    plot2.MarkerFaceAlpha = 0.1;
    xlim([DeltaDays1(1) DeltaDays1(end)]); axis square; 
    ylabel(labels{idx},'Interpreter','Latex'); 
    ax = gca; ax.XTick = setdiff(DeltaDays1,[-1 1]);
    ylim(ryaxes{idx}); ax.YTick = ryaxes{idx}; 
    ax.YTickLabel = {num2str(ryaxes{idx}(1)) num2str(ryaxes{idx}(2))}; 
    
    if idx==1
        title(['Mouse ' num2str(V.m1)], 'Interpreter', 'Latex');
    elseif idx==3
        xlabel('$\Delta$ Days','Interpreter','Latex'); 
    end 
    box off;
end 

% Right
loadname2 = ['m0' num2str(V.m2) '_bs' num2str(V.bs2) ...
             '_DecayMAE_StatLin_kfoldCV.mat'];
if exist([pwd '/Results/' loadname2],'file')==2 % load results
   O2 = load([pwd '/Results/' loadname2],'Output'); 
   O2 = O2.Output; 
else % perform computations 
   O2 = DecayRate_StatLin_kfoldCV('mouse',V.m2, 'BaseSession',V.bs2, ...
                                  'Diffs',V.Diffs2, 'viz',0); 
end

DeltaDays2 = O2.days - O2.days(3); 
sps        = [2 7 12]; 

for idx = 1:3 
    k = k_idxs(idx); 
    subplot(3,5,sps(idx)); % first mouse, forward position (y)
    yyaxis left;
    plot1 = shadedErrorBar(DeltaDays2, O2.mean_MAE_pdiff(k,:), ...
            O2.std_MAE_pdiff(k,:), 'lineprops', ...
            {'b-o','markerfacecolor','b',...
            'MarkerSize',4,'LineWidth',2});
    plot1.Color(4) = 0.1; 
    plot1.MarkerFaceAlpha = 0.1;
    
    yyaxis right; 
    plot2 = shadedErrorBar(DeltaDays2, O2.mean_MAE_abs(k,:), ...
            O2.std_MAE_abs(k,:), 'lineprops', ...
            {'r-o','markerfacecolor','r',...
            'MarkerSize',4,'LineWidth',2}); 
    plot2.Color(4) = 0.1; 
    plot2.MarkerFaceAlpha = 0.1;
    xlim([DeltaDays2(1) DeltaDays2(end)]); axis square; 
    ax = gca; ax.XTick = setdiff(DeltaDays2, [-1 1]);
    ylim(ryaxes{idx}); ax.YTick = ryaxes{idx}; 
    ax.YTickLabel = {num2str(ryaxes{idx}(1)) num2str(ryaxes{idx}(2))}; 
    
    if idx==1
        title(['Mouse ' num2str(V.m2)], 'Interpreter', 'Latex');
    elseif idx==3
        xlabel('$\Delta$ Days','Interpreter','Latex'); 
    end 
    box off;
end 

%% -- Panel (c): --

% MR20190722 Looks like all of these are precomputed and stored
% in the Results directory. Which script produced them? 
% "bs6" and "bs7" might stand for "base session", and refer to the
% first session in a series of consecutive days? 
% `StatLin_Concat_match.m` might be the matfile that wrote these?

% Load data:
m04_avg_file    = 'm04_bs6_AvgBestPerfkCV_behavs364.mat';
m03_avg_file    = 'm03_bs7_AvgBestPerfkCV_behavs364.mat';
m04_concat_file = 'm04_bs6_ConcatMatch_behavs23564.mat';
m03_concat_file = 'm03_bs7_ConcatMatch_behavs23564.mat';

m04_avg    = load([pwd '/Results/' m04_avg_file],'O'); 
m04_avg    = m04_avg.O.avg_bestdecodeperf_CV;
m03_avg    = load([pwd '/Results/' m03_avg_file],'O'); 
m03_avg    = m03_avg.O.avg_bestdecodeperf_CV;
m04_concat = load([pwd '/Results/' m04_concat_file],'Output'); 
m04_concat = m04_concat.Output.MAE_cache;
m03_concat = load([pwd '/Results/' m03_concat_file],'Output'); 
m03_concat = m03_concat.Output.MAE_cache;

% Forward position:
C_m04_pos = m04_concat(2,:); 
C_m03_pos = m03_concat(2,:); 
I_m04_pos = m04_avg(1,:); 
I_m03_pos = m03_avg(1,:); 
% Forward speed:
C_m04_fs = m04_concat(4,:); 
C_m03_fs = m03_concat(4,:); 
I_m04_fs = m04_avg(2,:); 
I_m03_fs = m03_avg(2,:); 
% View angle:
C_m04_va = m04_concat(5,:); 
C_m03_va = m03_concat(5,:); 
I_m04_va = m04_avg(3,:); 
I_m03_va = m03_avg(3,:); 

subplot(3,5,3); % position 
errorbar([0.5 2],[mean(C_m04_pos) mean(C_m03_pos)],...
         [std(C_m04_pos) std(C_m03_pos)],'go'); 
hold on;
errorbar([1 2.5],[mean(I_m04_pos) mean(I_m03_pos)],...
         [std(I_m04_pos) std(I_m03_pos)],'ko','MarkerFaceColor','k');
xlim([0 3]); axis square; ax = gca;
ax.XTick = [0.5 1 2 2.5];
ax.XTickLabel = {'C' 'I' 'C' 'I'}; 
ylabel('MAE Position (m)','Interpreter','Latex'); box off;  
xlabel('$\:\:\: m4 \:\:\:\:\:\:\:\:\:\:\:\: m3 \:\:\:$', ...
       'Interpreter','Latex');
ylim(ryaxes{1}); ax.YTick = ryaxes{1}; 
ax.YTickLabel = {num2str(ryaxes{1}(1)) num2str(ryaxes{1}(2))}; 

subplot(3,5,8); %forward speed
errorbar([0.5 2],[mean(C_m04_fs) mean(C_m03_fs)],...
         [std(C_m04_fs) std(C_m03_fs)],'go'); 
hold on;
errorbar([1 2.5],[mean(I_m04_fs) mean(I_m03_fs)],...
         [std(I_m04_fs) std(I_m03_fs)],'ko','MarkerFaceColor','k');
xlim([0 3]); axis square; ax = gca;
ax.XTick = [0.5 1 2 2.5];
ax.XTickLabel = {'C' 'I' 'C' 'I'}; 
ylabel('MAE Speed (m/s)','Interpreter','Latex'); box off; 
xlabel('$\:\:\: m4 \:\:\:\:\:\:\:\:\:\:\:\: m3 \:\:\:$', ...
       'Interpreter','Latex');
ylim(ryaxes{2}); ax.YTick = ryaxes{2}; 
ax.YTickLabel = {num2str(ryaxes{2}(1)) num2str(ryaxes{2}(2))}; 

subplot(3,5,13); %view angle
errorbar([0.5 2],[mean(C_m04_va) mean(C_m03_va)],...
         [std(C_m04_va) std(C_m03_va)],'go'); 
hold on;
errorbar([1 2.5],[mean(I_m04_va) mean(I_m03_va)],...
         [std(I_m04_va) std(I_m03_va)],'ko','MarkerFaceColor','k');
xlim([0 3]); axis square; ax = gca;
ax.XTick = [0.5 1 2 2.5];
ax.XTickLabel = {'C' 'I' 'C' 'I'}; 
ylabel('MAE View Angle ($^\circ$)','Interpreter','Latex'); box off; 
xlabel('$\:\:\: m4 \:\:\:\:\:\:\:\:\:\:\:\: m3 \:\:\:$', ...
       'Interpreter','Latex');
ylim(ryaxes{3}); ax.YTick = ryaxes{3}; 
ax.YTickLabel = {num2str(ryaxes{3}(1)) num2str(ryaxes{3}(2))}; 


% -- Format font etc. -- 
a=findobj(gcf); alltext=findall(a,'Type','text');
allaxes=findall(a,'Type','axes');
set(allaxes,'FontName','Helvetica','FontWeight','Normal','FontSize',14);
set(alltext,'FontName','Helvetica','FontWeight','Normal','FontSize',16);

end %main

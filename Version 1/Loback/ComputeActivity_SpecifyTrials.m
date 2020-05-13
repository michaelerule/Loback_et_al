function [Output] = ComputeActivity_SpecifyTrials(varargin)
global LoadDir

%__________________________________________________________________________
% Description: Function called by StaticLinearModel_Concat.m
% Use: Can then specify trial type and remove activity associated
%      with NaN trials.
% Written 18 May, 2018 by AL
%__________________________________________________________________________

%% -- Initializations: -- 
V.DataDir = LoadDir;                   %Load directory 
V.mouse   = 3;                         %Mouse identity 
V.session = 23;                        %Note: "_s" != day necessarily 
V         = parseargs(V, varargin{:}); 

%% -- Main Processing: -- 
loadname = sprintf('%s/m%02d_s%02d.mat',V.DataDir,V.mouse,V.session);
session_obj = load(loadname,'session_obj');
session_obj = session_obj.session_obj;

% -- Get time series specifying whether the current time t is an
% inter-trial interval (ITI) time bin or not: -- 

itis = session_obj.timeSeries.virmen.data(9,:);
d    = diff(itis); 
trialends   = find(d==1);
trialstarts = find(d==-1)+1;
if trialends(1)<trialstarts(1)
   trialends   = trialends(2:end);
   trialstarts = trialstarts(1:end-1);
end 
if length(trialends)~=length(trialstarts)
   if length(trialends)>length(trialstarts)
      trialends = trialends(1:end-1);
   end 
   if length(trialstarts)>length(trialends)
      if trialends(end)<trialstarts(end), trialstarts = trialstarts(1:end-1); end 
   end 
end 

if length(trialends)~=length(trialstarts)
   error('Size mistmatch\n');
end 

trial_times = [trialstarts(:) trialends(:)];

%% -- Only Keep Correct Trials: -- 
correct_trials         = find(session_obj.trials.correct==1); 
trialtimes_correct_all = trial_times(correct_trials,:); 

%% -- Compute by Trial Type: --
%Note: '2' indicates black right-turn trial, and 
%      '3' indicates white left-turn trial
BR_trials = intersect(find(session_obj.trials.trialType==2), correct_trials);
WL_trials = intersect(find(session_obj.trials.trialType==3), correct_trials); 
trialtimes_correct_BR = trial_times(BR_trials,:); 
trialtimes_correct_WL = trial_times(WL_trials,:); 

%% -- Cache: -- 
Output.trialtimes_correct_all = trialtimes_correct_all; 
Output.trialtimes_correct_BR  = trialtimes_correct_BR;
Output.trialtimes_correct_WL  = trialtimes_correct_WL;

end 

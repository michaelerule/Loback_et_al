function [] = MAEvsN_static_ranked(varargin)
%__________________________________________________________________________
% Description: Prepares data for extrapolation to assess performance  
%              of the *static decoder* to address the question:
%              How precisely could we decode with a larger population?
% Note:        This code performs the necessary analyses to generate 
%              the results for Fig 2e, i.e., using a rank ordering of
%              neurons so that the "best" subset of K neurons is used
%              for underlying data for the GP extrapolation (of 
%              MAE vs. number of PPC neurons). 
% Note:        The ranking is done separately for each kinematic variable.
%
% Written 22 March 2019 by AL
%__________________________________________________________________________

%% -- Initializations: -- 
V.mouse       = 4;          % mouse ID (can specify multiple IDs) 
V.Nthresh     = 200;        % choose all sessions that have >= this 
                            % threshold of high-confidence neurons           
V.behavvars   = [3 6 4];    % 3=forward position (units: m); 6=forward speed; 4=view angle
V.trial_type  = 'all';      % 2=black right-turn trial;
                            % 3=white left-turn trial;
                            % 'all' specifies both trial types
V.nxval       = 10;         % # of folds for cross-validation  
V             = parseargs(V, varargin{:});

Nvec = 5:5:250; 

%% -- Position: --
% Get ranking of "best" neurons and MAE vs. N for ranked neuron subsets
[ranking_pos, bestEV_pos, bestMAE_pos] = rankNeurons_explainedVar('mouse',V.mouse, ...
                                         'Nthresh',V.Nthresh, ...
                                         'behavvar',V.behavvars(1), ...
                                         'trial_type',V.trial_type, ...
                                         'nxval',V.nxval); %#ok                                      

%% -- Speed: --
[ranking_vel, bestEV_vel, bestMAE_vel] = rankNeurons_explainedVar('mouse',V.mouse, ...
                                         'Nthresh',V.Nthresh, ...
                                         'behavvar',V.behavvars(2), ...
                                         'trial_type',V.trial_type, ...
                                         'nxval',V.nxval); %#ok 

%% -- View Angle: --
[ranking_va, bestEV_va, bestMAE_va] = rankNeurons_explainedVar('mouse',V.mouse, ...
                                         'Nthresh',V.Nthresh, ...
                                         'behavvar',V.behavvars(3), ...
                                         'trial_type',V.trial_type, ...
                                         'nxval',V.nxval);%#ok 

%% -- Process Results for Exporting to Python: -- 

R_matrix = [Nvec(:) bestMAE_pos(Nvec)' bestMAE_vel(Nvec)' bestMAE_va(Nvec)'];

%% -- Save: -- 
savename = ['MAEvsN_static_ranked_m0' num2str(V.mouse) '_5.csv'];

% Put into .csv format to be read into Python: 
csvwrite(savename,R_matrix);

% Save in .mat format for backup:
save(strrep(savename,'.csv','.mat'),'R_matrix');

end
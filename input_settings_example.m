% Input settings
%-----------------------------------------------
% Note: All variables given/calculated within this file are accessible to the main program

% Time(s) at which the metric is determined
time_tag='tapp';    % Appearance times
%time_tag='last';    % Last time in the series (typically a steady state)
%time_tag='sel';     % Selected times (below)
%time_tag='all';     % All times

% Specific times at which to determine the metric, if the corresponding time option is chosen
if strcmp(time_tag,'sel')
    times_for_ddDp_ratio=90:10:150; % min
end

% In case appearance times are (roughly!) assessed, give here the definition of appearance as
% the relative increase compared to max. concentraton in bin; 0.5 is most often used
referenceC_fraction=0.5;

% Option to include weights in fits; comment out if not used
% Weights can be applied e.g. to assess effects of uncertainties in measured concentrations
%array_weight=nan;   % Weight for each size bin (can also be moved to the data input .m file, if more convenient)
% One can also give extra weight for closest bins when assessing the slope of the size distribution around given bin
%extra_weight=100;   % Weight (relative to default weight, i.e. 1 if no other weights are given)
%nadj_weight=2;      % Number of bins around fitted bin for which extra weight is used (on both sides, i.e. total 2*nadj_weight+1 bins)

% Threshold value for a low enough value for the metric (below which continuous models are assumed valid)
% This should be a small number, e.g. 0.05 (5%) is okay
small_th=0.05;

% Show the fits to examine how good they are (should alweays be done)
l_show_fits=1;  % Recommended; should be disabled only for speed after the fits have once been examined
l_pause=0;      % Pause between figures of different fits to see them better

% Figure styles
% Different data can be plotted in same figure by keeping the figure window open and changing the styles and/or number for color
Lstyle='-';
Mstyle='o';
nrun=1; % For colors

% Limits for surface plots of dNdlogDp (i.e. lower and upper limits for colormap)
% Note that this is given as the base-10 logarithm of dNdlogDp
lim_for_caxis=[-3 6]; % log10(cm^-3/log10Dp)

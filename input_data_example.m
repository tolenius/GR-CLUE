% Input data for nanoparticle size distribution
%-----------------------------------------------
% Note: All variables given/calculated within this file are accessible to the main program
%
% The size distribution data must be in the following format:
% dNdlogDp:         matrix with rows and columns corresponding to times and size bins, respectively; unit cm^-3/log10Dp
% time_vector:      column vector; leave out if only a steady-state distribution is given; unit min
% Dp_mean_bin:      row vector of bin mean mobility diameters; unit m
% Dp_limits_bin:    row vector of bin limits' mobility diameters; unit m
%
%-----------------------------------------------
%
% This example file demonstrates the usage by applying synthetic, model-generated data for the particle size distribution
% All calculations below can be modified according to the format of the applied data, as long as the required variables (above)
% are obtained
%
%-----------------------------------------------

% Load some constants needed for conversions etc.
get_constants;

%-----------------------------------------------

% Example file: total concentrations in given size classes; assumed format:
% Header:           "%d_mob" followed by mobility diameters at bin limits (here nm), e.g. "%d_mob 0.9 1.2 ..."
% 1st column:       Time since the start of the experiment (here h)
% Other columns:    Concentration in each bin (here cm^-3)

fn='time_conc_example.txt';

% Read in particle sizes from the header, and determine mean bin sizes and bin limits
fid=fopen(fn,'r');
s1=textscan(fid,'%s',1,'delimiter','\n');
fclose(fid);
s1=strtrim(s1);
s2=regexp(s1{1,1}{1},'\s+','split');
% Ignore the first character/word in the beginning of the header, and get size bin limits from
% the rest of the header elements
if strcmp(s2{1},'%d_mob')
    Dp_limits_bin=cellfun(@(x) str2double(x)*1e-9, s2(2:end)); % nm -> m
else
    error([s2{1},': unexpected format'])
end
% Assume linear midpoints for bin mean sizes (for simplicity; for the smallest nm-sizes it doesn't matter anyway)
% Mean sizes are used for the size distribution function; bin limits are used only for figure axes limits, legends etc.
Dp_mean_bin=(Dp_limits_bin(1:end-1)+Dp_limits_bin(2:end))/2;

% Load the data
data=load(fn);
time_vector=data(:,1)*60;   % h -> min
dN=data(:,2:end);           % cm^-3

% Divide dN by dlogDp to get dNdlogDp
dNdlogDp=bsxfun(@rdivide,dN,diff(log10(Dp_limits_bin)));

% Bins included for representing the distribution/applying fits (should be all, unless there's something wrong with some of the data)
% By default, all bins are included if no min/max indices are given here
nbin_min=1;
nbin_max=length(Dp_mean_bin)-6; % E.g. here the last bins are empty/uninteresting/may involve simulation boundary artefacts

% Bins analysed / plotted (if e.g. only a given size interval is of interest)
% By default, all bins are included if no min/max indices are given here
nbin_start=1;
nbin_end=length(Dp_mean_bin)-8; % E.g. here, leave out some last bins that have less neighbor bins for useful fits
%nbin_start=find(Dp_limits_bin >= 3e-9, 1, 'first'); % Or focus only on specific sizes
%nbin_end=find(Dp_limits_bin <= 4e-9, 1, 'last');

%-----------------------------------------------

% Option to give pre-calculated appearance times (in case the GR metric is determined at such times); otherwise a rough assessment is done
% (i.e. this option can be used for defining appearance times according to a preferred approach)
%apptime_bin=nan;

% Assumed molecular properties specific for the data set (used only for Dp / molec. conversions)
% One or more (for average properties) compounds can be assumed
mv=[500 300]*amu2kg;    % assumed ELVOC/LVOC masses; amu -> kg
rho=1400;               % kg m^-3

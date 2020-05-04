%% Assign electrodes to an ROI
% specifically to schaefer parcel, and then to a yeo 7 system


% @author JStiso jeni.stiso@gmail.com

clear
clc
close all

warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/BrainNetViewer_20171031/'))

%% load
% node data

top_dir = '/Volumes/bassett-data/Jeni/RAM/';
parcel_dir = '/Users/stiso/Documents/parcellations/';
win = 1;
detector = '';
betas = readtable([top_dir, 'group_analysis/win_',  num2str(win), '/node_stats', detector, '.csv']);
beta_names = betas.Properties.VariableNames(cellfun(@(x) contains(x, 'beta'), betas.Properties.VariableNames));
band_measures = unique(betas.band_measure);
parc = readtable([parcel_dir, 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv']);

% convert coords to numbers
betas.y = cellfun(@(x) str2double(x), betas.y);
betas.x = cellfun(@(x) str2double(x), betas.x);
betas.z = cellfun(@(x) str2double(x), betas.z);
        
% parse labels from table
label_chunks = cellfun(@(x) split(x, '_'), parc.ROIName, 'UniformOutput', false);
parc.sys = cellfun(@(x) x{3}, label_chunks, 'UniformOutput', false);
parc.reg = cellfun(@(x) x{4}, label_chunks, 'UniformOutput', false);
parc.hem = cellfun(@(x) x{2}, label_chunks, 'UniformOutput', false);

% RAS -> XYZ mapping: RAS = XYZ (confirmed by clicking through plots)

%% Get the labels

% get a subset of the data that only gives you uique identifiers
contacts = betas(:,[2, 19:end]);
contacts = unique(contacts, 'stable');

% get label for every electrode
% this is redunadant, could make more efficient
for e = 1:size(contacts,1)
   dists = rowfun(@(~, ~, x, y, z, ~, ~, ~) pdist([x y z; contacts.x(e), contacts.y(e), contacts.z(e)], 'euclidean'), parc, 'OutputVariableName', 'd');
   [min_dist, idx] = min(dists.d);
   contacts.parc(e) = parc.reg(idx);
   contacts.sys(e) = parc.sys(idx);
   contacts.hem(e) = parc.hem(idx);
   contacts.min_d(e) = min_dist;
end

% save
writetable(contacts, [top_dir, 'group_analysis/win_',  num2str(win), '/contact_sys', detector, '.csv'])

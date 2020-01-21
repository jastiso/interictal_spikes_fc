%% Get functional connectivity

clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath('/Users/stiso/Documents/MATLAB/arfit/')

%%

% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

% constants
load([top_dir, 'bad_datasets.mat'])

releases = ['1', '2', '3'];
table_names = [{'subj'}, {'exper'}, {'sess'}, {'time'}, {'power'}, {'fc_measure'}, {'band'}, {'str'},...
    {'str_soz'}, {'str_not_soz'}, {'str_spike'}, {'str_not_spike'}, {'elec'}, {'region'}, {'elec_in_soz'},...
    {'elec_has_spike'}, {'spike_num'}, {'age'}, {'gender'}, {'race'}, {'hand'}, {'x'}, {'y'}, {'z'}];
errors = struct('files', [], 'message', []);
%bands
freqs = unique(round(logspace(log10(4),log10(150),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70; 71, 150];
band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}, {'hg'}];

%fc measures
measure_names = [{'coh'}, {'plv'}, {'aec'}, {'xcorr'}, {'ar'}, {'pac'}];
%parameters
pmin = 1; pmax = 1; % order for AR model
spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds

parfor r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % for catching errors
    errors = struct('files', [], 'message', []);
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    functional_connectivity(protocols, release, release_dir, top_dir, errors, table_names, freqs, bands, ...
        band_names, measure_names, pmin, pmax, spike_win, win_length, bad_datasets)
end


% get and concatenate errors and warnings
errors_all = struct('files', [], 'message', []);
for r = 1:numel(releases)
    release = releases(r);
    release_dir = [top_dir, 'release', release '/'];
    
    curr_err = load([release_dir, 'protocols/fc_errors.mat']);
    
    errors_all= [errors_all, curr_err.errors];
end
save([top_dir, 'fc_errors.mat'], 'errors_all')


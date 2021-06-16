%% Count oscillations

clear
clc
close all
warning ON

restoredefaultpath
addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))

% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '';
win_length_ir = 400; % in ms
step = 50; %in ms
filter = 0; % with ot without lowpass alaising filter

perc_spike = [];
cnt = 1;

% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '_delphos_auto';

%bands
freqs = unique(round(logspace(log10(4),log10(150),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70];
band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}];
% constants
nBand = size(bands,1);
all_osc = {};
cnt = 1;

% subjects not to use
load([top_dir, 'bad_datasets.mat'])
errors = struct('files', [], 'message', []);

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    % just for testing that qsub fuction will work
    for p = 1:numel(protocols)
        protocol = protocols{p};
        
        % get global info struct
        fname = [release_dir 'protocols/', protocol, '.json'];
        fid = fopen(fname);
        raw = fread(fid);
        str = char(raw');
        fclose(fid);
        info = jsondecode(str);
        info = info.protocols.r1;
        
        % get subjects
        subjects = fields(info.subjects);
        parfor s = 1:numel(subjects)
            subj = subjects{s};

            get_win_osc(bands, detector, win_length, win_length_ir, step, filter, top_dir, release, protocol, subj)
        end
    end
end

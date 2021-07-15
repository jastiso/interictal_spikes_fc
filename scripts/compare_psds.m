%% Get PSDs for IED and normal windows

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
filter = 0; % without lowpass alaising filter

all_psd = [];
all_ied = [];
all_aper = [];
all_aperied = [];
cnt = 1;

% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '';

%bands
freqs = unique(round(logspace(log10(4),log10(150),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70];
band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}];
% constants
nBand = size(bands,1);
all_osc = {};
nFreq = 34; % number of frequencies in the IRASA


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
            
            psds(bands, detector, win_length, win_length_ir, spike_win, step, filter, top_dir, release, protocol, subj, nFreq)
        end
    end
end


%% Make plots

% load the data
cnt = 1;
group_aper = [];
group_psd = [];
group_iedaper = [];
group_iedpsd = [];
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
        for s = 1:numel(subjects)
            subj = subjects{s};
            if exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/slopes.mat'], 'file') && exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/psds.mat'], 'file')
                fprintf('Subj %s\n', subj)
                load([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/slopes.mat'], 'all_aper', 'all_aperied')
                load([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/psds.mat'], 'all_psd', 'all_ied')
                if ~isempty(all_aper)
                    try
                    group_aper(cnt) = mean(all_aper);
                    group_psd(cnt,:) = mean(all_psd,1);
                    group_iedaper(cnt) = mean(all_aperied);
                    group_iedpsd(cnt,:) = mean(all_ied,1);
                    cnt = cnt + 1;
                    catch
                        if ~exist(all_aperied)
                            %delete([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/psds.mat']);
                            %delete([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/slopes.mat'])
                        end
                    end
                else
                    %delete([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/psds.mat']);
                    %delete([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/slopes.mat'])
                end
            end
        end
    end
end
figure(1); clf
histogram(group_aper, 10); hold on
histogram(group_iedaper,10)
saveas(gca,[top_dir, 'img/slopes.png'], 'png')

figure(2); clf
for i = 1:size(group_psd,1)
    lh1 = plot(linspace(2,70,34),log(group_psd(i,:)), 'color', rgb('slategray')); hold on % freqs taken from get_IRASA_spec function
    lh2 = plot(linspace(2,70,34),log(group_iedpsd(i,:)), 'color', rgb('sandybrown'));
    lh1.Color = [lh1.Color 0.5];
    lh2.Color = [lh2.Color 0.5];
end
saveas(gca,[top_dir, 'img/psds.png'], 'png')


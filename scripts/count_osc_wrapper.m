%% Count oscillations

clear
clc
close all
warning ON

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

%%
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
            
            get_win_osc(bands, detector, win_length, win_length_ir, spike_win, step, filter, top_dir, release, protocol, subj)
        end
    end
end


%% Make plots
% want the average and variance of the percent of contacts with osc in
% each window

vari = [];
mu = [];
cnto = 1;
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
            release_dir = [top_dir, 'release', release '/'];
            
            % get global info struct
            fname = [release_dir 'protocols/', protocol, '.json'];
            fid = fopen(fname);
            raw = fread(fid);
            str = char(raw');
            fclose(fid);
            info = jsondecode(str);
            eval(['info = info.protocols.', protocol,';']);
            
            % subjects not to use
            load([top_dir, 'bad_datasets.mat'])
            errors = struct('files', [], 'message', []);
            
            % make subject directory
            subj_dir = [top_dir, 'FC/release',release, '/', protocol, '/', subj, '/'];
            if ~exist(subj_dir, 'dir')
                mkdir(subj_dir);
            end
            
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
                        
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1:numel(experiments)
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1:numel(sessions)
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    data_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    save_dir = [top_dir, 'FC/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    img_dir = [top_dir, 'img/FC/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    if ~exist(img_dir, 'dir')
                        mkdir(img_dir);
                    end
                    if ~exist(save_dir, 'dir')
                        mkdir(save_dir);
                    end
                    if ~exist([subj_dir, 'win_', num2str(win_length), '/'], 'dir')
                        mkdir([subj_dir, 'win_', num2str(win_length), '/']);
                    end
                    
                    if exist([data_dir, 'osc.mat'], 'file')
                        fprintf('\n******************************************Finishing oscillation counts for subject %s...\n', subj)
                        load([data_dir, 'osc.mat'])
                        nElec = size(osc,3);
                        
                        % get mean and var of the percent
                        for k = 1:nBand
                            mu(k,cnto) = mean(sum(osc(:,k,:),3)/nElec*100);
                            vari(k,cnto) = var(sum(osc(:,k,:),3)/nElec*100);
                        end
                        cnto = cnto + 1;
                    end
                end
            end
        end
    end
end


% plot
colors = ["#516888", "#FBE697", "#F3AE6D", "#C9DACA"];
figure(1); clf
for k = 1:nBand
    subplot(1,2,1)
    histogram((mu(5-k,:)), 'FaceColor', colors(5-k), 'EdgeColor', 'none', 'FaceAlpha', .9); hold on
    title('Mean Percentage of Contacts with Oscillation')
    subplot(1,2,2)
    histogram((vari(5-k,:)), 'FaceColor', colors(5-k), 'EdgeColor', 'none', 'FaceAlpha', .9); hold on
    title('Variance of Percentage of Contacts with Oscillations')
    legend(band_names)
end
saveas(gca, [img_dir, 'osc_hist.pdf'], 'pdf')

%% For subjects that don't have SOZ marked, look in other datasets
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT/'))

%%

% global variables and packages
top_dir = '/Volumes/bassett-data/Jeni/RAM/';
eval(['cd ', top_dir])

% parameters for eliminated spikes
min_chan = 3; % minimum number of channels that need to be recruited
win = 0.05; % size of the window to look for the minimum number of channels, in seconds
releases = ['1', '2', '3'];
discharge_tol=0.005; % taken from spike function
n_spikes = [];
% for catching errors
errors = struct('files', [], 'message', []);

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    for p = 1:numel(protocols)
        protocol = protocols{p};
        
        % get global info struct
        fname = [release_dir 'protocols/', protocol, '.json'];
        fid = fopen(fname);
        raw = fread(fid);
        str = char(raw');
        fclose(fid);
        info = jsondecode(str);
        eval(['info = info.protocols.', protocol,';']);
        
        % get subjects
        subjects = fields(info.subjects);
        for s = 1:numel(subjects)
            subj = subjects{s};
            
            % save command window
            %clc
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            %eval(['diary ', [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/log.txt']]);
            
            fprintf('\n******************************************\nStarting IED detection for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1:numel(experiments)
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1:numel(sessions)
                    %initialize
                    out_clean = [];
                    marker = [];
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    save_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    img_dir = [top_dir, 'img/spikes/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    if ~exist(img_dir, 'dir')
                        mkdir(img_dir);
                    end
                    
                    % check if this subect has clean data
                    if exist([save_dir, 'data_clean.mat'], 'file')
                        load([save_dir, 'channel_info.mat'])
                        orig_soz = soz;
                        orig_ictal = interictal_cont;
                         % if they dont have SOZ, look in other data for
                            % this subject
                            if isempty(orig_soz) || isempty(orig_ictal)
                                warning('The interictal and SOZ channels were not marked for this subject...checking other datasets')
                                
                                % get all paths to folders with
                                % "channel_info.mat"
                                save_flag = 0;
                                files = dir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/**/channel_info.mat']);
                                for i = 1:numel(files)
                                   try
                                       load([files(i).folder, '/channel_info.mat'], 'soz', 'interictal_cont')
                                       if isempty(orig_soz) && ~isempty(soz)
                                         save_flag = 1;
                                       else
                                           soz = orig_soz;
                                       end
                                       if isempty(orig_ictal) && ~isempty(interictal_cont)
                                         save_flag = 1;
                                       else
                                           interictal_cont = orig_ictal;
                                       end
                                       if save_flag
                                           print('Found a replacement file\n')
                                           save([save_dir, 'channel_info.mat'], 'channel_info', 'soz', 'interictal_cont', 'labels', 'chann_idx', 'regions') 
                                       end
                                   catch 
                                   end
                                end
                            end

                    end
                end
            end
        end
    end
end
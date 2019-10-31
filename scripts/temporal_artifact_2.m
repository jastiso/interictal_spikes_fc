%% Find Temporally SPecific Artifacts
% ~300 subjects, for a total of 1073 sessions of data.
% This script saves cleaned task-free data for each session

clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

%%

% global variables and packages
top_dir = '/Volumes/bassett-data/Jeni/RAM/';
eval(['cd ', top_dir])

% for marking artifacts
thr = 30000; % how big the derivative needs to be to flag it
releases = ['1', '2', '3'];
n_art = [];

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % for catching errors
    errors = struct('files', [], 'message', []);
    warnings = struct('files', [], 'message', []);
    
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
            
            fprintf('\n******************************************\nStarting artifact detection for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1:numel(experiments)
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1:numel(sessions)
                    %initialize
                    
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    save_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    img_dir = [top_dir, 'img/artifact/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    if ~exist(img_dir, 'dir')
                        mkdir(img_dir);
                    end
                    
                    % check if this subect has clean data
                    if exist([save_dir, 'data_clean.mat'], 'file')
                        load([save_dir, 'data_clean.mat'])
                        load([save_dir, 'header.mat'])
                        load([save_dir, 'channel_info.mat'])
                        
                        nTrial = numel(ft_data.trial);
                        nElec = numel(ft_data.label);
                        half = floor(nElec/2);
                        for j = 1:nTrial
                            curr = ft_data.trial{j};
                            curr_time = ft_data.time{j};
                            try
                                % get derivative of time series
                                td = diff(curr')'./diff(curr_time);
                                td = [ones(nElec,1), td]; % add a one to counteract size change from diff
                                
                                %check for 0 or very large
                                artifact_sharp = (td > thr );
                                artifact_flat = (td == 0);
                                artifacts = artifact_sharp | artifact_flat;
                                artifact_type = artifacts + artifact_sharp; % sharp derivative marked as 2, zeros as 1
                                n_chan = sum(artifacts);
                                % get only those in 1/3 of channels
                                artifact_idx = n_chan >= half;
                                artifact_type = artifact_type(artifact_idx);
                                
                                fprintf('%d artifacts found\n', sum(artifact_idx));
                                n_art{1,end+1} = subj; n_art{2,end} = exper;
                                n_art{3,end} = sess; n_art{4,end} = j;
                                n_art{5,end} = sum(artifact_idx);
                                
                                save([save_dir, 'artifact.mat'], 'artifact_idx', 'artifact_type')
                                
                                % plot
                                if any(artifact_idx)
                                    art = find(artifact_idx);
                                    art = art(1);
                                    if art > 501 && art < (size(curr,2) - 500)
                                        plot_lfp(curr(:,(art-500):(art+500)), header.sample_rate)
                                        close
                                    end
                                    saveas(gca, [img_dir, 'artifact.png'], 'png')
                                end
                                
                            catch ME
                                errors(end+1).files = [subj, '_', exper, '_', sess];
                                errors(end).message = ME.message;
                            end
                            
                        end
                        
                    end
                end
            end
        end
    end
    
end

errors = errors(2:end);
save([top_dir, 'n_artifacts.mat'], 'n_art')
save([top_dir, 'artifact_errors.mat'], 'errors')

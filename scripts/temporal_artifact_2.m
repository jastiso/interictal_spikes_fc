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
deriv_thr = 30000; % how big the derivative needs to be to flag it
std_thr = 300;
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
                        artifact_all = struct('idx', [], 'type', []);
                        for j = 1:nTrial
                            curr = ft_data.trial{j};
                            curr_time = ft_data.time{j};
                            try
                                % get derivative of time series
                                td = diff(curr')'./diff(curr_time);
                                td = [ones(nElec,1), td]; % add a placeholder to counteract size change from diff
                                % get variance (for finding flatlining)
                                sigma = std(td);
                                
                                %check for very large derivative
                                n_chan = sum(td > deriv_thr );
                                % make sure its present in at least half of
                                % electrodes
                                artifact_sharp = n_chan >= half;
                                % check for flatlining
                                artifact_flat = (sigma < std_thr);
                                artifact_flat(1) = 0; % accounting for the fact that the first time point is constant
                                % flatlineing artifacts have to have the
                                % next 5 seconds removed (5s comes from
                                % spike detection algorithm widnow size)
                                % flatlines create lots of spurious spikes
                                % because spikes are identified in
                                % comparison to the background
                                win = 5*header.sample_rate;
                                indices = find(artifact_flat);
                                prev_idx = 0;
                                if sum(artifact_flat > 0)
                                    for i = 1:sum(artifact_flat)
                                        if (indices(i) + win - 1) < numel(artifact_flat)
                                            artifact_flat(indices(i):(indices(i) + win - 1)) = 1;
                                        else
                                            artifact_flat(indices(i):end) = 1;
                                            break
                                        end
                                    end
                                end
                                
                                artifact_idx = artifact_sharp | artifact_flat;
                                artifact_type = artifact_idx + artifact_sharp; % sharp derivative marked as 2, zeros as 1
                                artifact_all(j).idx = artifact_idx;
                                artifact_all(j).type = artifact_type;

                                fprintf('%d sharp artifacts found and %d flat artifacts found\n', sum(artifact_sharp), sum(artifact_flat));
                                n_art{1,end+1} = subj; n_art{2,end} = exper;
                                n_art{3,end} = sess; n_art{4,end} = j;
                                n_art{5,end} = sum(artifact_sharp); 
                                n_art{6,end} = sum(artifact_flat);

                                
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
                        save([save_dir, 'artifact.mat'], 'artifact_all')
                        
                    end
                end
            end
        end
    end
    
end

errors = errors(2:end);
save([top_dir, 'n_artifacts.mat'], 'n_art')
save([top_dir, 'artifact_errors.mat'], 'errors')

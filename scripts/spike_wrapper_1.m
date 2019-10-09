%% Find Spikes
%

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

% parameters for eliminated spikes
min_chan = 3; % minimum number of channels that need to be recruited
win = 0.15; % size of the window to look for the minimum number of channels, in seconds
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
        
        % directories
        metadata_dir = dir([top_dir, 'release', release '/Release_Metadata*']);
        metadata_dir = [top_dir, 'release', release '/', metadata_dir.name, '/'];
        
        % get subjects
        subjects = fields(info.subjects);
        for s = 1:numel(subjects)
            subj = subjects{s};
            
            % save command window
            clc
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            %eval(['diary ', [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/log.txt']]);
            
            fprintf('******************************************\nStarting IED detection for subject %s...\n', subj)
            
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
                    img_dir = [top_dir, 'img/diagnostics/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    
                    % check if this subect has clean data
                    if exist([save_dir, 'data_clean.mat'], 'file')
                        load([save_dir, 'data_clean.mat'])
                        load([save_dir, 'header.mat'])
                        load([save_dir, 'channel_info.mat'])
                        
                        % combine soz and interictal
                        all_soz = unique([cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), soz, 'UniformOutput', false)); cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), interictal_cont, 'UniformOutput', false))]);
                        nChan = numel(ft_data.label);
                        
                        if isempty(all_soz)
                            warning('The interictal and SOZ channels were not marked for this subject')
                        end
                        
                        nTrial = numel(ft_data.trial);
                        for j = 1:nTrial
                            % run detection alg
                            try
                                [out,MARKER] = ...
                                    spike_detector_hilbert_v16_byISARG(ft_data.trial{j}', header.sample_rate);
                                
                                % eliminate some spikes
                                include_length = win;%.300*MARKER.fs;
                                nSamp = size(MARKER.d,1);
                                nSpike = numel(out.pos);
                                kept_spike = false(size(out.pos));
                                
                                % for each spike, check if there are other
                                % spikes within WIN samples in at least
                                % MIN_CHAN soz or interictal channels
                                for i = 1:nSpike
                                    curr_pos = out.pos(i);
                                    curr_chan = out.chan(i);
                                    
                                    if kept_spike(i) == 0
                                        win_spike = (out.pos > curr_pos & out.pos < curr_pos + win);
                                        win_chan = out.chan(win_spike);
                                        
                                        if ~isempty(soz)
                                            if sum(intersect(win_chan,all_soz)) >= min_chan
                                                kept_spike(win_spike) = true;
                                            end
                                        else
                                            if sum(unique(win_chan)) >= min_chan+1 % if not SOZ marked have slightly stricter cutoff
                                                kept_spike(win_spike) = true;
                                            end
                                        end
                                    end
                                end
                                % slect only good spikes
                                out_clean.pos = out.pos(kept_spike);
                                out_clean.dur = out.dur(kept_spike);
                                out_clean.chan = out.chan(kept_spike);
                                out_clean.weight = out.weight(kept_spike);
                                out_clean.con = out.con(kept_spike);
                                n_spikes = [n_spikes; sum(kept_spike)];
                                % get new M
                                m_clean = zeros(size(MARKER.M));
                                for i=1:size(out_clean.pos,1)
                                    m_clean(round(out_clean.pos(i)*MARKER.fs:out.pos(i)*MARKER.fs+discharge_tol*MARKER.fs),...
                                        out_clean.chan(i))=out_clean.con(i);
                                end
                                MARKER.m_clean = m_clean;
                                
                                
                                % save
                                save([save_dir, 'spike_info.mat'], 'out_clean', 'MARKER');
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

% remove empty entry
errors = errors(2:end);
save([release_dir, 'spike_errors.mat'], 'errors');

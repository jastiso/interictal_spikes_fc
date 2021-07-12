function [psds, aper] = psds(bands, detector, win_length, win_length_ir, spike_win, step, filter, top_dir, release, protocol, subj)
% constants
nBand = size(bands,1);
cnt = 1;
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

fprintf('\n******************************************\nStarting psd comparison for subject %s...\n', subj)

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
        
        if exist([data_dir, 'data_clean.mat'], 'file') && (exist([data_dir, 'spike_info', detector, '.mat'], 'file') || exist([data_dir, 'spike_info', num2str(win), '.mat'], 'file'))
            load([data_dir, 'data_clean.mat'])
            load([data_dir, 'header.mat'])
            load([data_dir, 'channel_info.mat'])
            if strcmp(detector, '')
                load([data_dir, 'spike_info_', num2str(spike_win), '.mat'])
            else
                load([data_dir, 'spike_info', detector, '.mat'])
            end
            load([data_dir, 'artifact.mat'])
            load([data_dir, 'demographics.mat'])
            %try
            % check if this subect has clean data
            reject = zeros(numel(ft_data.trial),1);
            for i = 1:numel(ft_data.trial)
                curr_ext = [subj, '_' exper, '_', sess, '_', num2str(i)];
                reject(i) = any(strcmp(curr_ext, bad_datasets));
            end
            fprintf('\nRejected %d datasets\n', sum(reject))
            
            ft_data.trial = ft_data.trial(~reject);
            ft_data.time = ft_data.time(~reject);
            out_clean = out_clean(~reject);
            artifact_all = artifact_all(~reject);
            ft_data.sampleinfo = ft_data.sampleinfo(~reject,:);
            
            if ~isempty(ft_data.trial)
                
                nElec = numel(ft_data.label);
                
                psds = zeros(nElec,nFreq);
                aper = zeros(nElec,1);
                for j = 1:nElec
                    spec = get_IRASA_spec(ft_data.trial{1}(j,:), 1, numel(ft_data.trial{1}(j,:)), ft_data.fsample, win_length_ir, step, filter);
                    psds(j,:) = mean(spec.mixd,2);
                    x = spec.freq;
                    y = mean(spec.frac,2);
                    P = polyfit(log10(x),log10(y),1);
                    aper(j) = P(1);
                end
                curr_clean = squeeze(mean(psds,1));
                all_psd(cnt,:) = curr_clean;
                all_aper(cnt) = mean(mean(aper));
                cnt = cnt + 1;
                %                             catch
                %                                 fprintf("\nError for subj %s", subj)
                %                             end
            end
        end
    end
end

save([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/slopes.mat'], 'all_aper')
save([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/psds.mat'], 'all_psd')
end


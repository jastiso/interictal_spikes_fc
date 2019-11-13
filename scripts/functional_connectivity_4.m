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
table_names = [{'subj'}, {'exper'}, {'sess'}, {'time'}, {'power'}, {'fc_measure'}, {'band'}, {'str_soz'}...
                {'str_not_soz'}, {'str_spike'}, {'str_not_spike'}, {'elec'}, {'region'}, {'elec_in_soz'},...
                {'elec_has_spike'}, {'spike_num'}, {'age'}, {'gender'}, {'race'}, {'hand'}];
%bands
freqs = unique(round(logspace(log10(4),log10(150),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70; 71, 150];
band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}, {'hg'}];
bands_cf = [1, 5; 2, 5; 3, 5]; % gives indices in bands variable for phase and amp bands
nBand = size(bands,1);
nBandCF = size(bands_cf,1);
%fc measures
measure_names = [{'coh'}, {'plv'}, {'aec'}, {'xcorr'}, {'ar'}, {'pac'}];
nMeasures = numel(measure_names);
%parameters
pmin = 1; pmax = 1; % order for AR model
spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds

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
            
            % initialize table
            power_vars = cellfun(@(x) ['power_', x], band_names, 'UniformOutput', false);
            coh_vars = cellfun(@(x) ['coh_', x], band_names, 'UniformOutput', false);
            plv_vars = cellfun(@(x) ['plv_', x], band_names, 'UniformOutput', false);
            aec_vars = cellfun(@(x) ['aec_', x], band_names, 'UniformOutput', false);

            fc_table = cell2table(cell(0,20), 'VariableNames', table_names);
            
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            
            fprintf('\n******************************************\nStarting functional connectivity for subject %s...\n', subj)
            
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
                    
                    if exist([data_dir, 'data_clean.mat'], 'file')
                        load([data_dir, 'data_clean.mat'])
                        load([data_dir, 'header.mat'])
                        load([data_dir, 'channel_info.mat'])
                        load([data_dir, 'spike_info_', num2str(spike_win), '.mat'])
                        load([data_dir, 'artifact.mat'])
                        load([data_dir, 'demographics.mat'])
                        
                        % check if this subect has clean data
                        reject = zeros(numel(ft_data.trial),1);
                        for i = 1:numel(ft_data.trial)
                            curr_ext = [subj, '_' exper, '_', sess, '_', num2str(i)];
                            reject(i) = any(strcmp(curr_ext, bad_datasets));
                        end
                        ft_data.trial = ft_data.trial(~reject);
                        ft_data.time = ft_data.time(~reject);
                        %ft_data.timeinfo = ft_data.timeinfo(~reject,:);
                        
                        if ~isempty(ft_data.trial)
                            tic
                           
                            nElec = numel(ft_data.label);
                            nPair = (nElec^2-nElec)/2;
                            
                            % constants
                            upper_tri = reshape(triu(true(nElec),1),[],1);
                            
                            % get window start times
                            trl = [];
                            spike_idx = [];
                            spike_num = [];
                            spike_chans = [];
                            time_vec = [];
                            cnt = 1;
                            for i = 1:numel(ft_data.trial)
                                idx = 1;
                                curr_data = ft_data.trial{i};
                                curr_spike = out_clean(i);
                                curr_artifact = artifact_all(i);
                                dur = size(curr_data,2);
                                trl_offset = ft_data.sampleinfo(i,1);
                                while (idx + (win_length*header.sample_rate)) <= dur
                                    st = round(idx); % gets rid of scientific notation
                                    en = round(st + (win_length*header.sample_rate));
                                    st_ms = st/header.sample_rate;
                                    en_ms = en/header.sample_rate;
                                    spike_flag = 0;
                                    % check if there is an artifact
                                    if ~any(curr_artifact.idx(st:en))
                                        % record if there is a spike, if so
                                        % move idx up
                                        if any((curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms))
                                            spike_flag = 1;
                                            % if there are spikes in this
                                            % window, move the start to the
                                            % first spike in the window
                                            st = min(curr_spike.pos((curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms)))*header.sample_rate - 1;
                                            en = st + (win_length*header.sample_rate - 1);
                                        end
                                        % check that we havent gone
                                        % past the end of the data
                                        if en <= dur
                                            % update spike idx
                                            spike_idx(cnt) = spike_flag;
                                            spike_num(cnt) = sum((curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms));
                                            spike_chan(cnt) = {curr_spike.chan((curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms))};
                                            % update trl
                                            trl(cnt,:) = [st, en, 0] + (trl_offset - 1);
                                            % add time
                                            time_vec(cnt) = (st + trl_offset - 1)/header.sample_rate;
                                            % update cnt
                                            cnt = cnt + 1;
                                        end
                                    end
                                    idx = en + 1;
                                end
                            end
                            % redefine trial
                            cfg = [];
                            cfg.trl = round(trl);
                            win_data = ft_redefinetrial(cfg,ft_data);
                            nTrial = numel(win_data.trial);
                            clear ft_data out_clean artifact_all
                            
                            % prewhiten
                            cfg = [];
                            cfg.derivative = 'yes';
                            ft_preprocessing(cfg, win_data);
                            
                            % band limited
                            fprintf('\nStarting multitaper FFT...\n')
                            % psd and csd - averging done over tapers
                            cfg = [];
                            cfg.method     = 'mtmfft';
                            cfg.taper      = 'dpss';
                            cfg.output     = 'powandcsd';
                            cfg.foi        = freqs;
                            cfg.tapsmofrq  = 4;
                            cfg.pad        = 'maxperlen';
                            cfg.keeptrials  = 'yes';
                            wave = ft_freqanalysis(cfg, win_data);
                            
                            % get labels for later
                            label = wave.label;
                            labelcmb = wave.labelcmb;
                            
                            % power
                            pow = zeros(nBand, nTrial, nElec);
                            fprintf('\npower...\n')
                            for i = 1:nBand
                                idx = wave.freq >= bands(i,1) & wave.freq <= bands(i,2);
                                pow(i,:,:) = mean(log10(wave.powspctrm(:,:,idx)),3);
                            end
                            
                            % coh
                            fprintf('\ncoherence...\n')
                            C = get_coh(wave,bands);
                            
                            % band limited, time resolved
                            fprintf('\nStarting Hilbert transform\n')
                            aec = zeros(nBand, nPair, nTrial);
                            plv = nan(nBand, nPair, nTrial);
                            bp_all_bands = cell(nBand,1);
                            for i = 1:nBand
                                fprintf('\n%s band...\n', band_names{i})
                                
                                curr_range = bands(i,:);
                                cfg = [];
                                cfg.bpfilter = 'yes';
                                cfg.bpfreq = curr_range;
                                cfg.bpfiltdf = 0.5; % bandpass transition width
                                cfg.bsfiltdf = 0.5; % bandstop transition width
                                cfg.bpfiltdev = 0.01; % bandpass max passband deviation
                                cfg.bsfiltdev = 0.05; % bandstp max passband deviation
                                cfg.bpfilttype = 'firws'; % or 'firls' (slower), but avoid the default 'but' (= not well-suited for hilbert phase estimate)
                                cfg.hilbert = 'complex';
                                
                                bp_data = ft_preprocessing(cfg, win_data);
                                bp_all_bands{i} = bp_data;
                                
                                % amp corr
                                fprintf('\namplitude envelope correlation...\n')
                                % do I need to orthogonalize?
                                for j = 1:nTrial
                                    full_corr = corr(abs(bp_data.trial{j}'));
                                    aec(i,:,j) = full_corr(upper_tri);
                                end
                                
                                % plv
                                fprintf('\nphaselocking value...\n')
                                % note that this not interpretable for wide
                                % band signals, like high gamma
                                if curr_range(1) < 60
                                    for j = 1:nTrial
                                        curr_phase = atan2(imag(bp_data.trial{j}),real(bp_data.trial{j}));
                                        cnt = 0;
                                        for k = 1:nElec
                                            for m = (k+1):nElec
                                                cnt = cnt + 1;
                                                dphase = curr_phase(k,:) - curr_phase(m,:);
                                                plv(i,cnt,j) = abs(mean(exp(1i*dphase)));
                                            end
                                        end
                                    end
                                end
                            end
                            
                            % cross freq coupling
                            % pac
                            fprintf('\nStarting PAC...\n')
                            pac = zeros(nBandCF, nPair, nTrial);
                            mi = zeros(nBandCF, nPair, nTrial);
%                             for i = 1:nBandCF
%                                 fprintf('\n%s band...\n', band_names{bands_cf(i,1)})
%                                 curr_phase = bp_all_bands{bands_cf(i,1)};
%                                 curr_amp = bp_all_bands{bands_cf(i,2)};
%                                 [pac(i,:,:), mi(i,:,:)] = get_pac(curr_phase, curr_amp);
%                             end
                            
                            % broadband
                            fprintf('\nStarting low-pass filter...\n')
                            % 4th order butterworth filter for 200 hz low pass
                            cfg = [];
                            cfg.lpfilt = 'yes';
                            cfg.lpfilttype = 'but';
                            cfg.lpfiltord = 4;
                            lfp = ft_preprocessing(cfg, win_data);
                            
                            % FC
                            % xcor
                            fprintf('\ncross-correlation...\n')
                            xcorr_lfp = zeros(nTrial, nPair);
                            for i = 1:nTrial
                                full_xcorr = max(xcorr(lfp.trial{i}','normalized'));
                                % get upper triangle
                                xcorr_lfp(i,:) = full_xcorr(upper_tri);
                            end
                            
                            % ar
                            fprintf('\nauto-regressive model...\n')
                            ar = zeros(nTrial,nElec^2);
                            labelcmb_dir = cell(nElec^2,2);
                            cnt = 0;
                            for i = 1:nElec:nElec^2
                                cnt = cnt + 1;
                                labelcmb_dir(i:(i+nElec-1),1) = label(cnt);
                                labelcmb_dir(i:(i+nElec-1),2) = label;
                            end
                            for i = 1:nTrial
                                [~, Aest] = arfit(lfp.trial{i}', pmin, pmax);
                                ar(i,:) = reshape(Aest, [], 1);
                            end
                            
                            fprintf('Done!\n'); toc
                            % save things
                            cell_ft = struct2cell(win_data);
                            fields = fieldnames(win_data);
                            keep_idx = ~(strcmpi(fields, 'time') | strcmpi(fields, 'trial'));
                            ft_header = cell2struct(cell_ft(keep_idx), fields(keep_idx));
                            
                            % add fields for mtmFFT
                            fc_header = ft_header;
                            fc_header.analysis_cfg = wave.cfg;
                            
                            % power
                            save([save_dir, 'power.mat'], 'pow', 'bands', 'label', 'fc_header', 'spike_idx')
                            
                            % coh
                            save([save_dir, 'mtcoherence.mat'], 'C', 'bands', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % add fields for hilbert
                            fc_header = ft_header;
                            fc_header.analysis_cfg = bp_data.cfg;
                            
                            % aec
                            save([save_dir, 'amp_env_corr.mat'], 'aec', 'bands', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % plv
                            save([save_dir, 'phase_lock_val.mat'], 'plv', 'bands', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % pac
                            save([save_dir, 'phase_amp_coupl.mat'], 'pac', 'mi', 'bands', 'bands_cf', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % add fields for low pass
                            fc_header = ft_header;
                            fc_header.analysis_cfg = lfp.cfg;
                            
                            % xcor
                            save([save_dir, 'crosscorr.mat'], 'xcorr_lfp', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % ar
                            save([save_dir, 'autoreg.mat'], 'ar', 'labelcmb_dir', 'fc_header', 'spike_idx')
                            
                            
                            % get strengths, and add to table
                            soz_idx = cellfun(@(x) any(strcmp(x, soz)), labelcmb(:,1)) |...
                                cellfun(@(x) any(strcmp(x, soz)), labelcmb(:,2));
                            spike_idx = cellfun(@(x) any(strcmp(x, interictal_cont)), labelcmb(:,1)) |...
                                cellfun(@(x) any(strcmp(x, interictal_cont)), labelcmb(:,2));

                            for i = 1:nMeasures
                                curr_measure = measure_names{i};
                                switch curr_measure
                                    case [{'coh'}, {'plv'}, {'aec'}]
                                        % get measure of interest
                                        if strcmp(curr_measure, 'coh')
                                            curr_data = C;
                                        elseif strcmp(curr_measure, 'plv')
                                            curr_data = plv;
                                        elseif strcmp(curr_measure, 'aec')
                                            curr_data = aec;
                                        end
                                        
                                        % will becoms columns in dataframe
                                        soz_str = nan(nElec*nTrial*nBand,1);
                                        not_soz_str = nan(nElec*nTrial*nBand,1);
                                        spike_str = nan(nElec*nTrial*nBand,1);
                                        not_spike_str = nan(nElec*nTrial*nBand,1);
                                        elec_order = cell(nElec*nTrial*nBand,1);
                                        band_order = cell(nElec*nTrial*nBand,1);
                                        elec_in_soz = nan(nElec*nTrial*nBand,1);
                                        elec_in_spike = nan(nElec*nTrial*nBand,1);
                                        region_order = cell(nElec*nTrial*nBand,1);
                                        spike_nums = nan(nElec*nTrial*nBand,1);
                                        time = nan(nElec*nTrial*nBand,1);
                                        control_pow = nan(nElec*nTrial*nBand,1);
                                        
                                        % add strengths for all bands
                                        % add strengths
                                        cnt = 1;
                                        for k = 1:nBand
                                            for j = 1:nElec
                                                % get indices and other stuff
                                                curr = label{j};
                                                region = regions{j};
                                                chan_idx = find(strcmp(label,label{j}));
                                                spike_flag = cellfun(@(x) any(chan_idx == x), spike_chan);
                                                elec_idx = cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,1)) |...
                                                    cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,2));
                                                % get strengths
                                                soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & soz_idx,:),2);
                                                not_soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & ~soz_idx,:),2);
                                                spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & spike_idx,:),2);
                                                not_spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & ~spike_idx,:),2);
                                                
                                                % get other elec vars
                                                elec_order(cnt:(cnt+nTrial-1)) = repmat({curr}, nTrial, 1);
                                                elec_in_soz(cnt:(cnt+nTrial-1)) = repmat(any(strcmp(curr, soz)), nTrial, 1);
                                                elec_in_spike(cnt:(cnt+nTrial-1)) = spike_flag;
                                                region_order(cnt:(cnt+nTrial-1)) = {region};
                                                
                                                % get other spike vars
                                                spike_nums(cnt:(cnt+nTrial-1)) = spike_num;
                                                
                                                %time vars
                                                time(cnt:(cnt+nTrial-1)) = time_vec;
                                                
                                                % power in relevant band
                                                control_pow(cnt:(cnt+nTrial-1)) = pow(k,:,chan_idx);
                                                
                                                %band name
                                                band_order(cnt:(cnt+nTrial-1)) = {band_names(k)};
                                                
                                                % update counter
                                                cnt = cnt+nTrial;
                                            end

                                        end
                                        
                                        % add things that are consistent
                                        meas_order = repmat({curr_measure}, nTrial*nElec*nBand, 1);
                                        subj_order = repmat({subj}, nTrial*nElec*nBand, 1);
                                        exper_order = repmat({exper}, nTrial*nElec*nBand, 1);
                                        sess_order = repmat({sess}, nTrial*nElec*nBand, 1);
                                        age_order = repmat({age}, nTrial*nElec*nBand, 1);
                                        gender_order = repmat({gender}, nTrial*nElec*nBand, 1);
                                        hand_order = repmat({hand}, nTrial*nElec*nBand, 1);
                                        race_order = repmat({race}, nTrial*nElec*nBand, 1);
                                        
                                        % add to table
                                        fc_table = [fc_table; table(subj_order, exper_order, sess_order, time,...
                                            control_pow, meas_order, band_order, soz_str, not_soz_str, spike_str,...
                                            not_spike_str, elec_order, region_order, elec_in_soz, elec_in_spike,...
                                            spike_nums, age_order, gender_order, race_order, hand_order, 'VariableNames', table_names)];
                                    case [{'xcorr'}, {'ar'}]
                                        % will becoms columns in dataframe
                                        soz_str = nan(nElec*nTrial,1);
                                        not_soz_str = nan(nElec*nTrial,1);
                                        spike_str = nan(nElec*nTrial,1);
                                        not_spike_str = nan(nElec*nTrial,1);
                                        elec_order = cell(nElec*nTrial,1);
                                        elec_in_soz = nan(nElec*nTrial,1);
                                        elec_in_spike = nan(nElec*nTrial,1);
                                        region_order = cell(nElec*nTrial,1);
                                        spike_nums = nan(nElec*nTrial,1);
                                        time = nan(nElec*nTrial,1);
                                        control_pow = nan(nElec*nTrial,1);
                                        
                                        if strcmp(curr_measure, 'xcorr')
                                            curr_data = xcorr_lfp;
                                        elseif strcmp(curr_measure, 'ar')
                                            curr_data = ar;
                                        end
                                        
                                        % add strengths
                                        cnt = 1;
                                        for j = 1:nElec
                                            % get indices and other stuff
                                            curr = label{j};
                                            region = regions{j};
                                            chan_idx = find(strcmp(label,label{j}));
                                            spike_flag = cellfun(@(x) any(chan_idx == x), spike_chan);
                                            elec_idx = cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,1)) |...
                                                cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,2));
                                            % get strengths
                                            soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & soz_idx),2);
                                            not_soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & ~soz_idx),2);
                                            spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & spike_idx),2);
                                            not_spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & ~spike_idx),2);
                                            
                                            % get other elec vars
                                            elec_order(cnt:(cnt+nTrial-1)) = repmat({curr}, nTrial, 1);
                                            elec_in_soz(cnt:(cnt+nTrial-1)) = repmat(any(strcmp(curr, soz)), nTrial, 1);
                                            elec_in_spike(cnt:(cnt+nTrial-1)) = spike_flag;
                                            region_order(cnt:(cnt+nTrial-1)) = {region};
                                            
                                            % get other spike vars
                                            spike_nums(cnt:(cnt+nTrial-1)) = spike_num;
                                            
                                            %time vars
                                            time(cnt:(cnt+nTrial-1)) = time_vec;
                                            
                                            % power in HG band
                                            control_pow(cnt:(cnt+nTrial-1)) = pow(end,:,chan_idx);
                                            
                                            % update counter
                                            cnt = cnt+nTrial;
                                        end
                                        
                                        % add things that are constant
                                        % across elecs
                                        band_order = repmat({'broadband'}, nTrial*nElec, 1);
                                        meas_order = repmat({curr_measure}, nTrial*nElec, 1);
                                        subj_order = repmat({subj}, nTrial*nElec, 1);
                                        exper_order = repmat({exper}, nTrial*nElec, 1);
                                        sess_order = repmat({sess}, nTrial*nElec, 1);
                                        age_order = repmat({age}, nTrial*nElec, 1);
                                        gender_order = repmat({gender}, nTrial*nElec, 1);
                                        race_order = repmat({race}, nTrial*nElec, 1);
                                        hand_order = repmat({hand}, nTrial*nElec, 1);
                                        
                                        % add to table
                                        fc_table = [fc_table; table(subj_order, exper_order, sess_order, time,...
                                            control_pow, meas_order, band_order, soz_str, not_soz_str, spike_str,...
                                            not_spike_str, elec_order, region_order, elec_in_soz, elec_in_spike,...
                                            spike_nums, age_order, gender_order, race_order, hand_order, 'VariableNames', table_names)];
                                    case 'pac'
                                        curr_data = pac;
                                end
                            end
                        end
                    end
                end
            end
            % save subject table
            writetable(fc_table, [save_dir, 'fc_data.csv'])
        end
    end
end
function [] = functional_connectivity_cfn(protocol, release, top_dir, subj, detector, spike_win, win_length)
% main function for functional connectivity - helps with paralelizing
% do you want to save all your FC mtrices? faster if no
save_flag = false;
release_dir = [top_dir, 'release', release '/'];

% get global info struct
fname = [release_dir 'protocols/', protocol, '.json'];
fid = fopen(fname);
raw = fread(fid);
str = char(raw');
fclose(fid);
info = jsondecode(str);
eval(['info = info.protocols.', protocol,';']);

% useful variables
table_names = [{'subj'}, {'exper'}, {'sess'}, {'time'}, {'power'}, {'fc_measure'}, {'band'}, {'str'}, {'ti'},...
    {'str_soz'}, {'str_not_soz'}, {'str_spike'}, {'str_not_spike'},{'str_grid'},{'str_depth'}, {'wm_depth'}, {'gm_depth'},...
    {'elec'}, {'region'}, {'wm'}, {'elec_in_soz'},...
    {'elec_has_spike'}, {'spike_num'}, {'spike_spread'}, {'age'}, {'gender'}, {'race'}, {'hand'}, {'x'}, {'y'}, {'z'}, {'type'}];
%bands
freqs = unique(round(logspace(log10(4),log10(150),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70; 71, 150];
band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}, {'hg'}];

%fc measures
measure_names = [{'coh'}, {'im_coh'}, {'plv'}, {'iplv'}, {'aec'}, {'aec_ortho'}, {'xcorr'}, {'ar'}, {'pac'}];
%parameters
pmin = 1; pmax = 1; % order for AR model
% constants
nBand = size(bands,1);
nMeasures = numel(measure_names);
% subjects not to use
load([top_dir, 'bad_datasets.mat'])
errors = struct('files', [], 'message', []);

% make subject directory
subj_dir = [top_dir, 'FC/release',release, '/', protocol, '/', subj, '/'];
if ~exist(subj_dir, 'dir')
    mkdir(subj_dir);
end

% initialize table
power_vars = cellfun(@(x) ['power_', x], band_names, 'UniformOutput', false);
coh_vars = cellfun(@(x) ['coh_', x], band_names, 'UniformOutput', false);
plv_vars = cellfun(@(x) ['plv_', x], band_names, 'UniformOutput', false);
aec_vars = cellfun(@(x) ['aec_', x], band_names, 'UniformOutput', false);

fc_table = cell2table(cell(0,32), 'VariableNames', table_names);

if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
    mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
end

% check that we need data for this subj
if ~exist([top_dir, 'FC/release',release, '/', protocol, '/', subj, '/', 'win_', num2str(win_length), '/fc_data', detector, '.csv'], 'file')
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
            if ~exist([subj_dir, 'win_', num2str(win_length), '/'], 'dir')
                mkdir([subj_dir, 'win_', num2str(win_length), '/']);
            end
            
            if exist([data_dir, 'data_clean.mat'], 'file') && exist([data_dir, 'spike_info_', num2str(spike_win), '.mat'], 'file')
                load([data_dir, 'data_clean.mat'])
                load([data_dir, 'header.mat'])
                load([data_dir, 'channel_info.mat'])
                if strcmp(detector, '_delphos')
                    load([data_dir, 'spike_info', detector, '.mat'])
                else
                    load([data_dir, 'spike_info_', num2str(spike_win), '.mat'])
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
                    nPair = (nElec^2-nElec)/2;
                    
                    % constants
                    lower_tri = reshape(tril(true(nElec),-1),[],1);
                    
                    % get window start times
                    trl = [];
                    spike_idx = [];
                    spike_spread = [];
                    spike_num = [];
                    spike_chan = {};
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
                                    % get all the spikes in the window
                                    curr_idx = (curr_spike.pos >= st_ms) & (curr_spike.pos <= en_ms);
                                    seqs = unique(curr_spike.seq(curr_idx));
                                    % update spike idx
                                    spike_idx(cnt) = spike_flag;
                                    spike_num(cnt) = numel(unique(curr_spike.seq(curr_idx)));
                                    if spike_flag
                                        spread = zeros(spike_num(cnt), 1);
                                        for m = 1:numel(spread)
                                            spread(m) = numel(unique(curr_spike.chan(curr_spike.seq == seqs(m) & curr_idx)));
                                        end
                                        spike_spread(cnt) = mean(spread);
                                    else
                                        spike_spread(cnt) = 0;
                                    end
                                    spike_chan(cnt) = {curr_spike.chan(curr_idx)};
                                    % update trl
                                    trl(cnt,:) = [st + (trl_offset - 1), en + (trl_offset - 1), 0];
                                    % add time
                                    time_vec(cnt) = (st + trl_offset - 1)/header.sample_rate;
                                    % update cnt
                                    cnt = cnt + 1;
                                end
                            end
                            idx = en + 1;
                        end
                    end
                    
                    % check that we found at least one window with
                    % a spike
                    if sum(spike_num) > 0
                        
                        % redefine trial
                        cfg = [];
                        cfg.trl = round(trl);
                        win_data = ft_redefinetrial(cfg,ft_data);
                        nTrial = numel(win_data.trial);
                        clear ft_data artifact_all
                        
                        % check if something went wrong in the trial
                        % definition, and it was getting the same data
                        % twice
                        
                        
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
                        C = get_coh(wave,bands, 'regular');
                        Ci = get_coh(wave,bands,'imaginary');
                        
                        % band limited, time resolved
                        fprintf('\nStarting Hilbert transform\n')
                        aec = zeros(nBand, nPair, nTrial);
                        aec_ortho = zeros(nBand, nPair, nTrial);
                        plv = nan(nBand, nPair, nTrial);
                        iplv = nan(nBand, nPair, nTrial);
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
                            % unorthogonalized
                            for j = 1:nTrial
                                full_corr = corr(abs(bp_data.trial{j}'));
                                aec(i,:,j) = full_corr(lower_tri);
                            end
                            
                            % orthogonalized Brookes et al., 2012, 2014
                            % if z1 and z2 are normalized so that
                            % mean(abs(zi)^2) = 1, then we replace
                            % z2 with z2 - R(c)*z1, where R(c) is
                            % the real part of coherence
                            for j = 1:nTrial
                                aec_ortho(i,:,j) = get_aec_ortho_cfn(bp_data.trial{j});
                            end
                            
                            % plv
                            fprintf('\nphaselocking value...\n')
                            % note that this not interpretable for wide
                            % band signals, like high gamma
                            if curr_range(1) < 60
                                for j = 1:nTrial
                                    curr_phase = bp_data.trial{j}./abs(bp_data.trial{j});
                                    curr_plv = abs(curr_phase*curr_phase')/(size(bp_data.trial{j},2));
                                    curr_iplv = imag(curr_phase*curr_phase')/size(bp_data.trial{j},2);
                                    plv(i,:,j) = curr_plv(lower_tri);
                                    iplv(i,:,j) = curr_iplv(lower_tri);
                                end
                            end
                        end
                        
                        % cross freq coupling
                        % pac
                        %fprintf('\nStarting PAC...\n')
                        %pac = zeros(nBandCF, nPair, nTrial);
                        %mi = zeros(nBandCF, nPair, nTrial);
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
                        cfg.lpfreq = 200;
                        cfg.lpfiltord = 4;
                        
                        lfp = ft_preprocessing(cfg, win_data);
                        
                        % FC
                        % xcor
                        fprintf('\ncross-correlation...\n')
                        xcorr_lfp = zeros(nTrial, nPair);
                        for i = 1:nTrial
                            full_xcorr = max(xcorr(lfp.trial{i}')./reshape(sqrt(sum(lfp.trial{i}.^2)*sum(lfp.trial{i}.^2)'),1,[]));
                            % get upper triangle
                            xcorr_lfp(i,:) = full_xcorr(lower_tri);
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
                        
                        fprintf('Done!\n');
                        if save_flag
                            % save things
                            cell_ft = struct2cell(win_data);
                            names_ft = fieldnames(win_data);
                            keep_idx = ~(strcmpi(names_ft, 'time') | strcmpi(names_ft, 'trial'));
                            ft_header = cell2struct(cell_ft(keep_idx), names_ft(keep_idx));
                            
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
                            save([save_dir, 'amp_env_corr_ortho.mat'], 'aec_ortho', 'bands', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % plv
                            save([save_dir, 'phase_lock_val.mat'], 'plv', 'bands', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % pac
                            %save([save_dir, 'phase_amp_coupl.mat'], 'pac', 'mi', 'bands', 'bands_cf', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % add fields for low pass
                            fc_header = ft_header;
                            fc_header.analysis_cfg = lfp.cfg;
                            
                            % xcor
                            save([save_dir, 'crosscorr.mat'], 'xcorr_lfp', 'labelcmb', 'fc_header', 'spike_idx')
                            
                            % ar
                            save([save_dir, 'autoreg.mat'], 'ar', 'labelcmb_dir', 'fc_header', 'spike_idx')
                        end
                        
                        
                        % get strengths, and add to table
                        soz_idx = cellfun(@(x) any(strcmp(x, soz)), labelcmb(:,1)) |...
                            cellfun(@(x) any(strcmp(x, soz)), labelcmb(:,2));
                        soz_idx_dir = cellfun(@(x) any(strcmp(x, soz)), labelcmb_dir(:,1)) |...
                            cellfun(@(x) any(strcmp(x, soz)), labelcmb_dir(:,2));
                        spike_idx = cellfun(@(x) any(strcmp(x, interictal_cont)), labelcmb(:,1)) |...
                            cellfun(@(x) any(strcmp(x, interictal_cont)), labelcmb(:,2));
                        spike_idx_dir = cellfun(@(x) any(strcmp(x, interictal_cont)), labelcmb_dir(:,1)) |...
                            cellfun(@(x) any(strcmp(x, interictal_cont)), labelcmb_dir(:,2));
                        % get grid and depth labels
                        grid = label(strcmp(elec_type, 'G'));
                        depth = label(strcmp(elec_type, 'D'));
                        wm_label = label(wm);
                        grid_idx = cellfun(@(x) any(strcmp(x, grid)), labelcmb(:,1)) |...
                            cellfun(@(x) any(strcmp(x, grid)), labelcmb(:,2));
                        depth_idx = cellfun(@(x) any(strcmp(x, depth)), labelcmb(:,1)) |...
                            cellfun(@(x) any(strcmp(x, depth)), labelcmb(:,2));
                        depth_wm_idx = cellfun(@(x) any(strcmp(x, depth)), labelcmb(:,1)) |...
                            cellfun(@(x) any(strcmp(x, depth)), labelcmb(:,2));
                        grid_idx_dir = cellfun(@(x) any(strcmp(x, grid)), labelcmb_dir(:,1)) |...
                            cellfun(@(x) any(strcmp(x, grid)), labelcmb_dir(:,2));
                        depth_idx_dir = cellfun(@(x) any(strcmp(x, depth)), labelcmb_dir(:,1)) |...
                            cellfun(@(x) any(strcmp(x, depth)), labelcmb_dir(:,2));
                        depth_wm_idx_dir = cellfun(@(x) any(strcmp(x, wm_label)), labelcmb_dir(:,1)) |...
                            cellfun(@(x) any(strcmp(x, wm_label)), labelcmb_dir(:,2));
                        
                        for i = 1:nMeasures
                            curr_measure = measure_names{i};
                            switch curr_measure
                                case [{'coh'}, {'im_coh'}, {'plv'}, {'iplv'}, {'aec'}, {'aec_ortho'}]
                                    % get measure of interest
                                    if strcmp(curr_measure, 'coh')
                                        curr_data = C;
                                    elseif strcmp(curr_measure, 'im_coh')
                                        curr_data = Ci;
                                    elseif strcmp(curr_measure, 'plv')
                                        curr_data = plv;
                                    elseif strcmp(curr_measure, 'iplv')
                                        curr_data = abs(iplv);
                                    elseif strcmp(curr_measure, 'aec')
                                        curr_data = abs(aec);
                                    elseif strcmp(curr_measure, 'aec_ortho')
                                        curr_data = abs(aec_ortho);
                                    end
                                    
                                    % will becoms columns in dataframe
                                    str = nan(nElec*nTrial*nBand,1);
                                    ti = nan(nElec*nTrial*nBand,1);
                                    soz_str = nan(nElec*nTrial*nBand,1);
                                    not_soz_str = nan(nElec*nTrial*nBand,1);
                                    spike_str = nan(nElec*nTrial*nBand,1);
                                    grid_str = nan(nElec*nTrial*nBand,1);
                                    depth_str = nan(nElec*nTrial*nBand,1);
                                    wm_str = nan(nElec*nTrial*nBand,1);
                                    gm_str = nan(nElec*nTrial*nBand,1);
                                    not_spike_str = nan(nElec*nTrial*nBand,1);
                                    elec_order = cell(nElec*nTrial*nBand,1);
                                    band_order = cell(nElec*nTrial*nBand,1);
                                    elec_in_soz = nan(nElec*nTrial*nBand,1);
                                    elec_in_spike = nan(nElec*nTrial*nBand,1);
                                    region_order = cell(nElec*nTrial*nBand,1);
                                    x = zeros(nElec*nTrial*nBand,1);
                                    y = zeros(nElec*nTrial*nBand,1);
                                    z = zeros(nElec*nTrial*nBand,1);
                                    type = cell(nElec*nTrial*nBand,1);
                                    wm_order = cell(nElec*nTrial*nBand,1);
                                    spike_nums = nan(nElec*nTrial*nBand,1);
                                    spike_spreads = nan(nElec*nTrial*nBand,1);
                                    time = nan(nElec*nTrial*nBand,1);
                                    control_pow = nan(nElec*nTrial*nBand,1);
                                    
                                    % add strengths for all bands
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
                                            str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx, :),2);
                                            ti(cnt:(cnt+nTrial-1)) = skewness(curr_data(k, elec_idx, :));
                                            if any(soz_idx)
                                                soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & soz_idx,:),2);
                                                not_soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & ~soz_idx,:),2);
                                            end
                                            if any(spike_idx)
                                                spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & spike_idx,:),2);
                                                not_spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & ~spike_idx,:),2);
                                            end
                                            if any(depth_wm_idx)
                                                wm_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & depth_wm_idx & depth_idx,:),2);
                                                gm_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & ~depth_wm_idx & depth_idx,:),2);
                                            end
                                            if any(grid_idx)
                                                grid_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & grid_idx,:),2);
                                            end
                                            if any(depth_idx)
                                                depth_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & depth_idx,:),2);
                                            end
                                            
                                            % get other elec vars
                                            elec_order(cnt:(cnt+nTrial-1)) = repmat({curr}, nTrial, 1);
                                            elec_in_soz(cnt:(cnt+nTrial-1)) = repmat(any(strcmp(curr, soz)), nTrial, 1);
                                            elec_in_spike(cnt:(cnt+nTrial-1)) = spike_flag;
                                            region_order(cnt:(cnt+nTrial-1)) = {region};
                                            x(cnt:(cnt+nTrial-1)) = mni_coords{j}(1);
                                            y(cnt:(cnt+nTrial-1)) = mni_coords{j}(2);
                                            z(cnt:(cnt+nTrial-1)) = mni_coords{j}(3);
                                            type(cnt:(cnt+nTrial-1)) = elec_type(j);
                                            wm_order(cnt:(cnt+nTrial-1)) = {wm(j)};
                                            
                                            % get other spike vars
                                            spike_nums(cnt:(cnt+nTrial-1)) = spike_num;
                                            spike_spreads(cnt:(cnt+nTrial-1)) = spike_spread;
                                            
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
                                    
                                case [{'xcorr'}, {'ar'}]
                                    % will becoms columns in dataframe
                                    str = nan(nElec*nTrial,1);
                                    ti = nan(nElec*nTrial,1);
                                    soz_str = nan(nElec*nTrial,1);
                                    grid_str = nan(nElec*nTrial,1);
                                    wm_str = nan(nElec*nTrial,1);
                                    gm_str = nan(nElec*nTrial,1);
                                    depth_str = nan(nElec*nTrial,1);
                                    not_soz_str = nan(nElec*nTrial,1);
                                    spike_str = nan(nElec*nTrial,1);
                                    not_spike_str = nan(nElec*nTrial,1);
                                    elec_order = cell(nElec*nTrial,1);
                                    elec_in_soz = nan(nElec*nTrial,1);
                                    elec_in_spike = nan(nElec*nTrial,1);
                                    region_order = cell(nElec*nTrial,1);
                                    x = zeros(nElec*nTrial,1);
                                    y = zeros(nElec*nTrial,1);
                                    z = zeros(nElec*nTrial,1);
                                    type = cell(nElec*nTrial,1);
                                    wm_order = cell(nElec*nTrial,1);
                                    spike_nums = nan(nElec*nTrial,1);
                                    spike_spreads = nan(nElec*nTrial,1);
                                    time = nan(nElec*nTrial,1);
                                    control_pow = nan(nElec*nTrial,1);
                                    
                                    if strcmp(curr_measure, 'xcorr')
                                        curr_data = xcorr_lfp;
                                        elec_idx = cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,1)) |...
                                            cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,2));
                                        
                                        % add strengths
                                        cnt = 1;
                                        for j = 1:nElec
                                            % get indices and other stuff
                                            curr = label{j};
                                            region = regions{j};
                                            chan_idx = find(strcmp(label,label{j}));
                                            spike_flag = cellfun(@(x) any(chan_idx == x), spike_chan);
                                            
                                            % get strengths
                                            str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx),2);
                                            ti(cnt:(cnt+nTrial-1)) = skewness(curr_data(:, elec_idx),1,2);
                                            if any(soz_idx)
                                                soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & soz_idx),2);
                                                not_soz_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & ~soz_idx),2);
                                            end
                                            if any(spike_idx)
                                                spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & spike_idx),2);
                                                not_spike_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & ~spike_idx),2);
                                            end
                                            if any(depth_wm_idx)
                                                wm_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & depth_wm_idx & depth_idx,:),2);
                                                gm_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & ~depth_wm_idx & depth_idx,:),2);
                                            end
                                            if any(grid_idx)
                                                grid_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & grid_idx,:),2);
                                            end
                                            if any(depth_idx)
                                                depth_str(cnt:(cnt+nTrial-1)) = mean(curr_data(:, elec_idx & depth_idx,:),2);
                                            end
                                            
                                            % get other elec vars
                                            elec_order(cnt:(cnt+nTrial-1)) = repmat({curr}, nTrial, 1);
                                            elec_in_soz(cnt:(cnt+nTrial-1)) = repmat(any(strcmp(curr, soz)), nTrial, 1);
                                            elec_in_spike(cnt:(cnt+nTrial-1)) = spike_flag;
                                            region_order(cnt:(cnt+nTrial-1)) = {region};
                                            x(cnt:(cnt+nTrial-1)) = mni_coords{j}(1);
                                            y(cnt:(cnt+nTrial-1)) = mni_coords{j}(2);
                                            z(cnt:(cnt+nTrial-1)) = mni_coords{j}(3);
                                            type(cnt:(cnt+nTrial-1)) = elec_type(j);
                                            wm_order(cnt:(cnt+nTrial-1)) = {wm(j)};
                                            
                                            % get other spike vars
                                            spike_nums(cnt:(cnt+nTrial-1)) = spike_num;
                                            spike_spreads(cnt:(cnt+nTrial-1)) = spike_spread;
                                            
                                            %time vars
                                            time(cnt:(cnt+nTrial-1)) = time_vec;
                                            
                                            % power in HG band
                                            control_pow(cnt:(cnt+nTrial-1)) = pow(end,:,chan_idx);
                                            
                                            % update counter
                                            cnt = cnt+nTrial;
                                        end
                                        
                                    elseif strcmp(curr_measure, 'ar')
                                        curr_data = ar;
                                        % get index for
                                        % directed connectivity
                                        elec_idx = cellfun(@(x) any(strcmp(x, curr)), labelcmb_dir(:,1)) |...
                                            cellfun(@(x) any(strcmp(x, curr)), labelcmb_dir(:,2));
                                        
                                        % add strengths
                                        cnt = 1;
                                        for j = 1:nElec
                                            % get indices and other stuff
                                            curr = label{j};
                                            region = regions{j};
                                            chan_idx = find(strcmp(label,label{j}));
                                            spike_flag = cellfun(@(x) any(chan_idx == x), spike_chan);
                                            
                                            % get strengths
                                            str(cnt:(cnt+nTrial-1)) = mean(abs(curr_data(:, elec_idx)),2);
                                            ti(cnt:(cnt+nTrial-1)) = skewness(abs(curr_data(:, elec_idx)),1,2);
                                            if any(soz_idx_dir)
                                                soz_str(cnt:(cnt+nTrial-1)) = mean(abs(curr_data(:, elec_idx & soz_idx_dir)),2);
                                                not_soz_str(cnt:(cnt+nTrial-1)) = mean(abs(curr_data(:, elec_idx & ~soz_idx_dir)),2);
                                            end
                                            if any(spike_idx_dir)
                                                spike_str(cnt:(cnt+nTrial-1)) = mean(abs(curr_data(:, elec_idx & spike_idx_dir)),2);
                                                not_spike_str(cnt:(cnt+nTrial-1)) = mean(abs(curr_data(:, elec_idx & ~spike_idx_dir)),2);
                                            end
                                            if any(depth_wm_idx_dir)
                                                wm_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & depth_wm_idx_dir & depth_idx_dir,:),2);
                                                gm_str(cnt:(cnt+nTrial-1)) = mean(curr_data(k, elec_idx & ~depth_wm_idx_dir & depth_idx_dir,:),2);
                                            end
                                            if any(grid_idx_dir)
                                                grid_str(cnt:(cnt+nTrial-1)) = mean(abs(curr_data(:, elec_idx & grid_idx_dir)),2);
                                            end
                                            if any(depth_idx_dir)
                                                depth_str(cnt:(cnt+nTrial-1)) = mean(abs(curr_data(:, elec_idx & depth_idx_dir)),2);
                                            end
                                            
                                            % get other elec vars
                                            elec_order(cnt:(cnt+nTrial-1)) = repmat({curr}, nTrial, 1);
                                            elec_in_soz(cnt:(cnt+nTrial-1)) = repmat(any(strcmp(curr, soz)), nTrial, 1);
                                            elec_in_spike(cnt:(cnt+nTrial-1)) = spike_flag;
                                            region_order(cnt:(cnt+nTrial-1)) = {region};
                                            x(cnt:(cnt+nTrial-1)) = mni_coords{j}(1);
                                            y(cnt:(cnt+nTrial-1)) = mni_coords{j}(2);
                                            z(cnt:(cnt+nTrial-1)) = mni_coords{j}(3);
                                            type(cnt:(cnt+nTrial-1)) = elec_type(j);
                                            wm_order(cnt:(cnt+nTrial-1)) = {wm(j)};
                                            
                                            % get other spike vars
                                            spike_nums(cnt:(cnt+nTrial-1)) = spike_num;
                                            spike_spreads(cnt:(cnt+nTrial-1)) = spike_spread;
                                            
                                            %time vars
                                            time(cnt:(cnt+nTrial-1)) = time_vec;
                                            
                                            % power in HG band
                                            control_pow(cnt:(cnt+nTrial-1)) = pow(end,:,chan_idx);
                                            
                                            % update counter
                                            cnt = cnt+nTrial;
                                        end
                                        
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
                                    hand_order = repmat({hand}, nTrial*nElec, 1);
                                    race_order = repmat({race}, nTrial*nElec, 1);
                                case 'pac'
                                    %curr_data = pac;
                            end
                            
                            % add to table
                            fc_table = [fc_table; table(subj_order, exper_order, sess_order, time,...
                                control_pow, meas_order, band_order, str, ti, soz_str, not_soz_str, spike_str,...
                                not_spike_str, grid_str, depth_str, wm_str, gm_str, elec_order, region_order, wm_order, elec_in_soz, elec_in_spike,...
                                spike_nums, spike_spreads, age_order, gender_order, race_order, hand_order, x, y, z, type, 'VariableNames', table_names)];
                        end
                    end
                    
                end
                %                 catch ME
                %                     errors(end+1).files = [subj, '_', exper, '_', sess];
                %                     errors(end).message = ME.message;
                %                 end
            end
        end
    end
    % save subject table
    if size(fc_table,1) > 0
        writetable(fc_table, [subj_dir, 'win_', num2str(win_length), '/fc_data', detector, '.csv'])
    end
end

save([subj_dir, 'fc_errors_win', num2str(win_length), detector, '.mat'], 'errors');
end


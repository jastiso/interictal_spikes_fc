clear
clc
restoredefaultpath
addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/IRASA/'))


% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';
img_dir = '/Users/stiso/Documents/Code/interictal_spikes_fc/img/';
releases = ['1'];
win_length = 3000; % in ms
step = 1000; %in ms
filter = 0; % with ot without lowpass alaising filter
nSim = 50;
factors = linspace(1,2,50);
freqs = unique(round(logspace(log10(4),log10(70),30)));
bands = [4, 8; 9, 15; 16 25; 36, 70];
band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}];
nBand = numel(band_names);
nPair = ((nSim+1)^2-(nSim + 1))/2;
lower_tri = reshape(tril(true(nSim+1),-1),[],1);
nRandPhase = 100; %number of random phase values
%%
for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    % just for testing that qsub fuction will work
    for p = 1
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
        for s = 1
            subj = subjects{s};
            % main function for functional connectivity - helps with paralelizing
            
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
            
            fprintf('\n******************************************\nStarting sim for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1 % this should be the same across all sessions
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    data_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    
                    
                    if exist([data_dir, 'channel_info.mat'], 'file')
                        load([data_dir, 'data_clean.mat'])
                        load([data_dir, 'channel_info.mat'])
                        
                        % pick 2 channels, and show IRASA
                        c1 = 1;
                        c2 = 2;
                        
                        spec1 = get_IRASA_spec(ft_data.trial{1}(c1,:), 1, numel(ft_data.trial{1}(c1,:)), ft_data.fsample, win_length, step, filter);
                        spec2 = get_IRASA_spec(ft_data.trial{1}(c2,:), 1, numel(ft_data.trial{1}(c2,:)), ft_data.fsample, win_length, step, filter);
                        figure(1); clf;
                        subplot(2,2,1);
                        loglog(spec1.freq,mean(spec1.mixd,2),'b',  'linewidth', 2); hold on
                        loglog(spec1.freq,mean(spec1.frac,2),'r', 'linewidth', 2);
                        subplot(2,2,2);
                        plot(spec1.freq, mean(spec1.osci,2), 'linewidth', 2); hold on
                        shade_plot(spec1.freq', mean(spec1.osci,2)', (std(spec1.osci,[],2))', rgb("slategrey"), 0.4);
                        subplot(2,2,3);
                        loglog(spec2.freq,mean(spec2.mixd,2),'b',  'linewidth', 2); hold on
                        loglog(spec2.freq,mean(spec2.frac,2),'r', 'linewidth', 2);
                        subplot(2,2,4);
                        plot(spec2.freq, mean(spec2.osci,2), 'linewidth', 2); hold on
                        shade_plot(spec2.freq', mean(spec2.osci,2)', (std(spec2.osci,[],2))', rgb("slategrey"), 0.4);
                        
                        nFreq = numel(spec1.freq);
                        x = spec1.freq;
                        y = mean(spec2.frac,2);
                        y1 = mean(spec1.frac,2);
                        P = polyfit(log10(x),log10(y1),1);
                        psds = zeros(nSim,nFreq);
                        
                        % systematically change 1/f
                        for k = 1:nSim
                            m_new = P(1)*factors(k);
                            b_new = log10(y(nFreq/2-1)) - m_new*log10(x(nFreq/2-1));
                            
                            psds(k,:) = 10.^(b_new + m_new*log10(x) + log10(mean(spec2.osci,2)));
                            figure(2); clf
                            plot(log10(spec2.freq),log10(mean(spec2.mixd,2)),'b',  'linewidth', 2); hold on
                            plot(log10(spec2.freq),log10(mean(spec2.frac,2)),'r', 'linewidth', 2);
                            plot(log10(spec2.freq),log10(psds(k,:)),'k',  'linewidth', 2); hold on
                            plot(log10(x),b_new + m_new*log10(spec2.freq),'g', 'linewidth', 2);
                        end
                        
                        % get out of frequency space (not exact timeseries,
                        % but with same characteristics
                        new_ts = zeros(nSim, length(psds(k,:))*2*nRandPhase);
                        for k = 1:nSim
                            dat = [];
                            for m = 1:nRandPhase
                                Sxx = [psds(k,:), fliplr(psds(k,:))];
                                curr = TimeseriesFromPSD(Sxx', ft_data.fsample);
                                dat = [curr;dat];
                            end
                            new_ts(k,:) = dat;
                        end
                        
                        % regression with power and FC (does slope predict
                        % FC, accounting for power)
                        % get into ft format
                        
                        data = [repmat(ft_data.trial{1}(c1,1:length(Sxx)),1,nRandPhase);new_ts];
                        labels = {'real'};
                        for k = 1:nSim
                            labels{numel(labels) + 1} = num2str(P(1)*factors(k));
                        end
                        trl = [[1:length(Sxx):length(Sxx)*nRandPhase]',...
                            [1:length(Sxx):length(Sxx)*nRandPhase]'+length(Sxx)-1, zeros(nRandPhase,1)];
                        sim_dat = fieldtrip_format(data, ft_data.fsample, labels, trl);
                        
                        % prewhiten
                        cfg = [];
                        cfg.derivative = 'yes';
                        ft_preprocessing(cfg, sim_dat);
                        
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
                        wave = ft_freqanalysis(cfg, sim_dat);
                        
                        % get labels for later
                        label = wave.label;
                        labelcmb = wave.labelcmb;
                        
                        % power
                        pow = zeros(nRandPhase, nBand, nSim+1);
                        fprintf('\npower...\n')
                        for i = 1:nBand
                            idx = wave.freq >= bands(i,1) & wave.freq <= bands(i,2);
                            pow(:,i,:) = squeeze(mean(log10(wave.powspctrm(:,:,idx)),3));
                        end
                        
                        % coh
                        fprintf('\ncoherence...\n')
                        Ci = get_coh(wave,bands,'imaginary');
                        
                        % band limited, time resolved
                        fprintf('\nStarting Hilbert transform\n')
                        aec_ortho = zeros(nRandPhase, nBand, nPair);
                        iplv = nan(nRandPhase, nBand, nPair);
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
                            
                            bp_data = ft_preprocessing(cfg, sim_dat);
                            
                            % amp corr
                            fprintf('\namplitude envelope correlation...\n')
                            % orthogonalized Brookes et al., 2012, 2014
                            % if z1 and z2 are normalized so that
                            % mean(abs(zi)^2) = 1, then we replace
                            % z2 with z2 - R(c)*z1, where R(c) is
                            % the real part of coherence
                            for m = 1:nRandPhase
                                aec_ortho(m,i,:) = get_aec_ortho_cfn(bp_data.trial{m});
                            end
                            % plv
                            fprintf('\nphaselocking value...\n')
                            % note that this not interpretable for wide
                            % band signals, like high gamma
                            for m = 1:nRandPhase
                                curr_phase = bp_data.trial{m}./abs(bp_data.trial{m});
                                curr_iplv = imag(curr_phase*curr_phase')/size(bp_data.trial{m},2);
                                iplv(m,i,:) = curr_iplv(lower_tri);
                            end
                        end
                        
                    end
                end
            end
        end
    end
end

% correlation with slope
iplv = abs(iplv);
aec_ortho = abs(aec_ortho);

% will becoms columns in dataframe
str_coh = nan(nSim*nBand*nRandPhase,1);
str_aec = nan(nSim*nBand*nRandPhase,1);
str_plv = nan(nSim*nBand*nRandPhase,1);
control_pow = nan(nSim*nBand*nRandPhase,1);
band_order = {};
phase_order = nan(nSim*nBand*nRandPhase,1);
slope = nan(nSim*nBand*nRandPhase,1);

% add strengths for all bands
cnt = 1;
for k = 1:nBand
    for m = 1:nRandPhase
        for j = 2:(nSim + 1)
            % get indices and other stuff
            curr = label{j};
            chan_idx = find(strcmp(label,label{j}));
            elec_idx = cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,1)) &...
                cellfun(@(x) any(strcmp('real', x)), labelcmb(:,2));
            % get strengths
            str_coh(cnt) = mean(Ci(k, elec_idx, m),2);
            str_aec(cnt) = mean(aec_ortho(m, k, elec_idx),2);
            str_plv(cnt) = mean(iplv(m, k, elec_idx),2);
            
            % power in relevant band
            control_pow(cnt) = pow(m, k, chan_idx);
            
            %band name
            band_order(cnt,1) = band_names(k);
            
            % phase
            phase_order(cnt,1) = m;
            
            % slope
            slope(cnt) = str2double(label{j});
            
            % update counter
            cnt = cnt+1;
        end
    end
end

measures = [{'aec'}, {'coh'}, {'plv'}];
reg_data = table(str_aec, str_coh, str_plv, control_pow, band_order, slope, phase_order, 'VariableNames',...
    [{'aec'}, {'coh'}, {'plv'}, {'pow'}, {'band'}, {'slope'}, {'phase'}]);

pvals = zeros(numel(measures),nBand,nRandPhase);
betas = zeros(numel(measures),nBand,nRandPhase);
for m = 1:nRandPhase
    for i = 1:numel(measures)
        for j = 1:nBand
            idx = strcmp(reg_data.band,band_names{j}) & (reg_data.phase == m);
            [b,dev,stats] = glmfit(table2array(reg_data(idx,[4,6])),...
                reg_data.(measures{i})(idx));
            pvals(i,j,m) = stats.p(3);
            betas(i,j,m) = stats.beta(3);
        end
    end
end
saveas(gca, [img_dir, 'p_async.png'], 'png')

figure(1); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(pvals(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        plot([0.05,0.05], [0,.4], 'r'); hold off
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'p_async.png'], 'png')

figure(2); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(betas(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'beta_async.png'], 'png')

%% Repeat with both channels changing slope

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    % just for testing that qsub fuction will work
    for p = 1
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
        for s = 1
            subj = subjects{s};
            % main function for functional connectivity - helps with paralelizing
            
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
            
            fprintf('\n******************************************\nStarting sim for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1 % this should be the same across all sessions
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    data_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    
                    
                    if exist([data_dir, 'channel_info.mat'], 'file')
                        load([data_dir, 'data_clean.mat'])
                        load([data_dir, 'channel_info.mat'])
                        
                        % pick 2 channels, and show IRASA
                        c1 = 1;
                        c2 = 2;
                        
                        spec1 = get_IRASA_spec(ft_data.trial{1}(c1,:), 1, numel(ft_data.trial{1}(c1,:)), ft_data.fsample, win_length, step, filter);
                        spec2 = get_IRASA_spec(ft_data.trial{1}(c2,:), 1, numel(ft_data.trial{1}(c2,:)), ft_data.fsample, win_length, step, filter);
                        
                        nFreq = numel(spec1.freq);
                        x = spec1.freq;
                        y1 = mean(spec1.frac,2);
                        y2 = mean(spec2.frac,2);
                        P = polyfit(log10(x),log10(y2),1);
                        psds = zeros(2,nSim,nFreq);
                        
                        % systematically change 1/f
                        for k = 1:nSim
                            m_new = P(1)*factors(k);
                            b1_new = log10(y1(nFreq/2-1)) - m_new*log10(x(nFreq/2-1));
                            b2_new = log10(y2(nFreq/2-1)) - m_new*log10(x(nFreq/2-1));
                            
                            psds(1,k,:) = 10.^(b1_new + m_new*log10(x)) + (mean(spec1.osci,2));
                            psds(2,k,:) = 10.^(b2_new + m_new*log10(x)) + (mean(spec2.osci,2));
                            
                            figure(2); clf
                            plot(log10(spec1.freq),log10(squeeze(psds(1,k,:))),'b',  'linewidth', 2); hold on
                            plot(log10(x),b1_new + m_new*log10(spec1.freq),'r', 'linewidth', 2)
                            plot(log10(spec2.freq),log10(squeeze(psds(2,k,:))),'k',  'linewidth', 2); hold on
                            plot(log10(x),b2_new + m_new*log10(spec2.freq),'g', 'linewidth', 2);
                        end
                        
                        % get out of frequency space (not exact timeseries,
                        % but with same characteristics
                        new_ts1 = zeros(nSim, length(psds(1,k,:))*2*nRandPhase);
                        for k = 1:nSim
                            dat = [];
                            for m = 1:nRandPhase
                                Sxx = [squeeze(psds(1,k,:))', fliplr(squeeze(psds(1,k,:)))'];
                                curr = TimeseriesFromPSD(Sxx', ft_data.fsample);
                                dat = [curr;dat];
                            end
                            new_ts1(k,:) = dat;
                        end
                        
                        new_ts2 = zeros(nSim, length(psds(2,k,:))*2*nRandPhase);
                        for k = 1:nSim
                            dat = [];
                            for m = 1:nRandPhase
                                Sxx = [squeeze(psds(2,k,:))', fliplr(squeeze(psds(2,k,:)))'];
                                curr = TimeseriesFromPSD(Sxx', ft_data.fsample);
                                dat = [curr;dat];
                            end
                            new_ts2(k,:) = dat;
                        end
                        
                        data = [new_ts1; new_ts2];
                        
                        % regression with power and FC (does slope predict
                        % FC, accounting for power)
                        % get into ft format
                        
                        labels = {};
                        for k = 1:nSim*2
                            if k <= 50
                                labels{k} = ['1_', num2str(P(1)*factors(k))];
                            else
                                labels{k} = ['2_', num2str(P(1)*factors(k-50))];
                            end
                        end
                        trl = [[1:length(Sxx):length(Sxx)*nRandPhase]',...
                            [1:length(Sxx):length(Sxx)*nRandPhase]'+length(Sxx)-1, zeros(nRandPhase,1)];
                        sim_dat = fieldtrip_format(data, ft_data.fsample, labels, trl);
                        
                        % prewhiten
                        cfg = [];
                        cfg.derivative = 'yes';
                        ft_preprocessing(cfg, sim_dat);
                        
                        % band limited
                        fprintf('\nStarting multitaper FFT...\n')
                        % psd and csd - averging done over tapers
                        aec_ortho = zeros(nRandPhase, nBand, nSim);
                        iplv = nan(nRandPhase, nBand, nSim);
                        for l = 1:nSim
                            cfg = [];
                            cfg.method     = 'mtmfft';
                            cfg.taper      = 'dpss';
                            cfg.output     = 'powandcsd';
                            cfg.foi        = freqs;
                            cfg.tapsmofrq  = 4;
                            cfg.pad        = 'maxperlen';
                            cfg.keeptrials  = 'yes';
                            cfg.channel = [labels(l), labels(l+50)];
                            wave = ft_freqanalysis(cfg, sim_dat);
                                                        
                            % power
                            pow = zeros(nRandPhase, nBand, nSim);
                            fprintf('\npower...\n')
                            for i = 1:nBand
                                idx = wave.freq >= bands(i,1) & wave.freq <= bands(i,2);
                                pow(:,i,l) = squeeze(mean(mean(log10(wave.powspctrm(:,:,idx)),3),2));
                            end
                            
                            % coh
                            fprintf('\ncoherence...\n')
                            Ci(:,l,:) = get_coh(wave,bands,'imaginary');
                            
                            % band limited, time resolved
                            fprintf('\nStarting Hilbert transform\n')
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
                                cfg.channel = [labels(l), labels(l+50)];
                                cfg.hilbert = 'complex';
                                
                                bp_data = ft_preprocessing(cfg, sim_dat);
                                
                                % amp corr
                                fprintf('\namplitude envelope correlation...\n')
                                % orthogonalized Brookes et al., 2012, 2014
                                % if z1 and z2 are normalized so that
                                % mean(abs(zi)^2) = 1, then we replace
                                % z2 with z2 - R(c)*z1, where R(c) is
                                % the real part of coherence
                                for m = 1:nRandPhase
                                    aec_ortho(m,i,l) = get_aec_ortho_cfn(bp_data.trial{m});
                                end
                                % plv
                                fprintf('\nphaselocking value...\n')
                                % note that this not interpretable for wide
                                % band signals, like high gamma
                                for m = 1:nRandPhase
                                    curr_phase = bp_data.trial{m}./abs(bp_data.trial{m});
                                    curr_iplv = imag(curr_phase*curr_phase')/size(bp_data.trial{m},2);
                                    iplv(m,i,l) = curr_iplv(2);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% correlation with slope
iplv = abs(iplv);
aec_ortho = abs(aec_ortho);

% will becoms columns in dataframe
str_coh = nan(nSim*nBand*nRandPhase,1);
str_aec = nan(nSim*nBand*nRandPhase,1);
str_plv = nan(nSim*nBand*nRandPhase,1);
control_pow = nan(nSim*nBand*nRandPhase,1);
band_order = {};
phase_order = nan(nSim*nBand*nRandPhase,1);
slope = nan(nSim*nBand*nRandPhase,1);

% add strengths for all bands
cnt = 1;
for k = 1:nBand
    for m = 1:nRandPhase
        for j = 1:nSim
            
            % get strengths
            str_coh(cnt) = Ci(k, j, m);
            str_aec(cnt) = aec_ortho(m, k, j);
            str_plv(cnt) = iplv(m, k, j);
            
            % power in relevant band
            control_pow(cnt) = pow(m, k, j);
            
            %band name
            band_order(cnt,1) = band_names(k);
            
            % phase
            phase_order(cnt,1) = m;
            
            % slope
            slope(cnt) = P(1)*factors(j);
            
            % update counter
            cnt = cnt+1;
        end
    end
end

measures = [{'aec'}, {'coh'}, {'plv'}];
reg_data = table(str_aec, str_coh, str_plv, control_pow, band_order, slope, phase_order, 'VariableNames',...
    [{'aec'}, {'coh'}, {'plv'}, {'pow'}, {'band'}, {'slope'}, {'phase'}]);

pvals = zeros(numel(measures),nBand,nRandPhase);
betas = zeros(numel(measures),nBand,nRandPhase);
for m = 1:nRandPhase
    for i = 1:numel(measures)
        for j = 1:nBand
            idx = strcmp(reg_data.band,band_names{j}) & (reg_data.phase == m);
            [b,dev,stats] = glmfit(table2array(reg_data(idx,[4,6])),...
                reg_data.(measures{i})(idx));
            pvals(i,j,m) = stats.p(3);
            betas(i,j,m) = stats.beta(3);
        end
    end
end

figure(1); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(pvals(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        plot([0.05,0.05], [0,.4], 'r'); hold off
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'p_sync.png'], 'png')

figure(2); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(betas(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'beta_sync.png'], 'png')

%% what if theres an oscillation

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    % just for testing that qsub fuction will work
    for p = 1
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
        for s = 1
            subj = subjects{s};
            % main function for functional connectivity - helps with paralelizing
            
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
            
            fprintf('\n******************************************\nStarting sim for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1 % this should be the same across all sessions
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    data_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    
                    
                    if exist([data_dir, 'channel_info.mat'], 'file')
                        load([data_dir, 'data_clean.mat'])
                        load([data_dir, 'channel_info.mat'])
                        
                        % pick 2 channels, and show IRASA
                        c1 = 1;
                        c2 = 2;
                        
                        spec1 = get_IRASA_spec(ft_data.trial{1}(c1,:), 1, numel(ft_data.trial{1}(c1,:)), ft_data.fsample, win_length, step, filter);
                        spec2 = get_IRASA_spec(ft_data.trial{1}(c2,:), 1, numel(ft_data.trial{1}(c2,:)), ft_data.fsample, win_length, step, filter);
                        
                        % add oscillation
                        spec2.osci(74,:) = spec2.osci(74,:) + 1;
                        figure(1); clf;
                        subplot(2,2,1);
                        loglog(spec1.freq,mean(spec1.mixd,2),'b',  'linewidth', 2); hold on
                        loglog(spec1.freq,mean(spec1.frac,2),'r', 'linewidth', 2);
                        subplot(2,2,2);
                        plot(spec1.freq, mean(spec1.osci,2), 'linewidth', 2); hold on
                        shade_plot(spec1.freq', mean(spec1.osci,2)', (std(spec1.osci,[],2))', rgb("slategrey"), 0.4);
                        subplot(2,2,3);
                        loglog(spec2.freq,mean(spec2.mixd,2),'b',  'linewidth', 2); hold on
                        loglog(spec2.freq,mean(spec2.frac,2),'r', 'linewidth', 2);
                        subplot(2,2,4);
                        plot(spec2.freq, mean(spec2.osci,2), 'linewidth', 2); hold on
                        shade_plot(spec2.freq', mean(spec2.osci,2)', (std(spec2.osci,[],2))', rgb("slategrey"), 0.4);
                        
                        
                        nFreq = numel(spec1.freq);
                        x = spec1.freq;
                        y = mean(spec2.frac,2);
                        y1 = mean(spec1.frac,2);
                        P = polyfit(log10(x),log10(y1),1);
                        psds = zeros(nSim,nFreq);
                        
                        % systematically change 1/f
                        for k = 1:nSim
                            m_new = P(1)*factors(k);
                            b_new = log10(y(nFreq/2-1)) - m_new*log10(x(nFreq/2-1));
                            
                            psds(k,:) = 10.^(b_new + m_new*log10(x) + log10(mean(spec2.osci,2)));
                            figure(2); clf
                            plot(log10(spec2.freq),log10(mean(spec2.mixd,2)),'b',  'linewidth', 2); hold on
                            plot(log10(spec2.freq),log10(mean(spec2.frac,2)),'r', 'linewidth', 2);
                            plot(log10(spec2.freq),log10(psds(k,:)),'k',  'linewidth', 2); hold on
                            plot(log10(x),b_new + m_new*log10(spec2.freq),'g', 'linewidth', 2);
                        end
                        
                        % get out of frequency space (not exact timeseries,
                        % but with same characteristics
                        new_ts = zeros(nSim, length(psds(k,:))*2*nRandPhase);
                        for k = 1:nSim
                            dat = [];
                            for m = 1:nRandPhase
                                Sxx = [psds(k,:), fliplr(psds(k,:))];
                                curr = TimeseriesFromPSD(Sxx', ft_data.fsample);
                                dat = [curr;dat];
                            end
                            new_ts(k,:) = dat;
                        end
                        
                        % regression with power and FC (does slope predict
                        % FC, accounting for power)
                        % get into ft format
                        
                        data = [repmat(ft_data.trial{1}(c1,1:length(Sxx)),1,nRandPhase);new_ts];
                        labels = {'real'};
                        for k = 1:nSim
                            labels{numel(labels) + 1} = num2str(P(1)*factors(k));
                        end
                        trl = [[1:length(Sxx):length(Sxx)*nRandPhase]',...
                            [1:length(Sxx):length(Sxx)*nRandPhase]'+length(Sxx)-1, zeros(nRandPhase,1)];
                        sim_dat = fieldtrip_format(data, ft_data.fsample, labels, trl);
                        
                        % prewhiten
                        cfg = [];
                        cfg.derivative = 'yes';
                        ft_preprocessing(cfg, sim_dat);
                        
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
                        wave = ft_freqanalysis(cfg, sim_dat);
                        
                        % get labels for later
                        label = wave.label;
                        labelcmb = wave.labelcmb;
                        
                        % power
                        pow = zeros(nRandPhase, nBand, nSim+1);
                        fprintf('\npower...\n')
                        for i = 1:nBand
                            idx = wave.freq >= bands(i,1) & wave.freq <= bands(i,2);
                            pow(:,i,:) = squeeze(mean(log10(wave.powspctrm(:,:,idx)),3));
                        end
                        
                        % coh
                        fprintf('\ncoherence...\n')
                        Ci = get_coh(wave,bands,'imaginary');
                        
                        % band limited, time resolved
                        fprintf('\nStarting Hilbert transform\n')
                        aec_ortho = zeros(nRandPhase, nBand, nPair);
                        iplv = nan(nRandPhase, nBand, nPair);
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
                            
                            bp_data = ft_preprocessing(cfg, sim_dat);
                            
                            % amp corr
                            fprintf('\namplitude envelope correlation...\n')
                            % orthogonalized Brookes et al., 2012, 2014
                            % if z1 and z2 are normalized so that
                            % mean(abs(zi)^2) = 1, then we replace
                            % z2 with z2 - R(c)*z1, where R(c) is
                            % the real part of coherence
                            for m = 1:nRandPhase
                                aec_ortho(m,i,:) = get_aec_ortho_cfn(bp_data.trial{m});
                            end
                            % plv
                            fprintf('\nphaselocking value...\n')
                            % note that this not interpretable for wide
                            % band signals, like high gamma
                            for m = 1:nRandPhase
                                curr_phase = bp_data.trial{m}./abs(bp_data.trial{m});
                                curr_iplv = imag(curr_phase*curr_phase')/size(bp_data.trial{m},2);
                                iplv(m,i,:) = curr_iplv(lower_tri);
                            end
                        end
                        
                    end
                end
            end
        end
    end
end

% correlation with slope
iplv = abs(iplv);
aec_ortho = abs(aec_ortho);

% will becoms columns in dataframe
str_coh = nan(nSim*nBand*nRandPhase,1);
str_aec = nan(nSim*nBand*nRandPhase,1);
str_plv = nan(nSim*nBand*nRandPhase,1);
control_pow = nan(nSim*nBand*nRandPhase,1);
band_order = {};
phase_order = nan(nSim*nBand*nRandPhase,1);
slope = nan(nSim*nBand*nRandPhase,1);

% add strengths for all bands
cnt = 1;
for k = 1:nBand
    for m = 1:nRandPhase
        for j = 2:(nSim + 1)
            % get indices and other stuff
            curr = label{j};
            chan_idx = find(strcmp(label,label{j}));
            elec_idx = cellfun(@(x) any(strcmp(x, curr)), labelcmb(:,1)) &...
                cellfun(@(x) any(strcmp('real', x)), labelcmb(:,2));
            % get strengths
            str_coh(cnt) = mean(Ci(k, elec_idx, m),2);
            str_aec(cnt) = mean(aec_ortho(m, k, elec_idx),2);
            str_plv(cnt) = mean(iplv(m, k, elec_idx),2);
            
            % power in relevant band
            control_pow(cnt) = pow(m, k, chan_idx);
            
            %band name
            band_order(cnt,1) = band_names(k);
            
            % phase
            phase_order(cnt,1) = m;
            
            % slope
            slope(cnt) = str2double(label{j});
            
            % update counter
            cnt = cnt+1;
        end
    end
end

measures = [{'aec'}, {'coh'}, {'plv'}];
reg_data = table(str_aec, str_coh, str_plv, control_pow, band_order, slope, phase_order, 'VariableNames',...
    [{'aec'}, {'coh'}, {'plv'}, {'pow'}, {'band'}, {'slope'}, {'phase'}]);

pvals = zeros(numel(measures),nBand,nRandPhase);
betas = zeros(numel(measures),nBand,nRandPhase);
for m = 1:nRandPhase
    for i = 1:numel(measures)
        for j = 1:nBand
            idx = strcmp(reg_data.band,band_names{j}) & (reg_data.phase == m);
            [b,dev,stats] = glmfit(table2array(reg_data(idx,[4,6])),...
                reg_data.(measures{i})(idx));
            pvals(i,j,m) = stats.p(3);
            betas(i,j,m) = stats.beta(3);
        end
    end
end

figure(1); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(pvals(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        plot([0.05,0.05], [0,.4], 'r'); hold off
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'p_async_osc.png'], 'png')
figure(2); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(betas(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'beta_async_osc.png'], 'png')

%% Repeat with both channels changing slope, and oscillations

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    % just for testing that qsub fuction will work
    for p = 1
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
        for s = 1
            subj = subjects{s};
            % main function for functional connectivity - helps with paralelizing
            
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
            
            fprintf('\n******************************************\nStarting sim for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1 % this should be the same across all sessions
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    data_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    
                    
                    if exist([data_dir, 'channel_info.mat'], 'file')
                        load([data_dir, 'data_clean.mat'])
                        load([data_dir, 'channel_info.mat'])
                        
                        % pick 2 channels, and show IRASA
                        c1 = 1;
                        c2 = 2;
                        
                        spec1 = get_IRASA_spec(ft_data.trial{1}(c1,:), 1, numel(ft_data.trial{1}(c1,:)), ft_data.fsample, win_length, step, filter);
                        spec2 = get_IRASA_spec(ft_data.trial{1}(c2,:), 1, numel(ft_data.trial{1}(c2,:)), ft_data.fsample, win_length, step, filter);
                        % add oscillation
                        spec2.osci(74,:) = spec2.osci(74,:) + 1;
                        spec1.osci(74,:) = spec2.osci(74,:) + 1;
                        
                        nFreq = numel(spec1.freq);
                        x = spec1.freq;
                        y1 = mean(spec1.frac,2);
                        y2 = mean(spec2.frac,2);
                        P = polyfit(log10(x),log10(y2),1);
                        psds = zeros(2,nSim,nFreq);
                        
                        % systematically change 1/f
                        for k = 1:nSim
                            m_new = P(1)*factors(k);
                            b1_new = log10(y1(nFreq/2-1)) - m_new*log10(x(nFreq/2-1));
                            b2_new = log10(y2(nFreq/2-1)) - m_new*log10(x(nFreq/2-1));
                            
                            psds(1,k,:) = 10.^(b1_new + m_new*log10(x)) + (mean(spec1.osci,2));
                            psds(2,k,:) = 10.^(b2_new + m_new*log10(x)) + (mean(spec2.osci,2));
                            
                            figure(2); clf
                            plot(log10(spec1.freq),log10(squeeze(psds(1,k,:))),'b',  'linewidth', 2); hold on
                            plot(log10(x),b1_new + m_new*log10(spec1.freq),'r', 'linewidth', 2)
                            plot(log10(spec2.freq),log10(squeeze(psds(2,k,:))),'k',  'linewidth', 2); hold on
                            plot(log10(x),b2_new + m_new*log10(spec2.freq),'g', 'linewidth', 2);
                        end
                        
                        % get out of frequency space (not exact timeseries,
                        % but with same characteristics
                        new_ts1 = zeros(nSim, length(psds(1,k,:))*2*nRandPhase);
                        for k = 1:nSim
                            dat = [];
                            for m = 1:nRandPhase
                                Sxx = [squeeze(psds(1,k,:))', fliplr(squeeze(psds(1,k,:)))'];
                                curr = TimeseriesFromPSD(Sxx', ft_data.fsample);
                                dat = [curr;dat];
                            end
                            new_ts1(k,:) = dat;
                        end
                        
                        new_ts2 = zeros(nSim, length(psds(2,k,:))*2*nRandPhase);
                        for k = 1:nSim
                            dat = [];
                            for m = 1:nRandPhase
                                Sxx = [squeeze(psds(2,k,:))', fliplr(squeeze(psds(2,k,:)))'];
                                curr = TimeseriesFromPSD(Sxx', ft_data.fsample);
                                dat = [curr;dat];
                            end
                            new_ts2(k,:) = dat;
                        end
                        
                        data = [new_ts1; new_ts2];
                        
                        % regression with power and FC (does slope predict
                        % FC, accounting for power)
                        % get into ft format
                        
                        labels = {};
                        for k = 1:nSim*2
                            if k <= 50
                                labels{k} = ['1_', num2str(P(1)*factors(k))];
                            else
                                labels{k} = ['2_', num2str(P(1)*factors(k-50))];
                            end
                        end
                        trl = [[1:length(Sxx):length(Sxx)*nRandPhase]',...
                            [1:length(Sxx):length(Sxx)*nRandPhase]'+length(Sxx)-1, zeros(nRandPhase,1)];
                        sim_dat = fieldtrip_format(data, ft_data.fsample, labels, trl);
                        
                        % prewhiten
                        cfg = [];
                        cfg.derivative = 'yes';
                        ft_preprocessing(cfg, sim_dat);
                        
                        % band limited
                        fprintf('\nStarting multitaper FFT...\n')
                        % psd and csd - averging done over tapers
                        aec_ortho = zeros(nRandPhase, nBand, nSim);
                        iplv = nan(nRandPhase, nBand, nSim);
                        for l = 1:nSim
                            cfg = [];
                            cfg.method     = 'mtmfft';
                            cfg.taper      = 'dpss';
                            cfg.output     = 'powandcsd';
                            cfg.foi        = freqs;
                            cfg.tapsmofrq  = 4;
                            cfg.pad        = 'maxperlen';
                            cfg.keeptrials  = 'yes';
                            cfg.channel = [labels(l), labels(l+50)];
                            wave = ft_freqanalysis(cfg, sim_dat);
                                                        
                            % power
                            pow = zeros(nRandPhase, nBand, nSim);
                            fprintf('\npower...\n')
                            for i = 1:nBand
                                idx = wave.freq >= bands(i,1) & wave.freq <= bands(i,2);
                                pow(:,i,l) = squeeze(mean(mean(log10(wave.powspctrm(:,:,idx)),3),2));
                            end
                            
                            % coh
                            fprintf('\ncoherence...\n')
                            Ci(:,l,:) = get_coh(wave,bands,'imaginary');
                            
                            % band limited, time resolved
                            fprintf('\nStarting Hilbert transform\n')
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
                                cfg.channel = [labels(l), labels(l+50)];
                                cfg.hilbert = 'complex';
                                
                                bp_data = ft_preprocessing(cfg, sim_dat);
                                
                                % amp corr
                                fprintf('\namplitude envelope correlation...\n')
                                % orthogonalized Brookes et al., 2012, 2014
                                % if z1 and z2 are normalized so that
                                % mean(abs(zi)^2) = 1, then we replace
                                % z2 with z2 - R(c)*z1, where R(c) is
                                % the real part of coherence
                                for m = 1:nRandPhase
                                    aec_ortho(m,i,l) = get_aec_ortho_cfn(bp_data.trial{m});
                                end
                                % plv
                                fprintf('\nphaselocking value...\n')
                                % note that this not interpretable for wide
                                % band signals, like high gamma
                                for m = 1:nRandPhase
                                    curr_phase = bp_data.trial{m}./abs(bp_data.trial{m});
                                    curr_iplv = imag(curr_phase*curr_phase')/size(bp_data.trial{m},2);
                                    iplv(m,i,l) = curr_iplv(2);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% correlation with slope
iplv = abs(iplv);
aec_ortho = abs(aec_ortho);

% will becoms columns in dataframe
str_coh = nan(nSim*nBand*nRandPhase,1);
str_aec = nan(nSim*nBand*nRandPhase,1);
str_plv = nan(nSim*nBand*nRandPhase,1);
control_pow = nan(nSim*nBand*nRandPhase,1);
band_order = {};
phase_order = nan(nSim*nBand*nRandPhase,1);
slope = nan(nSim*nBand*nRandPhase,1);

% add strengths for all bands
cnt = 1;
for k = 1:nBand
    for m = 1:nRandPhase
        for j = 1:nSim
            
            % get strengths
            str_coh(cnt) = Ci(k, j, m);
            str_aec(cnt) = aec_ortho(m, k, j);
            str_plv(cnt) = iplv(m, k, j);
            
            % power in relevant band
            control_pow(cnt) = pow(m, k, j);
            
            %band name
            band_order(cnt,1) = band_names(k);
            
            % phase
            phase_order(cnt,1) = m;
            
            % slope
            slope(cnt) = P(1)*factors(j);
            
            % update counter
            cnt = cnt+1;
        end
    end
end

measures = [{'aec'}, {'coh'}, {'plv'}];
reg_data = table(str_aec, str_coh, str_plv, control_pow, band_order, slope, phase_order, 'VariableNames',...
    [{'aec'}, {'coh'}, {'plv'}, {'pow'}, {'band'}, {'slope'}, {'phase'}]);

pvals = zeros(numel(measures),nBand,nRandPhase);
betas = zeros(numel(measures),nBand,nRandPhase);
for m = 1:nRandPhase
    for i = 1:numel(measures)
        for j = 1:nBand
            idx = strcmp(reg_data.band,band_names{j}) & (reg_data.phase == m);
            [b,dev,stats] = glmfit(table2array(reg_data(idx,[4,6])),...
                reg_data.(measures{i})(idx));
            pvals(i,j,m) = stats.p(3);
            betas(i,j,m) = stats.beta(3);
        end
    end
end

figure(1); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(pvals(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        plot([0.05,0.05], [0,.4], 'r'); hold off
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'p_sync_osc.png'], 'png')
figure(2); clf
cnt = 1;
for i = 1:numel(measures)
    for j = 1:nBand
        subplot(numel(measures),nBand,cnt)
        histogram(betas(i,j,:),20, 'Normalization', 'probability'); hold on
        title([band_names{j}, ' ', measures{i}])
        cnt = cnt + 1;
    end
end
saveas(gca, [img_dir, 'beta_sync_osc.png'], 'png')

clear
clc
close all
warning ON

%%

% global variables and packages
top_dir = 'smb://bassett-data.seas.upenn.edu/bassett-data/Jeni/RAM/';
eval(['cd ', top_dir])

% for removing electrodes
thr = 1.5;
releases = ['1', '2', '3'];

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    eval(['cd ', release_dir '/protocols'])
    
    % remove parent and hidden directories, then get protocols
    folders = dir(pwd);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    for p = 1:numel(protocols)
        protocol = protocols{p};
        
        % get global info struct
        fname = [protocol, '.json'];
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
            eval(['diary ', [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/log.txt']]);
            
            fprintf('******************************************\nStarting preprocessing for subject %s...\n', subj)
            
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
                    
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                     % folders
                    raw_dir = [release_dir, 'protocols/', protocol, '/subjects/', subj, '/experiments/', exper, '/sessions/', sess, '/ephys/current_processed/'];
                    save_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    img_dir = [top_dir, 'img/diagnostics/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    if ~exist(img_dir, 'dir')
                        mkdir(img_dir);
                    end
                    if ~exist(save_dir, 'dir')
                        mkdir(save_dir);
                    end

                    fprintf('\nExperiment %s session %s\n\n', exper, sess)

                    %% get channel info
                    
                    % get file info struct
                    fid = fopen([release_dir, curr_info.contacts]);
                    raw = fread(fid);
                    channel_info = jsondecode(char(raw'));
                    code = fields(channel_info); % sometimes this doen't match subject
                    eval(['channel_info = channel_info.',  code{1}, ';'])
                    fclose(fid);
                    
                    % get numbers
                    labels = fields(channel_info.contacts);
                    nChan = numel(labels);
                    chann_idx = zeros(nChan, 1);
                    for i = 1:numel(labels)
                        chann = labels{i};
                        eval(['chann_idx(i) = channel_info.contacts.', chann, '.channel;'])
                    end
                    % sort
                    [chann_idx, sort_idx] = sort(chann_idx, 'ascend');
                    labels = labels(sort_idx);
                    
                    
                    %% Get data
                    
                    % get header
                    fname = [raw_dir, 'sources.json'];
                    fid = fopen(fname);
                    raw = fread(fid,inf);
                    str = char(raw');
                    fclose(fid);
                    val = jsondecode(str);
                    
                    % get event info
                    if isfield(curr_info, 'all_events')
                        fid = fopen([release_dir, curr_info.all_events]);
                    else % some only have task events
                        fid = fopen([release_dir, curr_info.task_events]);
                    end
                    raw = fread(fid);
                    events = jsondecode(char(raw'));
                    eegfile = unique({events.eegfile});
                    eegfile = eegfile(cellfun(@(x) ~isempty(x), eegfile));

                    % check that there's only one file per session
                    if numel(eegfile) > 1
                       error('This session has multiple eegfiles') 
                    end
                    eval(['header = val.', eegfile{1} ';'])
                    
                    % get data
                    try
                        data_raw = zeros(nChan, header.n_samples-1);
                        for i = 1:nChan
                            chan = sprintf('%03d', chann_idx(i));
                            fid = fopen([raw_dir, 'noreref/', eegfile{1}, '.', chan]);
                            data_raw(i,:) = fread(fid, header.n_samples, header.data_format)';
                            fclose(fid);
                        end
                    catch
                        % sometimes the size of the file and header.n_samples
                        % don't match. If that is the case, keep the file
                        % format listed, and just don't specify a size. also
                        % print a warning
                        warning('The number of samples indicated inthe header and actual number of samples at the given precision do not match')
                        
                        data_raw = [];
                        for i = 1:nChan
                            chan = sprintf('%03d', chann_idx(i));
                            fid = fopen([raw_dir, 'noreref/', eegfile{1}, '.', chan]);
                            data_raw(i,:) = fread(fid, inf, header.data_format)';
                            fclose(fid);
                        end
                        fprintf('Updating header\n')
                        header.n_samples = size(data_raw,2);
                    end
                    
                    
                    %% Epoch?
                    
                    % get pre and post task data
                    data = [];
                    
                    switch exper(1:end-1)
                        case {'FR', 'CatFR', 'PAL'}
                            start_idx = find(strcmpi('sess_start', {events.type}));
                            end_idx = find(strcmpi('sess_end', {events.type}));

                        case {'YC', 'TH'}
                            % these tasks don't have "sess_start" and end,
                            % so we are just going to lok for the first and
                            % last task event (when eegoffset is greater
                            % than 0)
                            start_idx = find([events.eegoffset] > 1, 1);
                            end_idx = find([events.eegoffset] > 1, 1, 'last');
                    end
                    
                    % make sure there is enough data, and that the
                    % event was recorded
                    if ~isempty(start_idx)
                        if numel(1:(events(start_idx).eegoffset - 1)) > 1000
                            data{1} = data_raw(:, 1:(events(start_idx).eegoffset - 1));
                            fprintf('%d samples of pre-task data\n', (events(start_idx).eegoffset - 1))
                        else
                            fprintf('0 samples of pre-task data\n')
                        end
                    else
                        fprintf('0 samples of pre-task data\n')
                    end
                    
                    
                    % now the post task data
                    if ~isempty(end_idx)
                        if numel(data) > 0
                            if numel((events(end_idx).eegoffset + 1):header.n_samples) > 1000
                                data{2} = data_raw(:, (events(end_idx).eegoffset + 1):end);
                                fprintf('%d samples of post-task data\n', numel((events(end_idx).eegoffset + 1):header.n_samples))
                            else
                                fprintf('0 samples of post-task data\n')
                            end
                        else % is no start data, index at 1
                            if numel((events(end_idx).eegoffset + 1):header.n_samples) > 1000
                                data{1} = data_raw(:, (events(end_idx).eegoffset + 1):end);
                                fprintf('%d samples of post-task data\n', numel((events(end_idx).eegoffset + 1):header.n_samples))
                            else
                                fprintf('0 samples of post-task data\n')
                            end
                        end
                    else
                        fprintf('0 samples of post-task data\n')
                    end
                 
                    
                    %% Filter
                    fprintf('\nFiltering...')
                    
                    for j = 1:numel(data)
                        
                        % filter out 60 Hz harmonics
                        fprintf('60...')
                        [b,a] = butter(4, [59/(header.sample_rate/2), 61/(header.sample_rate/2)], 'stop');
                        data{j} = filtfilt(b,a,data{j}')';
                        
                        fprintf('120...')
                        [b,a] = butter(4, [119/(header.sample_rate/2), 121/(header.sample_rate/2)], 'stop');
                        data{j} = filtfilt(b,a,data{j}')';
                        
                        fprintf('180...')
                        [b,a] = butter(4, [179/(header.sample_rate/2), 181/(header.sample_rate/2)], 'stop');
                        data{j} = filtfilt(b,a,data{j}')';
                    end
                    
                    fprintf(' done!\n\n')
                    
                    %% Remove bad elecs
                    
                    % find bad elecs from documentation
                    fid = fopen([metadata_dir, 'electrode_categories/', subj, '_electrode_categories.txt'], 'r');
                    elec_cats = textscan(fid, '%s', 'Delimiter', '\n\n');
                    fclose(fid);
                    elec_cats = elec_cats{1}(~cellfun('isempty',elec_cats{1}));
                    % split by categories
                    soz_ind = find(cell2mat(cellfun(@(x) contains('seizure onset',x), elec_cats, 'UniformOutput', false)));
                    interictal_ind = find(cell2mat(cellfun(@(x) contains('interictal spikes',x), elec_cats, 'UniformOutput', false)));
                    bad_contact_ind = find(cell2mat(cellfun(@(x) contains('bad electrodes',x), elec_cats, 'UniformOutput', false)));
                    % check size
                    if (numel(soz_ind) > 1 || numel(interictal_ind) > 1 || numel(bad_contact_ind) > 1)
                        error('Something wen wrong with parsing the electrode categories');
                    end
                    soz = elec_cats((soz_ind + 1):(interictal_ind - 1));
                    ict_spike = elec_cats((interictal_ind + 1):(bad_contact_ind - 1));
                    bad_cont = elec_cats((bad_contact_ind + 1):end);
                    % actually removing the bad elecs
                    bad_cont_idx = cell2mat(cellfun(@(x) any(strcmp(x,bad_cont)), labels, 'UniformOutput',false));
                    for j = 1:numel(data)
                        data{j} = data{j}(~bad_cont_idx,:);
                    end
                    chann_idx = chann_idx(~bad_cont_idx);
                    labels = labels(~bad_cont_idx);
                    
                    % additional checking
                    fprintf('Removing contacts with kurtosis, or power greater than %d standard deviations above average, or line length greater than 3x the mean...', thr)
                    % get bad elecs across either session
                    rmv = false(size(labels));
                    for j = 1:numel(data)
                        rmv = rmv | reject_elecs(data{j}, thr, header.sample_rate);
                    end
                    for j = 1:numel(data)
                        data{j} = data{j}(~rmv, :);
                    end
                    chann_idx = chann_idx(~rmv);
                    labels = labels(~rmv);
                    
                    fprintf('done!\nRemoved %d contacts\n\n', sum(rmv))
                    
                    % make some plots pre CAR - if something looks weird might be best to
                    % remove this subject, because the extra noise will have affected the car
                    % (unless noise is local to electrode/grids)
                    for j = 1:numel(data)
                        dur = size(data{j},2)/header.sample_rate; % duration for plotting
                        
                        figure(1); clf;
                        fft_plot(data{j}', 1000, header.sample_rate);
                        saveas(gca, [img_dir, 'preCAR_FFT_', num2str(j), '.png'], 'png')
                        
                        cnt = 1;
                        ext = 1;
                        while cnt < dur
                            if (cnt+(100*header.sample_rate)-1) > dur
                                eegplot(data{j}(:, cnt:end), 'srate', header.sample_rate,  "winlength", 100);
                                saveas(gca, [img_dir, 'preCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            else
                                eegplot(data{j}(:, cnt:(cnt+(100*header.sample_rate)-1)), 'srate', header.sample_rate,  "winlength", 100);
                                saveas(gca, [img_dir, 'preCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            end
                            cnt = cnt + 100*header.sample_rate;
                            ext = ext + 1;
                        end
                        close all
                    end
                    
                    % save elec info
                    save([save_dir, 'channel_info.mat'], 'channel_info', 'soz', 'ict_spike', 'labels', 'chann_idx');

                    
                    %% CAR
                    
                    for j = 1:numel(data)
                        
                        data{j} = data{j} - mean(data{j},2); % demean
                        data{j} = detrend(data{j}')'; % as opposed to low pass filtering
                        data{j} = get_CAR(data{j}, labels); % CAR by group
                        
                        fprintf('Finished CAR!\n')
                        dur = size(data{j},2); % duration for plotting in samples
                        figure(2); clf;
                        fft_plot(data{j}', 1000, header.sample_rate);
                        saveas(gca, [img_dir, 'postCAR_FFT_', num2str(j), '.png'], 'png')
                        
                        % plot chunks of 100s
                        cnt = 1;
                        ext = 1;
                        while cnt < dur
                            if (cnt+(100*header.sample_rate)-1) > dur
                                eegplot(data{j}(:, cnt:end), 'srate', header.sample_rate,  "winlength", 100);
                                saveas(gca, [img_dir, 'postCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            else
                                eegplot(data{j}(:, cnt:(cnt+(100*header.sample_rate)-1)), 'srate', header.sample_rate,  "winlength", 100);
                                saveas(gca, [img_dir, 'postCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            end
                            cnt = cnt + 100*header.sample_rate;
                            ext = ext + 1;
                        end
                        close all
                        
                    end
                    
                    % save data
                    %save([save_dir, 'data_clean.mat'], 'data')
                    
                    %% Spikes
                    
                    
                    
                    % save other things
                    save([save_dir, 'header.mat'], 'header')
                end
            end
            diary off
        end
    end
end
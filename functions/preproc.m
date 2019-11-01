function [warnings] = preproc(thr, release_dir, top_dir, release, protocols, warnings, errors)
% main preprocessing script for RAM data - helps with parralelizing
eval(['cd ', release_dir '/protocols'])

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
        eval(['diary ', [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/log.txt']]);
        
        fprintf('******************************************\nStarting preprocessing for subject %s...\n', subj)
        
        % get experiements
        eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
        for e = 1:numel(experiments)
            exper = experiments{e};
            
            % get seesions
            eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
            for n = 1:numel(sessions)
                try
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
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
                    
                    % get file info struct - this file has names, channel index, and positions
                    % in space. It does not have categories (SOZ, interictal, etc)
                    fid = fopen([release_dir, curr_info.contacts]);
                    raw = fread(fid);
                    channel_info = jsondecode(char(raw'));
                    code = fields(channel_info); % sometimes this doen't match subject
                    eval(['channel_info = channel_info.',  code{1}, ';'])
                    fclose(fid);
                    
                    % get numbers so that you can load the correct files. Don't assume they are
                    % in the same order (they usually arent)
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
                    
                    % check if these channels match what is in the
                    % folder
                    all_eegfiles = dir([raw_dir, 'noreref/' subj, '_*']);
                    all_eegfiles = {all_eegfiles(:).name};
                    % get channels of files you actually have
                    all_channels = cellfun(@(x) strsplit(x, '.'), all_eegfiles, 'UniformOutput', false);
                    all_channels = cellfun(@(x) str2double(x{end}), all_channels);
                    % match them to labels
                    idx = false(size(chann_idx));
                    for k = 1:numel(idx)
                        curr = chann_idx(k);
                        % add only matching files, skip the other
                        % ones
                        idx(k) =  any(curr == all_channels);
                    end
                    % update
                    labels = labels(idx);
                    chann_idx = chann_idx(idx);
                    nChan = numel(labels);
                    
                    %% Get event info
                    
                    % load event info
                    if isfield(curr_info, 'all_events') && exist([release_dir, curr_info.all_events],'file')
                        fid = fopen([release_dir, curr_info.all_events]);
                    elseif isfield(curr_info, 'task_events') && exist([release_dir, curr_info.task_events],'file')% some only have tasks events
                        fid = fopen([release_dir, curr_info.task_events]);
                    elseif exist([release_dir, 'protocols/', protocol, '/subjects/', subj, '/experiments/', exper, '/sessions/', sess, '/behavioral/current_processed/all_events.json'],'file') %if info has the wrong file
                        fid = fopen([release_dir, 'protocols/', protocol, '/subjects/', subj, '/experiments/', exper, '/sessions/', sess, '/behavioral/current_processed/all_events.json']);
                    elseif exist([release_dir, 'protocols/', protocol, '/subjects/', subj, '/experiments/', exper, '/sessions/', sess, '/behavioral/current_processed/task_events.json'],'file') %if info has the wrong file
                        fid = fopen([release_dir, 'protocols/', protocol, '/subjects/', subj, '/experiments/', exper, '/sessions/', sess, '/behavioral/current_processed/task_events.json']);
                    else
                        error('Could not find an events file')
                    end
                    raw = fread(fid);
                    events = jsondecode(char(raw'));
                    
                    % get pre and post task data time points
                    switch lower(exper(1:end-1))
                        case {'fr', 'catfr', 'pal'}
                            start_idx = find(strcmpi('sess_start', {events.type}));
                            end_idx = find(strcmpi('sess_end', {events.type}));
                            
                        case {'yc', 'th'}
                            % these tasks don't have "sess_start" and end,
                            % so we are just going to look for the first and
                            % last task event (when eegoffset is greater
                            % than 0, if its less then this event is the end of the recording)
                            start_idx = find([events.eegoffset] > 1, 1);
                            end_idx = find([events.eegoffset] > 1, 1, 'last');
                    end
                    % if multiple, take the ends. Not sure this ever actually happens
                    if ~isempty(start_idx)
                        start_idx = start_idx(1);
                    end
                    if ~isempty(end_idx)
                        end_idx = end_idx(end);
                    end
                    
                    % make flags for pre and post task data - want to know if we have enough
                    % data to work with
                    pretask = ~isempty(start_idx); % we will check the length later
                    posttask = ~isempty(end_idx) && events(end_idx).eegoffset > 0; % we'll check the length of this one later
                    
                    % get eegfiles for pre and post task data
                    if pretask && posttask
                        pre_eeg = events(start_idx).eegfile;
                        post_eeg = events(end_idx).eegfile;
                        if strcmp(pre_eeg, post_eeg)
                            eegfile = pre_eeg;
                        else
                            % check that there's only one file per session - right now the code
                            % isn't equipped to handle multiple files
                            warning('This session has multiple eegfiles.')
                            warnings(end+1).files = [subj, '_', exper, '_', sess];
                            warnings(end).message = 'This session has multiple eegfiles';
                            
                            % pretask typically has more data, so we're
                            % going to pick that one
                            eegfile = pre_eeg;
                        end
                    elseif pretask
                        pre_eeg = events(start_idx).eegfile;
                        eegfile = pre_eeg;
                    elseif posttask
                        post_eeg = events(end_idx).eegfile;
                        eegfile = post_eeg;
                    end
                    
                    %% Get data
                    
                    % get header
                    fname = [raw_dir, 'sources.json'];
                    fid = fopen(fname);
                    raw = fread(fid,inf);
                    str = char(raw');
                    fclose(fid);
                    val = jsondecode(str);
                    
                    try
                        eval(['header = val.', eegfile ';'])
                    catch % sometimes the eegfile from the events is missing some strings (i.e. experiment, etc)
                        warning('The eegfile in events does not match the eegfile in the data')
                        warnings(end+1).files = [subj, '_', exper, '_', sess];
                        warnings(end).message = 'The eegfile in events does not match the eegfile in the data';
                        
                        eegfile = fields(val);
                        if numel(eegfile) == 1
                            eegfile = eegfile{1};
                        else
                            % check that there's only one file per session
                            error('This session has multiple eegfiles in header, and none of them are consistent with the events.')
                        end
                        eval(['header = val.', eegfile ';'])
                    end
                    min_seg = 5*header.sample_rate; % minimum duration to keep data - this is how much data the spike detection needs
                    
                    % get data
                    try
                        data_raw = zeros(nChan, header.n_samples-1);
                        % data is saved per channel
                        for i = 1:nChan
                            chan = sprintf('%03d', chann_idx(i));
                            fid = fopen([raw_dir, 'noreref/', eegfile, '.', chan]);
                            data_raw(i,:) = fread(fid, header.n_samples, header.data_format)';
                            fclose(fid);
                        end
                        header.n_samples = header.n_samples - 1;
                    catch
                        % sometimes the size of the file and header.n_samples
                        % don't match. If that is the case, keep the file
                        % format listed, and just don't specify a size. also
                        % print a warning
                        warning('The number of samples indicated in the header and actual number of samples at the given precision do not match')
                        warnings(end+1).files = [subj, '_', exper, '_', sess];
                        warnings(end).message = 'The number of samples indicated in the header and actual number of samples at the given precision do not match';
                        
                        data_raw = [];
                        for i = 1:nChan
                            chan = sprintf('%03d', chann_idx(i));
                            fid = fopen([raw_dir, 'noreref/', eegfile, '.', chan]);
                            data_raw(i,:) = fread(fid, inf, header.data_format)';
                            fclose(fid);
                        end
                        fprintf('Updating header\n')
                        header.n_samples = size(data_raw,2);
                    end
                    
                    
                    %% Epoch?
                    
                    % check that we have enought data, now
                    % that we know how long the data is
                    if posttask
                        posttask = numel((events(end_idx).eegoffset + 1):header.n_samples) > min_seg;
                    end
                    if pretask
                        pretask = numel(1:(events(start_idx).eegoffset - 1)) > min_seg;
                    end
                    
                    if ~(pretask || posttask)
                        error('This subject does not have enough task free data..stopping processing')
                    end
                    
                    
                    % make trl structre for fieldtrip: start, end,
                    % offset in samples
                    trl = [];
                    if pretask
                        trl(1, 1) = 1;
                        trl(1, 2) = events(start_idx).eegoffset - 1;
                        fprintf('%d samples of pre-task data\n', (events(start_idx).eegoffset - 1))
                    else
                        fprintf('0 samples of pre-task data\n')
                    end
                    
                    
                    % now the post task data
                    if posttask
                        if numel(trl) > 0
                            trl(2, 1) = events(end_idx).eegoffset + 1;
                            trl(2, 2) = header.n_samples;
                            fprintf('%d samples of post-task data\n', numel((events(end_idx).eegoffset + 1):header.n_samples))
                        else
                            trl(1, 1) = events(end_idx).eegoffset + 1;
                            trl(1, 2) = header.n_samples;
                            fprintf('%d samples of post-task data\n', numel((events(end_idx).eegoffset + 1):header.n_samples))
                        end
                    else
                        fprintf('0 samples of post-task data\n')
                    end
                    trl(:,3) = 0;
                    
                    % get in fieldtrip format
                    ft_data = fieldtrip_format(data_raw, header.sample_rate, labels, trl);
                    nTrial = numel(ft_data.trial);
                    
                    %% Downsample
                    % to 512, smallest number for all subjects
                    fprintf('\nResampling...\n')
                    cfg = [];
                    cfg.resamplefs = 512;
                    ft_data = ft_resampledata(cfg,ft_data);
                    header.sample_rate = 512;
                    
                    %% Filter
                    fprintf('\nFiltering...')
                    
                    for j = 1:nTrial
                        
                        % filter out 60 Hz harmonics
                        fprintf('60...')
                        [b,a] = butter(4, [59/(header.sample_rate/2), 61/(header.sample_rate/2)], 'stop');
                        ft_data.trial{j} = filtfilt(b,a,ft_data.trial{j}')';
                        
                        fprintf('120...')
                        [b,a] = butter(4, [119/(header.sample_rate/2), 121/(header.sample_rate/2)], 'stop');
                        ft_data.trial{j} = filtfilt(b,a,ft_data.trial{j}')';
                        
                        fprintf('180...')
                        [b,a] = butter(4, [179/(header.sample_rate/2), 181/(header.sample_rate/2)], 'stop');
                        ft_data.trial{j} = filtfilt(b,a,ft_data.trial{j}')';
                    end
                    
                    fprintf(' done!\n\n')
                    
                    %% Remove bad elecs
                    
                    % find bad elecs (as well as SOZ and interictal) from documentation in
                    % metadata
                    try
                        fid = fopen([metadata_dir, 'electrode_categories/', subj, '_electrode_categories.txt'], 'r');
                        elec_cats = textscan(fid, '%s', 'Delimiter', '\n\n');
                        fclose(fid);
                    catch % 1st release has a different naming scheme
                        fid = fopen([metadata_dir, 'electrode_categories/electrode_categories_', subj, '.txt'], 'r');
                        elec_cats = textscan(fid, '%s', 'Delimiter', '\n\n');
                        fclose(fid);
                    end
                    elec_cats = elec_cats{1}(~cellfun('isempty',elec_cats{1}));
                    
                    % split by categories - they are not uniform across sites and releases, but
                    % this should catch all of them
                    soz_ind = find(cell2mat(cellfun(@(x) contains(lower(x), [{'seizure'}, {'onset'}]), elec_cats, 'UniformOutput', false)));
                    interictal_ind = find(cell2mat(cellfun(@(x) contains(lower(x), 'interictal'), elec_cats, 'UniformOutput', false)));
                    bad_contact_ind = find(cell2mat(cellfun(@(x) ... % this one is complicated, some subjects have multiple categories of bad elecs
                        contains(lower(x), [{'bad'}, {'broken'} {'lesion'}]), elec_cats, 'UniformOutput', false)));
                    % check size - this is only a warning because some sites don't list SOZ,
                    % and some subjects don't have bad elecs marked...presumably because the
                    % recording was good. However this might also happen if something went
                    % wrong with the parsing
                    if (numel(soz_ind) ~= 1 || numel(interictal_ind) ~= 1 || numel(bad_contact_ind) < 1)
                        warning('Something went wrong with parsing the electrode categories, or this subject is missing categories.');
                        warnings(end+1).files = [subj, '_', exper, '_', sess];
                        warnings(end).message = 'Something went wrong with parsing the electrode categories, or this subject is missing categories.';
                    end
                    
                    % get contact names for each category
                    soz = elec_cats((soz_ind + 1):(interictal_ind - 1));
                    interictal_cont = elec_cats((interictal_ind + 1):(bad_contact_ind - 1));
                    if ~isempty(bad_contact_ind)
                        bad_cont = elec_cats((bad_contact_ind(1) + 1):end);
                        % remove empty categories and category names for bad electrodes
                        bad_cont(cell2mat(cellfun(@(x) ...
                            contains(lower(x), [{'none'}, {'bad'}, {'broken'}, {'lesion'}, {'-'}, {'*'}, {' '}] ),... % the space takes out notes
                            bad_cont, 'UniformOutput', false))) = [];
                    else
                        bad_cont = {};
                    end
                    
                    
                    % check that the order didn't get messed up - these sums should be the same
                    if (numel(soz_ind) + numel(interictal_ind) + numel(bad_contact_ind) + numel(soz) + numel(interictal_cont) + numel(bad_cont)) > (numel(elec_cats) - 1) % accounting for subj at the beginning
                        error('The order of electrodes categories on this subject was messed up')
                    end
                    
                    % actually removing the bad elecs
                    bad_cont_idx = cell2mat(cellfun(@(x) any(strcmp(x,bad_cont)), ft_data.label, 'UniformOutput',false));
                    for j = 1:nTrial
                        ft_data.trial{j} = ft_data.trial{j}(~bad_cont_idx,:);
                    end
                    chann_idx = chann_idx(~bad_cont_idx);
                    ft_data.label = ft_data.label(~bad_cont_idx);
                    
                    
                    % additional checking with my own algorthim. This will pick up elecs that
                    % are only bad during this particular chunk of recording. Also I don't know
                    % what their rejection criteria were. these are fairly strict
                    fprintf('Removing contacts with kurtosis, or power greater than %d standard deviations above average, or line length greater than 3x the mean...', thr)
                    % get bad elecs across either session
                    rmv = false(size(ft_data.label));
                    for j = 1:nTrial
                        rmv = rmv | reject_elecs(ft_data.trial{j}, thr, header.sample_rate);
                    end
                    for j = 1:nTrial
                        ft_data.trial{j} = ft_data.trial{j}(~rmv, :);
                    end
                    chann_idx = chann_idx(~rmv);
                    ft_data.label = ft_data.label(~rmv);
                    
                    fprintf('done!\nRemoved %d contacts\n\n', sum(rmv))
                    
                    % make some plots pre CAR - if something looks really weird might be best to
                    % remove this subject, because the extra noise will have affected the car
                    % (unless noise is local to specific electrode/grids)
                    for j = 1:nTrial
                        dur = size(ft_data.trial{j},2)/header.sample_rate; % duration for plotting
                        
                        figure(1); clf;
                        fft_plot(ft_data.trial{j}', 1000, header.sample_rate);
                        saveas(gcf, [img_dir, 'preCAR_FFT_', num2str(j), '.png'], 'png')
                        
                        cnt = 1;
                        ext = 1;
                        while cnt < dur
                            if (cnt+(100*header.sample_rate)-1) > dur
                                plot_lfp(ft_data.trial{j}(:, cnt:end), header.sample_rate, 2);
                                saveas(gcf, [img_dir, 'preCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            else
                                plot_lfp(ft_data.trial{j}(:, cnt:(cnt+(100*header.sample_rate)-1)), header.sample_rate, 2);
                                saveas(gcf, [img_dir, 'preCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            end
                            cnt = cnt + 100*header.sample_rate;
                            ext = ext + 1;
                        end
                        close all
                    end
                    
                    % save elec info
                    save([save_dir, 'channel_info.mat'], 'channel_info', 'soz', 'interictal_cont', 'labels', 'chann_idx');
                    
                    
                    %% CAR
                    
                    for j = 1:nTrial
                        
                        ft_data.trial{j} = ft_data.trial{j} - mean(ft_data.trial{j},2); % demean
                        ft_data.trial{j} = detrend(ft_data.trial{j}')'; % as opposed to low pass filtering
                        ft_data.trial{j} = get_CAR(ft_data.trial{j}, ft_data.label); % CAR by group
                        
                        fprintf('Finished CAR!\n')
                        dur = size(ft_data.trial{j},2); % duration for plotting in samples
                        figure(2); clf;
                        fft_plot(ft_data.trial{j}', 1000, header.sample_rate);
                        saveas(gcf, [img_dir, 'postCAR_FFT_', num2str(j), '.png'], 'png')
                        
                        % plot chunks of 100s
                        cnt = 1;
                        ext = 1;
                        while cnt < dur
                            if (cnt+(100*header.sample_rate)-1) > dur
                                plot_lfp(ft_data.trial{j}(:, cnt:end), header.sample_rate, 2);
                                saveas(gcf, [img_dir, 'postCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            else
                                plot_lfp(ft_data.trial{j}(:, cnt:(cnt+(100*header.sample_rate)-1)), header.sample_rate, 2);
                                saveas(gcf, [img_dir, 'postCAR_LFP_', num2str(j), '_', num2str(ext), '.png'], 'png')
                            end
                            cnt = cnt + 100*header.sample_rate;
                            ext = ext + 1;
                        end
                        close all
                        
                    end
                    
                    % save data and other things
                    save([save_dir, 'data_clean.mat'], 'ft_data')
                    save([save_dir, 'header.mat'], 'header')
                catch ME
                    errors(end+1).files = [subj, '_', exper, '_', sess];
                    errors(end).message = ME.message;
                end
            end
        end
        diary off
    end
end


% remove empty entry
errors = errors(2:end);
warnings = warnings(2:end);
save([release_dir, 'errors.mat'], 'errors');
save([release_dir, 'warnings.mat'], 'warnings');


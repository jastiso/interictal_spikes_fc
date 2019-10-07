clear
clc
close

% global variables and packages
addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

top_dir = '/Volumes/bassett-data/Jeni/RAM/';
eval(['cd ', top_dir])

% for removing electrodes
thr = 1.5;
releases = ['1', '2', '3'];

for r = 1%:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    eval(['cd ', release_dir '/protocols'])
    
    % remove parent and hidden directories, then get protocols
    folders = dir(pwd);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    
    for p = 1%:numel(protocols)
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
        for s = 1%:numel(subjects)
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
            for e = 1%:numel(experiments)
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1%:numel(sessions)
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
                    if numel(all_eegfiles) < numel(labels) % we only care if there are fewer files - having extra is not a problem and this happens frequently
                        % get channels of files you actually have
                        all_channels = cellfun(@(x) strsplit(x, '.'), all_eegfiles, 'UniformOutput', false);
                        all_channels = cellfun(@(x) str2double(x{end}), all_channels);
                        % match them to labels
                        idx = zeros(size(all_channels));
                        for k = 1:numel(idx)
                            curr = all_channels(k);
                            % add only matching files, skip the other
                            % ones
                            if any(chann_idx == curr)
                                idx(k) = find(chann_idx == curr);
                            end
                        end
                        % remove zeros from skipped channels
                        idx(idx == 0) = [];
                        % update
                        labels = labels(idx);
                        chann_idx = chann_idx(idx);
                        nChan = numel(labels);
                    end
                    
                    %% Get event info
                    
                    % load event info
                    if isfield(curr_info, 'all_events')
                        fid = fopen([release_dir, curr_info.all_events]);
                    else % some only have tasks events
                        fid = fopen([release_dir, curr_info.task_events]);
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
                    pretask = ~isempty(start_idx) && numel(1:(events(start_idx).eegoffset - 1)) > 1000;
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
                            warning('This session has multiple eegfiles')
                            
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
                        
                        eegfile = fields(val);
                        if numel(eegfile) == 1
                            eegfile = eegfile{1};
                        else
                            % check that there's only one file per session
                            warning('This session has multiple eegfiles in header')
                        end
                        eval(['header = val.', eegfile ';'])
                    end
                    
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
                    catch
                        % sometimes the size of the file and header.n_samples
                        % don't match. If that is the case, keep the file
                        % format listed, and just don't specify a size. also
                        % print a warning
                        warning('The number of samples indicated in the header and actual number of samples at the given precision do not match')
                        
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
                    
                    % check that we have enought post task data, now
                    % that we know how long the data is
                    if posttask
                        posttask = numel((events(end_idx).eegoffset + 1):header.n_samples) > 1000;
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
                    soz_ind = find(cell2mat(cellfun(@(x) contains(lower(x), 'seizure') | contains(lower(x), 'onset'), elec_cats, 'UniformOutput', false)));
                    interictal_ind = find(cell2mat(cellfun(@(x) contains(lower(x), 'interictal'), elec_cats, 'UniformOutput', false)));
                    bad_contact_ind = find(cell2mat(cellfun(@(x) ... % this one is complicated, some subjects have multiple categories of bad elecs
                        contains(lower(x), 'bad electrodes') | contains(lower(x), 'broken') | contains(lower(x), 'lesion'), elec_cats, 'UniformOutput', false)));
                    % check size - this is only a warning because some sites don't list SOZ,
                    % and some subjects don't have bad elecs marked...presumably because the
                    % recording was good. However this might also happen if something went
                    % wrong with the parsing
                    if (numel(soz_ind) ~= 1 || numel(interictal_ind) ~= 1 || numel(bad_contact_ind) < 1)
                        warning('Something went wrong with parsing the electrode categories, or this subject is missing categories.');
                    end
                    
                    % get contact names for each category
                    soz = elec_cats((soz_ind + 1):(interictal_ind - 1));
                    interictal_cont = elec_cats((interictal_ind + 1):(bad_contact_ind - 1));
                    bad_cont = elec_cats((bad_contact_ind(1) + 1):end);
                    % remove empty categories and category names for bad electrodes
                    bad_cont(cell2mat(cellfun(@(x) ...
                        contains(lower(x), 'none') | contains(lower(x), 'bad') | contains(lower(x), 'broken') | contains(lower(x), 'lesion'),...
                        bad_cont, 'UniformOutput', false))) = [];
                    
                    
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
                    
                    %% CAR
                    
                    for j = 1:nTrial
                        
                        ft_data.trial{j} = ft_data.trial{j} - mean(ft_data.trial{j},2); % demean
                        ft_data.trial{j} = detrend(ft_data.trial{j}')'; % as opposed to low pass filtering
                        ft_data.trial{j} = get_CAR(ft_data.trial{j}, ft_data.label); % CAR by group
                        
                        fprintf('Finished CAR!\n')
                        
                    end
                    
                end
            end
        end
    end
end

%%

[out,MARKER,envelope,background,discharges,envelope_pdf] = spike_detector_hilbert_v16_byISARG(ft_data.trial{1}', header.sample_rate);

out
all_soz = unique([cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), soz, 'UniformOutput', false)); cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), interictal_cont, 'UniformOutput', false))]);
num_soz = sum(sum(MARKER.M(:,cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), soz, 'UniformOutput', false)))))
num_ii = sum(sum(MARKER.M(:,cell2mat(cellfun(@(x) find(strcmpi(x, ft_data.label)), interictal_cont, 'UniformOutput', false)))))
numel(out.pos) - (num_soz + num_ii)


%% Try to weed out false positives
% find all spikes that occur within  120ms(?) of each other, and are
% present in at least XX of SOZ channels

%numbers = 2:7;
%wins = 0.05:0.05:.5;
numbers = 4;
wins = 0.2;
num_soz_params = zeros(numel(numbers), numel(wins));
num_other = zeros(numel(numbers), numel(wins));

for k = 1:numel(numbers)
    num_soz = numbers(k);%round(numel(soz)*0.33);
    for j = 1:numel(wins)
        win = ceil(wins(j)*MARKER.fs); % in samples, window size for finding spikes
        include_length = win;%.300*MARKER.fs;
        nSamp = size(MARKER.d,1);
        
        gold_spikes = zeros(size(MARKER.d));
        nWin = ceil((nSamp - win)./win);
        cnt = 0;
        for i = 1:win:(nSamp - win)
            cnt = cnt + 1;
            curr_win = MARKER.M(i:(i+win-1),:);
            curr_soz = curr_win(:,all_soz);
            if sum(any(curr_soz)) > num_soz
                gold_spikes(i:(i+win+include_length),:) = MARKER.M(i:(i+win+include_length),:);
            end
        end
        
        MARKER.M_clean = gold_spikes;
        sum(gold_spikes)
        num_soz_params(k,j) = sum(sum(MARKER.M_clean(:,all_soz)))
        num_other(k,j) = sum(sum(MARKER.M_clean)) - sum(sum(MARKER.M_clean(:,all_soz)))
    end
    
end

figure(1); clf
imagesc(wins, numbers, num_soz_params); colorbar
xlabel('Win Size')
ylabel('Num Chan')
title('SOZ')

figure(2); clf
imagesc(wins, numbers, num_other); colorbar
xlabel('Win Size')
ylabel('Num Chan')
title('Other')

%% Plot everything

eegplot(MARKER.d', 'srate', MARKER.fs)
eegplot(MARKER.M', 'srate', MARKER.fs)
eegplot(MARKER.M_clean', 'srate', MARKER.fs)

%% Plot

l = 2000;
for i = 1:numel(chann_idx)
    curr = i;
    
    st = find(MARKER.M(:,curr) ~= 0);
    if ~isempty(st)
        st = st(1);
        st = st - 100;
    else
        st = 1;
    end
    
    figure(4); clf
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    subplot(2,1,1)
    plot(1:numel(st:(st+l)), MARKER.d(st:(st+l),curr), 'linewidth', 1.5); hold on
    plot(find(MARKER.M(st:(st+l),curr) ~= 0), MARKER.d(MARKER.M(st:(st+l),curr) ~= 0,curr), 'rx', 'linewidth', 2)
    title('Voltage')
    xlim([1,l])
    subplot(2,1,2)
    plot(1:numel(st:(st+l)), envelope(st:(st+l),curr), 'linewidth', 1.5); hold on
    plot(find(MARKER.M(st:(st+l),curr) ~= 0), envelope(MARKER.M(st:(st+l),curr) ~= 0,curr), 'rx', 'linewidth', 2)
    title('Envelope')
    xlim([1,l])
    hold off
    pause(0.5)
    %saveas(gca, [top_dir, 'img/spike_examples/',num2str(i),'.png'], 'png')
end



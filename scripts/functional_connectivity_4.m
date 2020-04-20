%% Get functional connectivity

clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath('/Users/stiso/Documents/MATLAB/arfit/')

%% Run this on a remote machine (main FC calculation)


% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '_delphos';

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
            
    % just for testing that qsub fuction will work
    for p = 1:numel(protocols)
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
        for s = 1:numel(subjects)
            subj = subjects{s};
            
            functional_connectivity(protocol, release, top_dir, subj, detector, spike_win, win_length)
        end
    end
end

%% combine subjects with duplicates in release 1 and release 3

protocol = 'r1';
all_subjects = [];
for r = 1:numel(releases)
    rel = releases(r);
    
    release_dir = [top_dir, 'release', rel '/'];
    
    % for catching errors
    errors = struct('files', [], 'message', []);
    warnings = struct('files', [], 'message', []);
    
    % remove parent and hidden directories, then get protocols
    folders = dir([release_dir '/protocols']);
    folders = {folders([folders.isdir]).name};
    protocols = folders(cellfun(@(x) ~contains(x, '.'), folders));
    protocol = protocols{1}; % theres only 1 protocol;
    
    % get global info struct
    fname = [release_dir 'protocols/', protocol, '.json'];
    fid = fopen(fname);
    raw = fread(fid);
    str = char(raw');
    fclose(fid);
    info = jsondecode(str);
    eval(['info = info.protocols.', protocol,';']);
    
    % get subjects
    release_vect = cell(numel(fields(info.subjects)),1);
    [release_vect{:}] =  deal(r);
    all_subjects = [all_subjects; [fields(info.subjects), release_vect]];
end

% get duplicates
[U, I] = unique(all_subjects(:,1), 'first');
idx = 1:length(all_subjects(:,1));
idx(I) = [];
duplicates = all_subjects(idx,1);
demographics = [{'race'}, {'hand'}, {'gender'}, {'age'}];

% loop through duplicates and combine csv files
for i = 1:numel(duplicates)
    subj = duplicates{i};
    
    % find all the releases with this subject
    curr_rels = all_subjects(strcmpi(all_subjects(:,1), subj),2);
    
    % load base case
    top_idx = 1;
    flag = 0;
    top_rel = curr_rels{1};
    while top_idx < numel(curr_rels) && ~flag
        if exist([top_dir, 'FC/release',num2str(top_rel), '/', protocol, '/', subj, '/win_', num2str(win_length), '/fc_data.csv'],'file')
            fc_table = readtable([top_dir, 'FC/release',num2str(top_rel), '/', protocol, '/', subj, '/win_', num2str(win_length), '/fc_data.csv']);
            % for demographics checks later
            if fc_table.age(1) == 0
                fc_table.age = NaN(numel(fc_table.age), 1);
            end
            curr_rels = curr_rels((top_idx + 1):end);
            flag = 1;
            nObs = size(fc_table,1);
        end
        top_idx = top_idx + 1;
    end
    
    if exist('fc_table', 'var')
        % cycle through extras add thme onti the end, with the release appended
        %to session
        for r = 1:numel(curr_rels)
            if exist([top_dir, 'FC/release',num2str(curr_rels{r}), '/', protocol, '/', subj, '/win_', num2str(win_length), '/fc_data', detector, '.csv'], 'file')
                new_table = readtable([top_dir, 'FC/release',num2str(curr_rels{r}), '/', protocol, '/', subj, '/win_', num2str(win_length), '/fc_data', detector, '.csv']);
                if new_table.age(1) == 0
                    new_table.age = NaN(numel(new_table.age), 1);
                end
                % check that demgraphics line up
                for f = 1:numel(demographics)
                    demo = demographics{f};
                    try
                        dem_flag = new_table.(demo)(1) ~= fc_table.(demo)(1);
                    catch
                        dem_flag = ~strcmp(new_table.(demo)(1), fc_table.(demo)(1));
                    end
                    if dem_flag
                        % get unique indicator for pattern of nans
                        try
                            nan_new = isnan(new_table.(demo)(1));
                        catch
                            nan_new = cellfun(@isnan,new_table.(demo)(1));
                        end
                        try
                            nan_fc = isnan(fc_table.(demo)(1));
                        catch
                            nan_fc = cellfun(@isnan,fc_table.(demo)(1));
                        end
                        nan_flag = [nan_new, nan_fc]*[2;1];
                        % if one is nan, replace it with the one with a real
                        % value. otherwise set all to nan
                        if nan_flag == 2
                            new_table.(demo) = repmat(fc_table.(demo)(1), numel(new_table.(demo)), 1);
                        elseif nan_flag == 1
                            fc_table.(demo) = repmat(new_table.(demo)(1), numel(fc_table.(demo)), 1);
                        else
                            fc_table.(demo) = NaN(numel(fc_table.(demo)), 1);
                            new_table.(demo) = NaN(numel(new_table.(demo)), 1);
                        end
                    end
                end
                fc_table = [fc_table; new_table];
                
                %when a new file is added, delete it
                eval(['delete ', top_dir, 'FC/release',num2str(curr_rels{r}), '/', protocol, '/', subj, '/win_', num2str(win_length), '/fc_data', detector, '.csv']);
            end
        end
        %save new file in top directory
        if size(fc_table,1) > nObs
            writetable(fc_table, [top_dir, 'FC/release',num2str(top_rel), '/', protocol, '/', subj, '/win_', num2str(win_length), '/fc_data', detector, '.csv']);
        end
    end
end


%% get and concatenate errors and warnings

errors_all = struct('files', [], 'message', []);
for r = 1:numel(releases)
    release = releases(r);
    release_dir = [top_dir, 'release', release '/'];
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
            subj_dir = [top_dir, 'FC/release',release, '/', protocol, '/', subj, '/'];
            curr_err = load([release_dir, 'protocols/fc_errors.mat']);
            errors_all= [errors_all, curr_err.errors];
        end
    end
end
save([top_dir, 'fc_errors', detector, '.mat'], 'errors_all')


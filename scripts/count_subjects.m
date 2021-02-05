% global variables
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

releases = ['1', '2', '3'];

spike_win = 0.05; %for loading spike data
win_length = 1; % in seconds
detector = '_delphos';
cnt = 0
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
            
            % useful variables
            table_names = [{'subj'}, {'exper'}, {'sess'}, {'time'}, {'power'}, {'fc_measure'}, {'band'}, {'str'}, {'ti'},...
                {'str_soz'}, {'str_not_soz'}, {'str_spike'}, {'str_not_spike'}, {'elec'}, {'region'}, {'elec_in_soz'},...
                {'elec_has_spike'}, {'spike_num'}, {'spike_spread'}, {'age'}, {'gender'}, {'race'}, {'hand'}, {'x'}, {'y'}, {'z'}, {'type'}];
            %bands
            freqs = unique(round(logspace(log10(4),log10(150),30)));
            bands = [4, 8; 9, 15; 16 25; 36, 70; 71, 150];
            band_names = [{'theta'}, {'alpha'}, {'beta'}, {'gamma'}, {'hg'}];
            
            %fc measures
            measure_names = [{'coh'}, {'im_coh'}, {'plv'}, {'aec'}, {'aec_ortho'}, {'xcorr'}, {'ar'}, {'pac'}];
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
            
            fc_table = cell2table(cell(0,27), 'VariableNames', table_names);
            
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            
            %
            if exist([top_dir, 'FC/release',release, '/', protocol, '/', subj, '/', 'win_', num2str(win_length), '/fc_data', detector, '.csv'], 'file')
                fprintf('\n******************************************\nStarting functional connectivity for subject %s...\n', subj)
                cnt = cnt + 1
                subjects{cnt} = subj;
            end
        end
    end
    
    
end

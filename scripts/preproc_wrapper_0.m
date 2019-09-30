clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

%%

% global variables and packages
top_dir = '/Volumes/bassett-data/Jeni/RAM/';
eval(['cd ', top_dir])

% for removing electrodes
thr = 1.5;
releases = ['1', '2', '3'];

% for catching errors
errors = struct('files', [], 'message', []);
warnings = struct('files', [], 'message', []);

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
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};

                    try
                        warnings = preproc(thr, release_dir, info, metadata_dir, top_dir, release, protocol, subj, exper, sess, warnings);
                        
                     catch ME
                         errors(end+1).files = [subj, '_', exper, '_', sess];
                         errors(end).message = ME.message;
                     end
                end
            end
            diary off
        end
    end
end

% remove empty entry
errors = errors(2:end);
warnings = warnings(2:end);
save([top_dir, 'errors.mat'], 'errors');
save([top_dir, 'warnings.mat'], 'warnings');

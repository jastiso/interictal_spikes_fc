%% Remove datasets
% find data sets with a ton of artifacts, or really bad PSDs
% this will require some manual inspection of some diagnostic files saved
% in earlier steps

clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
top_dir = '/Volumes/bassett-data/Jeni/RAM/';

load([top_dir, 'n_artifacts.mat'])
releases = ['1', '2', '3'];

%% get subj with a lot of artifacts

problem_data = cell2mat(n_art(5,:)) > 1000;
sum(problem_data)

% go look at the these lfp plots
problem_art = n_art(:,problem_data);
problem_ext = cell(sum(problem_data),1);

for i = 1:numel(problem_ext)
    problem_ext{i} = [problem_art{1,i}, '_', problem_art{2,i}, ...
        '_', problem_art{3,i}, '_', num2str(problem_art{4,i})];
end

% for now I'm going to throw all of these out because I think I can afford
% to

%% Look at FFTs

% for FFT
win_len = 500;
max_freq = 200;
min_freq = 3;
sse = [];
names = {};
curves = zeros(1,max_freq*2-1);

for r = 1:numel(releases)
    release = releases(r);
    
    release_dir = [top_dir, 'release', release '/'];
    
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
            
            % save command window
            %clc
            if ~exist([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/'], 'dir')
                mkdir([top_dir, 'processed/release',release, '/', protocol, '/', subj, '/']);
            end
            %eval(['diary ', [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/log.txt']]);
            
            fprintf('\n******************************************\nStarting FFT fit for subject %s...\n', subj)
            
            % get experiements
            eval(['experiments = fields(info.subjects.' subj, '.experiments);'])
            for e = 1:numel(experiments)
                exper = experiments{e};
                
                % get seesions
                eval(['sessions = fields(info.subjects.' subj, '.experiments.', exper, '.sessions);'])
                for n = 1:numel(sessions)
                    %initialize
                    out_clean = [];
                    marker = [];
                    
                    sess = sessions{n};
                    sess = strsplit(sess, 'x');
                    sess = sess{end};
                    % get the path names for this session, loaded from a json file
                    eval(['curr_info = info.subjects.' subj, '.experiments.' exper, '.sessions.x', sess, ';'])
                    
                    % folders
                    save_dir = [top_dir, 'processed/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    img_dir = [top_dir, 'img/diagnostics/release',release, '/', protocol, '/', subj, '/', exper, '/', sess, '/'];
                    if ~exist(img_dir, 'dir')
                        mkdir(img_dir);
                    end
                    
                    % check if this subect has clean data
                    if exist([save_dir, 'data_clean.mat'], 'file')
                        load([save_dir, 'data_clean.mat'])
                        load([save_dir, 'header.mat'])
                        
                        nTrial = numel(ft_data.trial);
                        for j = 1:nTrial
                             % get avg fft
                             [pxx,f] = pwelch(ft_data.trial{j}',hanning(win_len),0, [], header.sample_rate);
                             pxx = mean(10*log10(pxx),2);
                             % select only freqs in bands of interest
                             f = f(f > min_freq & f < max_freq);
                             pxx = pxx(f > min_freq & f < max_freq);
                             % exclude notch filters
                             exc = (f < 10) | (f > 59 & f < 61) | (f > 119 & f < 121) | (f > 179 & f < 181); % filter width taken from preproc function
                             
                             % fit exponenetial
                             [curve, goodness] = fit(f,pxx, 'smoothingspline', 'SmoothingParam', 0.01, 'Exclude', exc);
                             
                             figure(1); clf
                             plot(curve,f,pxx);
                             title(num2str(goodness.sse/goodness.dfe));
                             %saveas(gca, [img_dir, 'polyfit_', num2str(j), '.png'], 'png')
                             
                             sse = [sse, goodness.sse/goodness.dfe];
                             curves(end+1,:) = curve(1:0.5:(max_freq));
                             names{end+1} = [subj, '_', exper, '_', sess, '_', num2str(j)];
                             fprintf('The MSE for this dataset was %d\n', goodness.sse/goodness.dfe)
                        end
                    end
                end
            end
        end
    end
end


%% Remove subjects

names = names(2:end);

% spikes - sse
sse_thr = quantile(sse, .9);

% shape of psd
curves = curves(2:end,:); % get rid of empty entry
psd_sim = mean(corr(curves', 'type', 'spearman'));
psd_thr = quantile(psd_sim, .1);

bad_datasets = names(psd_sim < psd_thr | sse > sse_thr | problem_data); 

save([top_dir, 'bad_datasets.mat'], 'bad_datasets')

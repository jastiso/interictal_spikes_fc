%% Plot BNV

clear
clc
close all
warning ON

addpath(genpath('/Users/stiso/Documents/Code/interictal_spikes_fc/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath('/Users/stiso/Documents/MATLAB/arfit/')

%% load 
% node data

top_dir = '/Volumes/bassett-data/Jeni/RAM/';

betas = readtable([top_dir, 'group_analysis/node_stats.csv']);
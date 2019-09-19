function [data_CAR] = get_CAR(data, elec_labels)
% Common average reference by electrode
%   Basically taked the average signal for each electrode and subtracts it
% data should be demeaned and detrended already. Bad electrodes should
% already be removed
% 
% Inputs:
% data              N x T (channels by time) matrix
% elec_labels       cell vector of electrode labels (clinical). The point
%                   is to group them by the same physical electrode, not the same region in
%                   MNI space

% Heuristic for creating groups is to match all the letters in the elec
% labels

% @author JStiso 07/19

% get groups
prefix = cellfun(@(x) x(isstrprop(x,'alpha')),elec_labels, 'UniformOutput', 0);
groups = unique(prefix);
nGroup = numel(groups);

fprintf('Starting CAR with %s groups...\n', num2str(nGroup))

data_CAR = zeros(size(data));
for g = 1:nGroup
    curr_group = groups{g};
    group_idx = cellfun(@(x) strcmp(x, curr_group), prefix);
    data_CAR(group_idx, :) = data(group_idx,:) - mean(data(group_idx,:));
end




end


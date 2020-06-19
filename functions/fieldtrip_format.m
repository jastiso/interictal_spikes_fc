function [ft_data] = fieldtrip_format(data, srate, elec_labels, trl)
%get data into fieldtrip format

ft_data = [];
ft_data.trial{1,1} = data;
ft_data.time{1,1} = (1/srate):(1/srate):(size(data,2)/srate);
ft_data.label = elec_labels;
ft_data.fsample = srate;
ft_data.nSamples = size(data,2);


if ~isempty(trl)
    ft_data.sampleinfo = [trl(1,1), trl(end,2)];
    cfg = [];
    cfg.trl = trl;
    ft_data = ft_redefinetrial(cfg, ft_data);
else
    ft_data.sampleinfo = [1, size(data,2)];
    cfg = [];
    ft_data = ft_preprocessing(cfg,ft_data);
end

end


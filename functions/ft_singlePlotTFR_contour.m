function ft_singlePlotTFR_contour( cfg,data )
%ft_singlePlotTFR_contour plots a single channel or ROI on a log scale with
%contour
%   use as: ft_singleplotTFR_contour(cfg,data)
%   in which cfg = channels -> cell array containing the channels to plot
%   data is. 
%  v17012016 - berryvdberg@gmail.com

chanIdx  = ismember(data.label,cfg.channels);

switch data.dimord
    case 'subj_chan_freq_time'
         data.powspctrm = squeeze(mean(mean(data.powspctrm(:,chanIdx,:,:),1),2));
    case 'chan_freq_time'
        data.powspctrm = squeeze(mean(data.powspctrm(chanIdx,:,:),1));
    otherwise
        error('dimord not supported')
end


% make the actual contour plot
contourf(data.time,1:numel(data.freq),data.powspctrm,cfg.ncontours);
yLabels = get(gca, 'YTickLabel');
if ~iscell(yLabels)
    yLabels = str2num(yLabels);
end


% find index and create new labels in a cell of strings
newYLabels = [];
for i=1:numel(yLabels)
    if  ~iscell(yLabels)
        ind = yLabels(i);
    else
        ind = str2num(yLabels{i});
    end
  newYLabels{i} = num2str(data.freq(ind));
end

set(gca, 'YTickLabel', newYLabels')


% set xlim if field exists
if isfield(cfg,'xlim')
    xlim(cfg.xlim)
end

% set clim if field exists
if isfield(cfg,'zlim')
    set(gca, 'Clim', cfg.zlim)
end


end
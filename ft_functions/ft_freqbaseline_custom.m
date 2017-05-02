function [freq] = ft_freqbaseline_custom(cfg, freq)

% FT_FREQBASELINE performs baseline normalization for time-frequency data
%
% Use as
%    [freq] = ft_freqbaseline(cfg, freq)
% where the freq data comes from FT_FREQANALYSIS and the configuration
% should contain
%   cfg.baseline     = [begin end] (default = 'no')
%   cfg.baselinetype = 'absolute', 'relchange', 'relative', or 'db' (default = 'absolute')
%   cfg.parameter    = field for which to apply baseline normalization, or
%                      cell array of strings to specify multiple fields to normalize
%                      (default = 'powspctrm')
%
% See also FT_FREQANALYSIS, FT_TIMELOCKBASELINE, FT_FREQCOMPARISON,
% FT_FREQGRANDAVERAGE

% Undocumented local options:
%   cfg.inputfile  = one can specifiy preanalysed saved data as input
%   cfg.outputfile = one can specify output as file to save to disk

% Copyright (C) 2004-2006, Marcel Bastiaansen
% Copyright (C) 2005-2006, Robert Oostenveld
% Copyright (C) 2011, Eelke Spaak
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_freqbaseline.m 9520 2014-05-14 09:33:28Z roboos $

revision = '$Id: ft_freqbaseline.m 9520 2014-05-14 09:33:28Z roboos $';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar freq

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% check if the input data is valid for this function
freq = ft_checkdata(freq, 'datatype', 'freq', 'feedback', 'yes');

% update configuration fieldnames
cfg              = ft_checkconfig(cfg, 'renamed', {'param', 'parameter'});

% set the defaults
cfg.baseline     =  ft_getopt(cfg, 'baseline', 'no');
cfg.baselinetype =  ft_getopt(cfg, 'baselinetype', 'absolute');
cfg.parameter    =  ft_getopt(cfg, 'parameter', 'powspctrm');

% check validity of input options
cfg =               ft_checkopt(cfg, 'baseline', {'char', 'doublevector'});
cfg =               ft_checkopt(cfg, 'baselinetype', 'char', {'absolute', ...
    'relative', 'relchange','db','db_trial','ztransformSession','logpower'});
cfg =               ft_checkopt(cfg, 'parameter', {'char', 'charcell'});

% make sure cfg.parameter is a cell array of strings
if (~isa(cfg.parameter, 'cell'))
  cfg.parameter = {cfg.parameter};
end

% is input consistent?
if ischar(cfg.baseline) && strcmp(cfg.baseline, 'no') && ~isempty(cfg.baselinetype)
  warning('no baseline correction done');
end

% process possible yes/no value of cfg.baseline
if ischar(cfg.baseline) && strcmp(cfg.baseline, 'yes')
  % default is to take the prestimulus interval
  cfg.baseline = [-inf 0];
elseif ischar(cfg.baseline) && strcmp(cfg.baseline, 'no')
  % nothing to do
  return
end

% check if the field of interest is present in the data
if (~all(isfield(freq, cfg.parameter)))
  error('cfg.parameter should be a string or cell array of strings referring to (a) field(s) in the freq input structure')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output structure
freqOut        = [];
freqOut.label  = freq.label;
freqOut.freq   = freq.freq;
freqOut.dimord = freq.dimord;
freqOut.time   = freq.time;

if isfield(freq, 'grad')
  freqOut.grad = freq.grad;
end
if isfield(freq, 'elec')
  freqOut.elec = freq.elec;
end
if isfield(freq, 'trialinfo')
  freqOut.trialinfo = freq.trialinfo;
end

% loop over all fields that should be normalized
% X=nanmean(nanmean(nanmean(freq.powspctrm,2),3),4);
% zX=zscore(X); rej=zX>4; 
% freqOut.rej=or(~cfg.inc_trials,rej); clear X zX;
freqOut.rej=~cfg.inc_trials; clear X zX;

for k = 1:numel(cfg.parameter)
  par = cfg.parameter{k};
  
  if strcmp(freq.dimord, 'chan_freq_time_X')
    
    freqOut.(par) = ...
      performNormalization(freq.time, freq.(par), cfg.baseline, cfg.baselinetype);
    
  elseif strcmp(freq.dimord, 'subj_chan_freq_time') || ...
          strcmp(freq.dimord, 'rpt_chan_freq_time') || ...
          strcmp(freq.dimord, 'chan_chan_freq_time')

    freqOut.(par) = zeros(size(freq.(par)));
    %=====================================================================%
    % Z transform custom code
    %=====================================================================%
    if (strcmp(cfg.baselinetype,'logpower'))
        
        selectedData = freq.(par);
        LselectedData=log10(selectedData);
        
        % mu=squeeze(nanmean(LselectedData(cfg.inc_trials,:,:,:),1));
        % sigma=squeeze(nanstd(LselectedData(cfg.inc_trials,:,:,:),1));
        for l = 1:size(freq.(par), 1)
            freqOut.(par)(l,:,:,:) = (squeeze(LselectedData(l,:,:,:))); %-mu)./sigma;
        end  
    elseif (strcmp(cfg.baselinetype,'ztransformSession'))
        
        selectedData = freq.(par);
        LselectedData=log10(selectedData);
        
        mu=squeeze(nanmean(LselectedData(cfg.ztrials,:,:,:),1));
        sigma=squeeze(nanstd(LselectedData(cfg.ztrials,:,:,:),1));
        for l = 1:size(freq.(par), 1)
            freqOut.(par)(l,:,:,:) = (squeeze(LselectedData(l,:,:,:))-mu)./sigma;
        end  
    elseif (strcmp(cfg.baselinetype,'CI'))
        % Compute average baseline conversion purposes (of entire set!)
        baselineTimes = (freq.time >= cfg.baseline(1) & freq.time <= cfg.baseline(2));
        data=squeeze(nanmean(freq.(par)(cfg.inc_trials,:,:,:)));
%         meanVals = repmat(nanmean(data(:,:,baselineTimes), 3), [1 1 size(data, 3)]);
        clear data;
        for l = 1:size(freq.(par), 1)
            tfdata = freq.(par)(l,:,:,:);
            siz    = size(tfdata);
            tfdata = reshape(tfdata, siz(2:end));
            % Normalize based upon single trials (off)
            meanVals = repmat(nanmean(tfdata(:,:,baselineTimes), 3), [1 1 size(tfdata, 3)]);
            freqOut.(par)(l,:,:,:) = ...
                performNormalization(tfdata,cfg.baselinetype,meanVals);
        end % Trial by trial loop
    else
        % Compute average baseline conversion purposes (of entire set!)
        baselineTimes = (freq.time >= cfg.baseline(1) & freq.time <= cfg.baseline(2));
        data=squeeze(nanmean(freq.(par)(cfg.inc_trials,:,:,:)));
%         meanVals = repmat(nanmean(data(:,:,baselineTimes), 3), [1 1 size(data, 3)]);
        clear data;
        for l = 1:size(freq.(par), 1)
            tfdata = freq.(par)(l,:,:,:);
            siz    = size(tfdata);
            tfdata = reshape(tfdata, siz(2:end));
            % Normalize based upon single trials (off)
            meanVals = repmat(nanmean(tfdata(:,:,baselineTimes), 3), [1 1 size(tfdata, 3)]);
            freqOut.(par)(l,:,:,:) = ...
                performNormalization(tfdata,cfg.baselinetype,meanVals);
        end % Trial by trial loop
    end      
  else
    error('unsupported data dimensions: %s', freq.dimord);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output scaffolding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous freq

% rename the output variable to accomodate the savevar postamble
freq = freqOut;

ft_postamble history freq
ft_postamble savevar freq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that actually performs the normalization on an arbitrary quantity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = performNormalization(data, baselinetype, meanVals)

if length(size(data)) ~= 3,
  error('time-frequency matrix should have three dimensions (chan,freq,time)');
end

if (strcmp(baselinetype, 'absolute'))
  data = data - meanVals;
elseif (strcmp(baselinetype, 'relative'))
  data = data ./ meanVals;
elseif (strcmp(baselinetype, 'relchange'))
  data = (data - meanVals) ./ meanVals;
elseif (strcmp(baselinetype, 'db'))
  data = 10*log10(data ./ meanVals);
else
  error('unsupported method for baseline normalization: %s', baselinetype);
end
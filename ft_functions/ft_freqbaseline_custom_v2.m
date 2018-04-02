function [freq_out] = ft_freqbaseline_custom_v2(cfg, freq)
% FT_FREQBASELINE performs baseline normalization for time-frequency data
%
% BRG NOTES: This scrip has been re-written Z-transform and log transform 
% the data only! It expects freq to have .powspctrm with the dimensions 
% being {trial X channel X freq X time}
%
% cfg.inc_trials -> include these trials
% cfg.ztrials    -> Z-score mean/mu values (subset of inc_trials)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq_out=freq;

selectedData = freq.powspctrm;       % Rename selected data
LselectedData=log10(selectedData);   % Log10 transform
%=====================================================================%
% Z transform custom code
%=====================================================================%
if (strcmp(cfg.baselinetype,'logpower'))
    freq_out.powspctrm=LselectedData;
elseif (strcmp(cfg.baselinetype,'ztransformSession'))
    OUT = zeros(size(freq.powspctrm));   % Initialize output
    mu=squeeze(nanmean(LselectedData(cfg.ztrials,:,:,:),1));
    sigma=squeeze(nanstd(LselectedData(cfg.ztrials,:,:,:),1));
    for l = 1:size(freq.powspctrm, 1) % Number of trials
        OUT(l,:,:,:) = (squeeze(LselectedData(l,:,:,:))-mu)./sigma;
    end  
    freq_out.powspctrm=OUT;
else
    error('unsupported transformation function');
end


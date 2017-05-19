% B.R. Geib (Winter 2015)
% Function file
%
% ft_plot(ctrast,ID,scale,save_as,interactive)
%
% Inputs:
%   ctrast{}        -> cell array of contrast to plot, up to a maximum of
%                      four. Plotted in order: 1=blue, 2=red, 3=green,
%                      4=black. Keep this at 2 for topo-maps. Will plot
%                      ctrast{1} - ctrast{2}
%   ID              -> Subject to plot (0 = all subjects)
%   scale           -> scale of ERP maps
%   save_as         -> save name, by default all images are saved to
%                      RUN.dir.QAL
%   interactive     -> interactive contrast maps, if turned off, topos are
%                      made and saved. by defualt, topos are set to be
%                      between [-2.5uV and 2.5uV]
% Outputs:
%   None
% Description: Plots the difference between two conditions and two time
% points. Also prefroms a t-test on the difference and displays weather the
% result is significant or not.


function Ga=ft_plot_v3(ctrast,cfg)

%=========================================================================%
% Assign Variables
%=========================================================================%
Ga=cfg.ga; 
cfg.subj=logical(cfg.subj); % Ensure cast correctly

% Can't recall why this template was setup
% cfg_erp2.layout = template{1}; 
% cfg_erp2.showlabels  = 'yes'; 
% if ~isempty(bs), cfg_erp2.baseline = bs; end
% cfg_erp2.ylim=cfg.scale;
% cfg_erp2.linewidth=2;
% cfg_erp2.xlim=[-0.3 0.3];
cfg_filt.lpfilter = 'yes';
cfg_filt.lpfreq = 10;
%=========================================================================%
% Generate Contrast
%=========================================================================%
for ii=1:length(ctrast)
    f1=ctrast{ii}; 
    if ~isfield(Ga.(f1),'avg'), Ga.(f1).avg=squeeze(mean(Ga.(f1).individual(cfg.subj,:,:),1)); end
    if ~isfield(Ga.(f1),'dimord'),Ga.(f1).dimord='chan_time'; end
    plot{ii}=Ga.(f1);
    plot{ii}.individual=Ga.(f1).individual(cfg.subj,:,:);
    if cfg.filt==1   
        plot{ii}=ft_preprocessing(cfg_filt,plot{ii}); 
        try
            plot{ii}.individual=plot{ii}.trial;
        catch err
            % If only one subject is selected
            plot{ii}.individual=plot{ii}.avg;
        end
    end
end    

%=========================================================================%
% Plot
%=========================================================================%
cfg_erp = [];
cfg_erp.showlabels  = 'yes'; 
% cfg_erp.baseline = RUN.pre.baseline;
cfg_erp.ylim=cfg.scale;
cfg_erp.linewidth=2;
cfg_erp.parameter='avg';
cfg_erp.layout = cfg.template;

figure(1); 
% set(gcf,'position',[0 0 1280 1024]);
switch length(ctrast)
    case 1, ft_multiplotER(cfg_erp,plot{1});  
    case 2, ft_multiplotER(cfg_erp,plot{1},plot{2});  
    case 3, ft_multiplotER(cfg_erp,plot{1},plot{2},plot{3});  
    case 4, ft_multiplotER(cfg_erp,plot{1},plot{2},plot{3},plot{4});  
end







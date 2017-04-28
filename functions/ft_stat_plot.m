function ft_stat_plot(ga,f1,f2,template)
global RUN;

f1='stm_p_1_R1hit'; f2='stm_p_1_R1miss';
A=ga.(f1);
B=ga.(f2);

% Discover Neighbors
cfg=[];
cfg.method='distance';
cfg.template=RUN.template.plt; % Actual distance
cfg.elec=template{2}.elec;
cfg.feedback='yes';
cfg.neighbourdist=0.6; % 0.4dm is default
neighbours=ft_prepare_neighbours(cfg);

cfg = [];
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = [-.2 1];      % [0 1]s
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 1000;

Nsub = size(A.individual,1);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_timelockstatistics(cfg,A,B);


cfg = [];
cfg.highlightsymbolseries = ['*','*','.','.','.'];
cfg.layout = template{2};
cfg.contournum = 0;
cfg.markersymbol = '.';
cfg.alpha = 0.05;
cfg.parameter='stat';
cfg.zlim = [-5 5];
ft_clusterplot(cfg,stat);
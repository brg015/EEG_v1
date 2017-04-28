function fq_plot(ga,c1,c2,template,bs)

cfg = []; cfg.layout = template; cfg.showlabels  = 'yes'; cfg.baseline = bs;
plot1=ga.(c1); plot1.powspctrm=plot1.powspctrm-ga.(c2).powspctrm;
ft_multiplotTFR(cfg,plot1);

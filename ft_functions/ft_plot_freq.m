function ft_plot_freq(f1,f2,ga_freq,cfg)

cfg_erp = [];
cfg_erp.showlabels  = 'yes'; 
cfg_erp.linewidth=1;
cfg_erp.layout = cfg.template;
cfg_erp.zlim=cfg.topo_lim;
cfg_erp.param='avg';

plot1 = ga_freq.(f1);
plot2 = ga_freq.(f2);

if ~strcmp(f1,f2)
    plot3=plot1;
    plot3.powspctrm=plot1.powspctrm-plot2.powspctrm;
    plot3.avg=squeeze(mean(plot3.powspctrm));
%     [~,~,~,a]=ttest(plot1.powspctrm-plot2.powspctrm);
%     plot3.avg=a.tstat;
%     plot3.avg(and(plot3.avg>-2,plot3.avg<2))=NaN;
%     plot3.powspctrm=plot3.avg;
else
    plot3=plot1;   
end

if ~strcmp(f1,f2)
    ft_multiplotTFR(cfg_erp,plot3); colorbar;
else
   ft_multiplotTFR(cfg_erp,plot1); colorbar;
end



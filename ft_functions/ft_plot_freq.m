function ft_plot_freq(f1,f2,ga_freq)
global RUN;

template=RUN.template2;

fpre='';

if ~strcmp(f1,f2), topo_lim=[-.1 .1]; else topo_lim=[-1 1]; end


cfg_topo = []; cfg_erp = [];
cfg_topo.commentpos='title';
cfg_topo.comment='no';
cfg_topo.layout = template{2}; 

cfg_erp.showlabels  = 'yes'; 
cfg_erp.linewidth=1;

% if ~interactive
    cfg_erp.layout = template{1}; 
    cfg_erp.layout.fscale=[2 2 1.5 1.5];
    cfg_erp.layout.scale=[2 2];
    cfg_erp.layout.scale_offset=[0 0.13];
    cfg_erp.layout.cfg.skipscale='yes';
    cfg_erp.layout.cfg.skipcomnt='yes';
    cfg_erp.layout.pos(66,:)=[-.8 -.8];
    cfg_erp.layout.pos(67,:)=[.8 -.8];
% else
%     cfg_erp.layout = template{2};
% end
% cfg_erp.ylim=[4 40];
% cfg_erp.zlim=[-.5 .5];
%-------------------------------------------------------------------------%



% save1=fullfile(save_dir,[fs '_' fpre '_both.png']);
%-------------------------------------------------------------------------%
plot1 = ga_freq.(f1);
plot2 = ga_freq.(f2);

cfg_erp.zlim=topo_lim;
cfg_erp.param='avg';

if ~strcmp(f1,f2)
    plot3=plot1;
    plot3.powspctrm=plot1.powspctrm-plot2.powspctrm;
    plot3.avg=squeeze(mean(plot3.powspctrm));
%     [~,~,~,a]=ttest(plot1.powspctrm-plot2.powspctrm);
%     plot3.avg=a.tstat;
%     plot3.avg(and(plot3.avg>-1,plot3.avg<1))=0;
%     plot3.powspctrm=plot3.avg;
else
    plot3=plot1;   
end

if ~strcmp(f1,f2)
    ft_multiplotTFR(cfg_erp,plot3); colorbar;
else
   ft_multiplotTFR(cfg_erp,plot1); colorbar;
end

set(gcf,'position',[0 0 1280 1024]);


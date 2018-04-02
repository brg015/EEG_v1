function ft_node_plot(cfg,T,p)

factor=cfg.factor;
for ii=1:length(T)
    if T(ii)>0, col='r.'; else col='b.'; end
    ms=ceil(abs(T(ii))*factor);
    plot(cfg.template(ii,1),cfg.template(ii,2),col,'markersize',ms); hold on;
    if ~isempty(p)
        if p(ii)<0.05
            plot(cfg.template(ii,1),cfg.template(ii,2),'k+','markersize',8); hold on;
        end
    end
end
set(gca,'xlim',[-1.2 1.2]); set(gca,'ylim',[-.8 .8]);
% subplot(1,2,2); hist(T);
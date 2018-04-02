% BRG 2017 Summer
% cfg.template(x,y)    -> electrode position
% cfg.template_label{} -> electrode names
% cfg.Ci(mod)          -> module assignment
% cfg.nxn(elec X elec) -> connectivity map
% cfg.factor           -> for dynamic line width
function ft_net_plot(cfg)
%-------------------------------------------------------------------------%
% Simple setup
%-------------------------------------------------------------------------%
I=find(cfg.network.time>=cfg.time(1) & cfg.network.time<=cfg.time(2));
time=cfg.network.time(I);
c=1; for jj=I
    a=squeeze(cfg.network.(cfg.contrast{1})(:,jj,:,:));
    b=squeeze(cfg.network.(cfg.contrast{2})(:,jj,:,:));
    [~,~,~,stat]=ttest(a-b); 
    RT(c,:,:)=squeeze(stat.tstat); c=c+1;
end
cfg.nxn=RT;

if ~isfield(cfg,'thresh'), cfg.thresh=0.025; end
% Unless of crazy small numbers this works
if ~isfield(cfg,'factor'), cfg.factor=1e-6; end
% Just set Ci == 1 for all if it don't exist
if ~isfield(cfg,'Ci'), cfg.Ci=ones(size(cfg.nxn,1),size(cfg.nxn,2)); end
%-------------------------------------------------------------------------%
% Plot Network
%-------------------------------------------------------------------------%
figure; set(gcf,'color','w');
Cset={'k.','r.','g.','b.','c.','m.','y.'}; % easy start

c=1; for tt=I
    
% Plot community
subplot(ceil(sqrt(length(time))),ceil(sqrt(length(time))),c);
for jj=1:size(unique(cfg.Ci(c,:)))
    plot(cfg.template(cfg.Ci(c,:)==jj,1),cfg.template(cfg.Ci(c,:)==jj,2),Cset{jj},'markersize',5);
    hold on;
end
set(gca,'xlim',[-1 1]); set(gca,'ylim',[-1 1]); axis square;
title(num2str(time(c)));

NXN=squeeze(cfg.nxn(c,:,:)); NXN(logical(eye(size(NXN))))=0;
NXN(logical(triu(NXN)))=0;
% Take top X
NXN_lin=sort(reshape(abs(NXN(logical(tril(NXN)))),1,[]));

% LIM=NXN_lin(floor(length(NXN_lin)*(1-cfg.thresh))); 
LIM=cfg.thresh;

NXN(abs(NXN)<LIM)=0;
for ii=1:length(NXN)
    for jj=1:length(NXN)
        if (jj<ii && abs(NXN(ii,jj))>0)
            % Linewdith is ceil of coh*10
            POS=cfg.template([ii,jj],:);
            if NXN(ii,jj)<0
                plot(POS(:,1),POS(:,2),'b','linewidth',ceil(abs(NXN(ii,jj))*cfg.factor)); hold on;
            else
                plot(POS(:,1),POS(:,2),'r','linewidth',ceil(abs(NXN(ii,jj)*cfg.factor))); hold on;
            end
        end
    end
end

for jj=1:size(unique(cfg.Ci(c,:)))
    plot(cfg.template(cfg.Ci(c,:)==jj,1),cfg.template(cfg.Ci(c,:)==jj,2),Cset{jj},'markersize',5);
    hold on;
end
set(gca,'xlim',[-1 1]); set(gca,'ylim',[-1 1]); axis square;
title(num2str(time(c)));
c=c+1;

end

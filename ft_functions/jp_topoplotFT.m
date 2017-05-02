function h = jp_topoplotFT(cfg,varargin)

% JP_TOPOPLOT_FT plots topographic maps using linear griddata interpolation.
%
% use as JP_TOPOPLOT(cfg,data,....)
%
% cfg.latency = -1:cfg.latencyInc:1.6;cfg.latency = [cfg.latency' cfg.latency'+cfg.latencyInc];
% cfg.foi = [8 12];
% cfg.perspective = {'top','back'};
% cfg.clim = [-0.07 0.07];
% cfg.layout = RUN.layoutmw64.elec;
%
% Joonkoo Park (joonkoo@umass.edu) 4-30-2015
% addapted for fieldtrip and added more functionality
% berry van den berg (berryv.dberg@gmail.com)
% default values when no options are given
%todo:
% add left,right
%-------------------------------------------------------------------------%
% BG edits
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
cfg.funcname = mfilename;
if nargin > 1
    cfg.dataname = {inputname(2)};
    for k = 3:nargin
        cfg.dataname{end+1} = inputname(k);
    end
end

Ndata = numel(varargin);

% Granularity
grn = 0.01;    % granularity of line boundaries

if isfield(cfg,'granularity')
    meshgrn = cfg.granularity;    % too fine granularity will slow down the plotting
else
    meshgrn = 0.05;
end
% number of topo maps (plus an additional one for colorbar)
numtopo = size(cfg.latency,1);
% number of bins
numbin  = Ndata; % only available for 1 bin
h = figure('Position',[100 100 200*numtopo 250*numbin],'color','w');

for ibin = 1 : numbin
    %% get the xyz coordinates here
    xyz = zeros(length(varargin{ibin}.label),3);
    % Minus one to avoid plotting right mastoid
    for i = 1:length(varargin{ibin}.label)-1;
        xyz(i,:) = cfg.layout.pnt(strcmp(varargin{ibin}.label{i},cfg.layout.label),:);
    end
    
    % Extract coordinate information from ERP struct and normalize
    chanlabel = varargin{ibin}.label;
    xyz = bsxfun(@rdivide,xyz,sqrt(sum(xyz.^2,2)));
    %%
    for itopo = 1 : numtopo
        subplot(numbin,numtopo,(ibin-1)*numtopo+itopo);
        idxtime = varargin{ibin}.time >= cfg.latency(itopo,1) & varargin{ibin}.time <= cfg.latency(itopo,2);
        
        if strcmp(varargin{ibin}.dimord, 'subj_chan_time');
            v = mean(mean(varargin{ibin}.individual(:,:,idxtime),1),3);
        elseif strcmp(varargin{ibin}.dimord, 'subj_chan_freq_time');
            idxfreq = varargin{ibin}.freq >= cfg.foi(1) & varargin{ibin}.freq <= cfg.foi(2);
            v = mean(mean(mean(varargin{ibin}.powspctrm(:,:,idxfreq,idxtime),1),3),4);
        elseif strcmp(varargin{ibin}.dimord, 'chan_time');
            v =mean(varargin{ibin}.individual(:,idxtime),2);
            
        else
            error('no support for dimord')
        end
        
        
        %%%%%%%%%%%%%%from here on will be the different perspectives%%%%%%%%%%%%%
        switch cfg.perspective{ibin}
            case 'back'
                % Back view is where X < 0
                chanidx = find(xyz(:,1)<=0.01);
                xyz_b = xyz(chanidx,:);
                chanlabelb = chanlabel(chanidx);
                x = -xyz_b(:,2); % *-1 is because you are looking from the back
                y = xyz_b(:,3);
                
                [xq,yq] = meshgrid(-1:meshgrn:1, -1:meshgrn:1);
                
                vq = griddata(x,y,v(chanidx),xq,yq,'v4');
                
                % trim outside the unit circle and the minimum y val
                % vq(yq < min(y)) = NaN;
                vq((xq.^2 + yq.^2) > 1) = NaN;
                vq((.15 .* xq.^2 + yq + .6) < 0) = NaN;
                
%                 contourf(xq,yq,vq,cfg.ncontours,'fill','on','LineStyle','-');
%                 surfc(xq,yq,vq); view(0,90); shading interp; hold on;
                
                pcolor(xq,yq,vq); shading interp; hold on;
                contourf(xq,yq,vq,cfg.ncontours,'fill','off','LineStyle','-');

                if isempty(cfg.clim)
                    maxval = ceil(max(abs(v)) * 10) / 10;
                    caxis([-maxval maxval]);
                else
                    caxis(cfg.clim(:)); % add support for different clim per bin?
                end
%                 %colorbar;
%                 hold on;
%                 % scatter(-xyz_b(:,2),xyz_b(:,3),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
%                 scatter(-xyz_b(:,2),xyz_b(:,3),12,'o','MarkerFaceColor','k','MarkerEdgeColor','none');
                
                % head & neck (hn_)
                hn_k = 0.45;     % head-neck joint
                hn_b = 1.15;     % mid-point of neck
                hn_x1 = -1:grn:hn_k;
                hn_x2 = hn_k:grn:1.2;
                hn_y1 = sqrt(1-hn_x1.^2);
                hn_a = (-hn_k / sqrt(1-hn_k.^2)) / (2*(hn_k-hn_b));     % curve of parabola
                hn_c = sqrt(1-hn_k.^2) - hn_a .* (hn_k - hn_b) .^ 2;
                hn_y2 = hn_a * (hn_x2 - hn_b) .^ 2 + hn_c;
                
                plot([hn_y1, hn_y2],-[hn_x1, hn_x2],'Color',[.2 .2 .2],'LineWidth',1.2/numtopo);
                plot(-[hn_y1, hn_y2],-[hn_x1, hn_x2],'Color',[.2 .2 .2],'LineWidth',1.2/numtopo);
                axis equal;
                % plot the channel labels
                if isfield(cfg,'plotchannellab')
                    switch cfg.plotchannellab
                        case 'yes'  
                            for i = 1 : length(chanlabelb)
                                text(-xyz_b(i,2)+.01,xyz_b(i,3)+.075,chanlabelb{i});
                            end
                            
                        otherwise
                    end
                end
                
                % little hack to remove some of the white space
                hold off;
                
                axis square off;
                set(gca,'XLim',[-1.2 1.2],'YLim',[-1.2 1.2]);
                
                sub_pos = get(gca,'position');
                
                if numtopo>1
                    set(gca,'Position',[sub_pos(1)-0.1 sub_pos(2)-0 sub_pos(3)*1.4 sub_pos(4)*1.4])
                end
                
                text(0,-1.4,sprintf('%.2f to %.2f',cfg.latency(itopo,:)),'HorizontalAlignment','center');
                
                %colorbar;
                hold on;
                % scatter(-xyz_b(:,2),xyz_b(:,3),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
                scatter(-xyz_b(:,2),xyz_b(:,3),8,'o','MarkerFaceColor','k','MarkerEdgeColor','none');
            case 'top'
                % Back view is where z > 0
                chanidx = find(xyz(:,3)>=-0.01);
                xyz_b = xyz(chanidx,:);
                
                
                x = -xyz_b(:,2);
                y = xyz_b(:,1);
                [xq,yq] = meshgrid(-1:meshgrn:1, -1:meshgrn:1);
                vq = griddata(x,y,v(chanidx),xq,yq,'v4');
                
                % trim outside the unit circle and the minimum y val
                vq((xq.^2 + yq.^2) > 1) = NaN;
                
                % trim outside the unit circle and the minimum y val
                % vq(yq < min(y)) = NaN;
                vq((xq.^2 + yq.^2) > 1) = NaN;
                
                contourf(xq,yq,vq,cfg.ncontours,'fill','on','LineStyle','-');
                if isempty(cfg.clim)
                    maxval = ceil(max(abs(v)) * 10) / 10;
                    caxis([-maxval maxval]);
                else
                    caxis(cfg.clim(:)); % add support for different clim per bin?
                end
                %colorbar;
                hold on;
%                 scatter(-xyz_b(:,2),xyz_b(:,3),'filled','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none');
                scatter(-xyz_b(:,2),xyz_b(:,1),2.5,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
                
                % head
                hn_x = -1:.01:1;
                hn_y = sqrt(1-hn_x.^2);
                
                % nose
                ns_x1 = -.14:grn:0;
                ns_x2 = 0:grn:.14;
                ns_y1 = 1.5*ns_x1 + 1.2;
                ns_y2 = -1.5*ns_x2 + 1.2;
                
                plot(hn_x, hn_y,'Color',[.2 .2 .2],'LineWidth',1.2/numtopo);
                plot(hn_x,-hn_y,'Color',[.2 .2 .2],'LineWidth',1.2/numtopo);
                plot(ns_x1, ns_y1,'Color',[.2 .2 .2],'LineWidth',1.2/numtopo);
                plot(ns_x2, ns_y2,'Color',[.2 .2 .2],'LineWidth',1.2/numtopo);
                axis equal;
                
                hold off;
                set(gca,'XLim',[-1.2 1.2],'YLim',[-1.2 1.2]);
                axis square off;
                sub_pos = get(gca,'position');
                
                % this helps to utilize more of the  space for each plot
                set(gca,'Position',[sub_pos(1)-0.1 sub_pos(2)-0 sub_pos(3)*1.4 sub_pos(4)*1.4])
                text(0,-1.2,sprintf('%.2f to %.2f',cfg.latency(itopo,:)),'HorizontalAlignment','center');
                
        end
    end
end
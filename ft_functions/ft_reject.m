% Remove lateral eyemovements

function [lrej,lv,trej,tv,Cdat] = ft_reject(data,cfg)
% sr = 500ms

% distance is based upon # points
% cfg.sstep             -> step distance (half)
% cfg.smax              -> step max
% cfg.tmax              -> jump max
% cgf.win(min max)      -> window to scan through
% cfg.type              -> 'lat' or 'jump'
% ----------------------------------------------------------------------- %
% Setup Filter
% ----------------------------------------------------------------------- %
sstep=cfg.sstep; % must be an even number
smax=cfg.smax;
tmax=cfg.tmax;

% Identify time window to scan in
twin=cfg.win;
I=find(data.time{1}>=twin(1) & data.time{1}<=twin(2));
% ----------------------------------------------------------------------- %
sfunc=[zeros(1,sstep/2) ones(1,sstep/2)]; % Step function
tfunc=[zeros(1,6) ones(1,6)]; % 6*2*4 48ms step function
% Loop through all trials
for ii=1:length(data.trial)
    x1=data.trial{ii}(21,I);
    x2=data.trial{ii}(22,I);
    x12=x1-x2; % Lateral frontal electrode difference
    
    C=abs(conv(sfunc,x12))/(sstep/2);
    Cdat(ii,:)=C;
    % *2 camptures the 
    lv(ii)=max(C);
    if max(C)>smax
        lrej(ii)=1;
    else
        lrej(ii)=0;
    end
    
    % Arupt data jumps
    for jj=1:64
        t=abs(diff(data.trial{ii}(jj,I)));
        if max(t)>tmax
            trej(ii,jj)=1;
        else
            trej(ii,jj)=0;
        end
        tv(ii,jj)=max(t);
    end 
end

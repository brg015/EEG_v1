% Remove lateral eyemovements

function [lrej,lv,trej,tv,Cdat] = ft_reject(data,cfg)

% Generate step filter
% sr = 250 Hz
% 50 pt step (200ms)
% 30 uv step (30uv)
sstep=cfg.sstep;
smax=cfg.smax;
tmax=cfg.tmax;

twin=cfg.win;
I=find(data.time{1}>=twin(1) & data.time{1}<=twin(2));


sfunc=[ones(1,sstep).*-1 ones(1,sstep)];

% Loop through all trials
for ii=1:length(data.trial)
    x1=data.trial{ii}(21,I);
    x2=data.trial{ii}(22,I);
    x12=x1-x2; % Lateral frontal electrode difference
    
    C=abs(conv(sfunc,x12)./50);
    Cdat(ii,:)=C;
    lv(ii)=max(C);
    if max(C)>smax
        lrej(ii)=1;
    else
        lrej(ii)=0;
    end
    
    for jj=1:64
        t=abs(diff(data.trial{ii}(jj,:)));
        if max(t)>tmax
            trej(ii,jj)=1;
        else
            trej(ii,jj)=0;
        end
        tv(ii,jj)=max(t);
    end 
end

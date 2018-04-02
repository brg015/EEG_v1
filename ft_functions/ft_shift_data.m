function data=ft_shift_data(cfg,data)
% BRG function file (Fall 2017)
% 
% cfg.RT         -> RTs to shift timecourse by
% cfg.rej        -> Rejected trials
% cfg.tlim       -> 
%-------------------------------------------------------------------------%
% Generate null temporal trial
T=(cfg.tlim(1):1/data.fsample:cfg.tlim(end));
rejRT=false(1,length(cfg.RT));               % Rejected RTs
rejRT(cfg.RT>cfg.RTmax);
for itrial=1:length(data.trial)
   % Pull out proper temporal vector
   if cfg.RT(itrial)<cfg.RTmax, RT=cfg.RT(itrial); else RT=0; end
   % Anchor RT to the closest timepoint to keep vectors = length
   t1=data.time{itrial}>=(RT+cfg.tlim(1));
   t2=data.time{itrial}<=(RT+cfg.tlim(2));
   t=(t1 & t2);
   if sum(t)~=length(T)
       % Indicates that vector is off by a tiny bit, gotta be padded at the
       % beginning or end (4ms)
       I=find(t==1); t(I(end)+1)=1;
       if sum(t)~=length(T)
           I=find(t==1); t(I(1)-1)=1;
       end
   end
   
   if sum(t)~=length(T)
       keyboard;
   end
   % Pull the original baseline [-100 0ms]
   t11=data.time{itrial}>=-.1;
   t12=data.time{itrial}<=0;
   t13=(t11 & t12);
   uV=data.trial{itrial}(:,t13);
   clear t11 t12 t13;
   % Now pull out the actual data
   trial=data.trial{itrial}(:,t);
   data.trial{itrial}=trial-mean(uV,2);
   data.time{itrial}=T;
   clear t1 t2 t trial;
end
data.reject=(data.reject | rejRT);
data.RT_reject=rejRT;
%-------------------------------------------------------------------------%

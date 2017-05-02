function ITC = ft_ITC(booleans,power)

fn=fieldnames(booleans);
% Unpackage power into ITC structure
ITC.label=power.label;
ITC.freq=power.freq;
ITC.time=power.time;
ITC.elec=power.elec;
ITC.cfg=power.cfg;
for ii=1:length(fn)
    % Logic from http://www.fieldtriptoolbox.org/faq/itc
    F=power.fourierspctrm(logical(booleans.(fn{ii})),:,:,:); % data
    N=size(F,1); % number of trials
    ITC.itpc.(fn{ii})=F./abs(F); % divide by amplitude
    ITC.itpc.(fn{ii})=sum(ITC.itpc.(fn{ii}),1); % sum angles
    ITC.itpc.(fn{ii})=abs(ITC.itpc.(fn{ii}))/N; 
    ITC.itpc.(fn{ii})=squeeze(ITC.itpc.(fn{ii}));
end
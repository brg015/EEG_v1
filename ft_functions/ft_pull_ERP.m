function out=ft_pull_ERP(data,elec,time,rdata)
out=[];

if isempty(rdata)
    for ii=1:length(data.trial), rdata(ii,:,:)=data.trial{ii}; end
end

rtime=data.time{1};
relec=data.elec.label;

% I1 -> electrode index
% I2 -> time index
[~,~,I1]=intersect(elec,relec); 
I2=find((rtime>=time(1) & rtime<=time(2)));

% Pull wanted values
if length(I1)>1
    t1=rdata(:,I1,I2);
    out=squeeze(mean(mean(t1,2),3)); clear t1;
else
    t1=rdata(:,I1,I2);
    out=squeeze(mean(mean(t1,2),3)); clear t1;
end



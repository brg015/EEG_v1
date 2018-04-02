function out=ft_pull_FRQ(rdata,elec,time,freq)

rtime=rdata.freqlock.time;
relec=rdata.freqlock.label;
rfreq=rdata.freqlock.freq;

% I1 -> electrode index
% I2 -> time index
% I3 -> freq index
[~,~,I1]=intersect(elec,relec); 
I2=find((rtime>=time(1) & rtime<=time(2)));
I3=find((rfreq>=freq(1) & rfreq<=freq(2)));

% frej=rdata.frej;
out=squeeze(mean(mean(mean(rdata.freqlock.powspctrm(:,I1,I3,I2),2),3),4));

% if isfield(ST.dat{ii},'mir')
%     [~,~,I1]=intersect(ST.dat{ii}.mir,elec); 
%     data2=squeeze(mean(mean(mean(rdata.fdata.powspctrm(:,I1,I3,I2),2),3),4));
% end
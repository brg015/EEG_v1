function freq_method_check(D,freqlock)

ch=38; f=8; 

spct=D.powspctrm;
spct=squeeze(spct(~D.rej,ch,f,:));
plot(D.time,mean(spct)); xlabel('time'); ylabel('value'); grid; grid minor;

display('time');
u=nanmean(nanmean(spct)); display(['u = ' num2str(u)]);
s=nanstd(nanmean(spct)); display(['std = ' num2str(s)]);

display('trial');

a1=squeeze(freqlock.powspctrm(~D.rej,ch,f,:));
a2=log10(a1);
a3=zscore(a2);

for ii=1:length(a1)
    a4(ii,:)=(a2(ii,:)-nanmean(a2))./nanstd(a2);
end
function [A1,A2,A3]=EEG_ERS_Plot(E1,E2,M1,M2,t)

if ~isempty(E1)
figure(1);
subplot(2,2,1);
plot(t,mean(E1),'b'); hold on; plot(t,mean(E2),'r');
ylabel('norm(uV)'); xlabel('Time(s)'); grid;
subplot(2,2,2);
[h,p,ci,stats]=ttest(E1,E2,'tail','right');
plot(t,stats.tstat); 
ylabel('T-value'); xlabel('Time(s)'); grid;

subplot(2,2,3);
plot(t,mean(M1),'b'); hold on; plot(t,mean(M2),'r');
ylabel('Mean R2'); xlabel('Time(s)'); grid;
subplot(2,2,4);
[h,p,ci,stats]=ttest(M1,M2,'tail','right');
plot(t,stats.tstat); 
ylabel('T-value'); xlabel('Time(s)'); grid;

else
%     subplot(2,2,1);
%     plot(t,mean(M1),'b');
%     ylabel('Mean R2'); xlabel('Time(s)'); grid;
%     
    [h,p,ci,stats]=ttest(M1,0,'tail','right');
%     subplot(2,2,2);
%     plot(t,stats.tstat); 
%     ylabel('T-value'); xlabel('Time(s)'); grid;
%     
%     subplot(2,2,3);
%     plot(t,p); 
%     ylabel('p-value'); xlabel('Time(s)'); grid;
%     
%     subplot(2,2,4);
%     plot(t,h); 
%     ylabel('p-value'); xlabel('Time(s)'); grid;
    
    A1=mean(M1);
    A2=stats.tstat;
    A3=p;
end
function dt_plot(data,EEG,trialStruct,ch1,ch2)

ST=1;

if isempty(EEG), eo=1; else eo=0; end

I1=find(trialStruct.rProfile(:,4)==1); % Hits
I2=find(trialStruct.rProfile(:,5)==1); % Misses
IC=find(trialStruct.rejthresh==1);
IE=find(trialStruct.event~=2);
I1=setdiff(I1,IC); I1=setdiff(I1,IE);
I2=setdiff(I2,IC); I2=setdiff(I2,IE);

t1=zeros(1,1000); % EEG hits ch1
t2=zeros(1,1000); % FT  hits ch2
t3=zeros(1,1000); % EEG miss ch1
t4=zeros(1,1000); % FT  miss ch2

c1=1; for ii=I1'
    if eo~=1, t1=t1+squeeze(EEG.data(ch1,:,ii)); end
    t2=t2+data.trial{ii}(ch2,:);
    c1=c1+1;
end
c2=1; for ii=I2'
    if eo~=1, t3=t3+squeeze(EEG.data(ch1,:,ii)); end
    t4=t4+data.trial{ii}(ch2,:);
    c2=c2+1;
end

t1=t1./(c1-1);
t2=t2./(c1-1);
t3=t3./(c2-1);
t4=t4./(c2-1);

figure(10); title(['ch1:' num2str(ch1) ' - ch2:' num2str(ch2)]);
if eo~=1, plot(t1,'b+'); hold on; end
plot(t2,'b'); hold on;
if eo~=1, plot(t3,'r+'); hold on; end
plot(t4,'r'); 

figure(11); 
c=1; for ii=I1'
    V(c,:)=data.trial{ii}(ch2,1:500)-mean(data.trial{ii}(ch2,1:200)); c=c+1;
    figure; plot(V(c-1,:)); pause; close all;
end
figure(3); plot(abs(std(V)),'k')
hold on; plot(mean(V)-std(V),'k--');
hold on; plot(mean(V)+std(V),'k--');    
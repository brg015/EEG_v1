function return_val=lin3D(X,i1,i2,i3,dat,cstr)
% X   => 3D Matrix
% i1  => channel
% i2  => freq
% i3  => time
% dat => data info
%=========================================================================%
dat_file='eeg';
switch dat_file
    case 'eeg'
        c=dat.label;
        t=dat.time;
        f=dat.freq;
    otherwise     
end
%=========================================================================%
% Pull channel information
v1=[];
for ii=1:length(i1)
    switch class(i1{ii})
        case 'double', I=num2str(i1{ii});
        otherwise I=i1{ii};
    end
    v1=[v1, find(strcmp(I,c))];
end
% Reduce X
if length(v1)>1
    X1=squeeze(mean(X(v1,:,:)));
else
    X1=squeeze(X(v1,:,:));
end
clear I v1 X;
%=========================================================================%
% Frequency Info
I=find(f>=i2(1) & f<=i2(2));
if length(I)>1
    X2=squeeze(mean(X1(I,:)));
else
    X2=squeeze(X1(I,:));
end
clear I X1;
%=========================================================================%
% Time
I=find(t>=i3(1) & t<=i3(2));
X3=X2(I); 
clear X2;
plot(t(I),X3,cstr,'linewidth',3);
set(gcf,'color','w');
set(gca,'fontsize',14);
grid;

return_val=X3;
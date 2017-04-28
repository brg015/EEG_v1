function fancy_eeg_adj(X,S,fband,tband,chband,save_name)
% Input variables include
% X => Coherence infor
% S => Spectral data
% fband => Frequ band to include
% tband => timeperioed to incldue (in ms)

% Declare global RUN
global RUN;

% Define save profile and skip if run
save_file=fullfile(RUN.adv.save,[save_name '.mat']);
% if exist(save_file,'file'), return; end

% Determine size of data (trial chan freq time)
s=size(S);

% determine time and freq
if isempty(fband)
    f=ones(1,s(3));
else
    f=and(X.freq>=fband(1),X.freq<=fband(2));
end

if isempty(tband)
    t=ones(1,s(4));
else
    t=and(X.time>=tband(1)/1000,X.time<=tband(2)/1000);
end

if isempty(chband)
    ch=1:s(2)-1;
else
    ch=chband;
end

% Prep S for weights
% 1) logsig transform
for ii=1:s(1),
    for jj=1:s(2),
        for kk=1:s(3),
            Ssigm(ii,jj,kk,:)=logsig(S(ii,jj,kk,:));
        end
    end
end
clear S;
% Never comibne logical and index
mS=squeeze(mean(Ssigm));
smS=mS(ch,find(f),find(t));
L=X.labelcmb;        % Get Combos
C=unique(L);         % Get Labels
C=C(ch);
for ii=1:size(L,1) % All combos (ignore RM)
    x=strcmp(L{ii,1},C);
    y=strcmp(L{ii,2},C);
    if ~isempty(x) && ~isempty(y)
        V=squeeze(X.cohspctrm(ii,:,:));
        c1=1;
        for ff=find(f)
            c2=1;
            for tt=find(t)
                XR{c1,c2}(x,y)=V(ff,tt);
                XR{c1,c2}(y,x)=V(ff,tt);
                c2=c2+1;
            end
            c1=c1+1; 
        end
    end
end

%=========================================================================%
% Now the hard part
%=========================================================================%
% 1) We know that the labels match, both came from freqlock 
% 2) We can also treat time and space as equal
% 3) But distance can't exceed 1...

% XR  {f,t}(ch X ch)
% smS [ch X f X t]    z -> sigma transformed

% Intialize R & Box
R=zeros(length(ch)*numel(XR));
Box=ones(1,numel(XR))*length(ch);
% Reshape XR
c=1;
for ii=1:size(XR,1)
    for jj=1:size(XR,2)
        XRL{c}=XR{ii,jj}; 
        smSL(:,c)=squeeze(smS(:,ii,jj));
        I(c,:)=[ii,jj];
        c=c+1;
    end
end
clear XR smS;
% frequency then time info... this makes the diagnol straightforward

R=zeros(size(R));

% => Gonna loop time then space
for ii=1:length(Box)       % Row
    ct=1; % Time counter
    for jj=1:length(Box)   % Col
        
        % Gotta determine block-type
        [r,c]=cell_block(ii,jj,Box);
        md=I(ii,:);   % Diag match
        m=I(ct,:);    % Current
        %======================%
        % On Fat Diag
        %======================%
        if ii==jj
            R(r,c)=XRL{ct};
        end
        %======================%
        % Time change
        %======================%
        % if freq is same and time is one forward
        if ((md(1)==m(1)) && (m(2)-md(2)==1))
           a=eye(length(r),length(c));
           % Project power of prior match forward in time
           a(a==1)=smSL(:,ct-1);
           R(r,c)=a;
           clear a;
        end
        %======================%
        % Freq change
        %======================%
        % if time is equal and freq is one forward
        if ((md(2)==m(2) && (m(1)-md(1)==1)))
           a=eye(length(r),length(c));
           % Project power of prior match forward in freq
           cnct(1)=md(1)+1; cnct(2)=md(2);
           linker=ismember(I,cnct,'rows');
           a(a==1)=smSL(:,linker);
           R(r,c)=a;
           clear a linker cnct;
        end
        
        % if time is equal and freq is one back
        if ((md(2)==m(2) && (md(1)-m(1)==1)))
           a=eye(length(r),length(c));
           % Project power of prior match back in freq
           cnct(1)=md(1)-1; cnct(2)=md(2);
           linker=ismember(I,cnct,'rows');
           a(a==1)=smSL(:,linker);
           R(r,c)=a;
           clear a;
        end
        
        clear r c md m;
        ct=ct+1;
    end
end

% Save that output
dat.R=R;
dat.freq=X.freq(f);
dat.time=X.time(t);
dat.I=I;
dat.Box=Box;

save(save_file,'dat');

















function COH=eeg_coherence_network(fn,iConditions,gcohlock)

Data=gcohlock.(fn{iConditions});
L=Data{1}.labelcmb; % Get Combos
C=unique(L);        % Get Labels

for s=1:length(Data)
    for ii=1:size(L,1)
        x=strcmp(L{ii,1},C);
        y=strcmp(L{ii,2},C);
        V=squeeze(Data{s}.cohspctrm(ii,:,:));
        for f=1:length(Data{s}.freq)
            for t=1:length(Data{s}.time)
                COH.V{s,f,t}(x,y)=V(f,t);
            end
        end
    end
    if s==1,
        COH.time=Data{s}.time;
        COH.freq=Data{s}.freq;
        COH.ch=C;
    end
end
% Remove lateral eyemovements

function [rej,M] = EEGartstep(EEG)

% Generate step filter
% sr = 250 Hz
% 50 pt step (200ms)
% 30 uv step (30uv)
sstep=50;
smax=20;
sfunc=[ones(1,sstep).*-1 ones(1,sstep)];

% Loop through all trials
for ii=1:size(EEG.data,3)
    x1=EEG.data(21,:,ii);
    x2=EEG.data(22,:,ii);
    x12=x1-x2; % Lateral frontal electrode difference
    
    C=abs(conv(sfunc,x12)./50);
    M(ii)=max(C);
    if max(C)>smax
        rej(ii)=1;
    else
        rej(ii)=0;
    end
end

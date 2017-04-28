function [data]  = contraMinusIpsi(layoutmw64, dataLeft, dataRight, powVSerp)
%-------- channel index for left and right ----%

if isempty(dataLeft)
    V=dataRight;
else
    V=dataLeft;
end

chanIdx = zeros(length(layoutmw64.channels.conIpsi),2);
for iChans = 1:length(layoutmw64.channels.conIpsi)
    for i = 1:2
        chanIdx(iChans,i) = find(strcmp(V.label,layoutmw64.channels.conIpsi(iChans,i)));
    end
end
%-------- channel index for middle electrodes ----%
chanIdxMiddle = zeros(length(layoutmw64.channels.middle),1);
for i = 1:length(layoutmw64.channels.middle)
    chanIdxMiddle(i) = find(strcmp(V.label,layoutmw64.channels.middle(i)));
end;


%----- put it in a more convenient structure-----%
electrodes.left =  chanIdx(:,1);
electrodes.right =  chanIdx(:,2);
electrodes.middle = chanIdxMiddle;

switch powVSerp
    case 'power', fn='powspctrm';
    case 'erp', fn='individual';
end
            
if (~isempty(dataRight) && ~isempty(dataLeft))
    conIpsi = dataLeft;
    % R = [(L,L + R,R) - (R,L + L,R)] / 2
    % L = [(R,L + L,R) - (L,L + R,R)] / 2
    % Mid = 0;
    conIpsi.(fn)(:,electrodes.right,:) = 0.5 *((dataLeft.(fn)(:,electrodes.left,:) + ...
        dataRight.(fn)(:,electrodes.right,:))-...
        (dataRight.(fn)(:,electrodes.left,:)+ ...
        dataLeft.(fn)(:,electrodes.right,:)));

    conIpsi.(fn)(:,electrodes.left,:) = 0.5 *((dataRight.(fn)(:,electrodes.left,:)+ ...
        dataLeft.(fn)(:,electrodes.right,:))-...
        (dataLeft.(fn)(:,electrodes.left,:)+ ...
        dataRight.(fn)(:,electrodes.right,:)));

    conIpsi.(fn)(:,electrodes.middle,:) = dataLeft.(fn)(:,electrodes.middle,:)- ...
        dataLeft.(fn)(:,electrodes.middle,:);
elseif isempty(dataRight)
    % R (ipsi) = [L - R]
    % L (contra) = [R - L]
    conIpsi=dataLeft;
    conIpsi.(fn)(:,electrodes.right,:) = ((dataLeft.(fn)(:,electrodes.left,:) - ...
        dataLeft.(fn)(:,electrodes.right,:)));

    conIpsi.(fn)(:,electrodes.left,:) = ((dataLeft.(fn)(:,electrodes.right,:) - ...
        (dataLeft.(fn)(:,electrodes.left,:))));

    % Central calc is the same i.e. just zeroing
    conIpsi.(fn)(:,electrodes.middle,:) = dataLeft.(fn)(:,electrodes.middle,:)- ...
        dataLeft.(fn)(:,electrodes.middle,:);
elseif isempty(dataLeft)
    conIpsi = dataRight;

    conIpsi.(fn)(:,electrodes.right,:) = ((dataRight.(fn)(:,electrodes.right,:)) - ...
        (dataRight.(fn)(:,electrodes.left,:)));

    conIpsi.(fn)(:,electrodes.left,:) = 0.5 *((dataRight.(fn)(:,electrodes.left,:) - ...
        dataRight.(fn)(:,electrodes.right,:)));

    conIpsi.(fn)(:,electrodes.middle,:) = dataRight.(fn)(:,electrodes.middle,:)- ...
        dataRight.(fn)(:,electrodes.middle,:);
end

selchan = ft_channelselection({'all' '-LM' '-LVEOG'},V.label); % '-RM' 
conIpsi = ft_selectdata(conIpsi,'channel',selchan);

data = conIpsi;
 
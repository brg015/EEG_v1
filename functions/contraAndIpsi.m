function [dataContra, dataIpsi]  = contraAndIpsi(layoutmw64, dataLeft, dataRight, powVSerp)
%-------- channel index for left and right ----%


chanIdx = zeros(length(layoutmw64.channels.conIpsi),2);
for iChans = 1:length(layoutmw64.channels.conIpsi)
    for i = 1:2
        chanIdx(iChans,i) = find(strcmp(dataLeft.label,layoutmw64.channels.conIpsi(iChans,i)));
    end
end
%-------- channel index for middle electrodes ----%
chanIdxMiddle =zeros(length(layoutmw64.channels.middle),1);
for i = 1:length(layoutmw64.channels.middle)
    chanIdxMiddle(i) = find(strcmp(dataLeft.label,layoutmw64.channels.middle(i)));
end;


%----- put it in a more convenient structure-----%
electrodes.left =  chanIdx(:,1);
electrodes.right =  chanIdx(:,2);
electrodes.middle = chanIdxMiddle;

switch powVSerp
    case 'power'
        
        contra = dataLeft;
        contra.powspctrm(:,electrodes.right,:,:) = 0.5 * ((dataRight.powspctrm(:,electrodes.left,:,:) + dataLeft.powspctrm(:,electrodes.right,:,:)));
        contra.powspctrm(:,electrodes.middle,:,:) = dataLeft.powspctrm(:,electrodes.middle,:,:) - dataLeft.powspctrm(:,electrodes.middle,:,:);
        contra.powspctrm(:,electrodes.left,:,:) = 0.5 * ((dataRight.powspctrm(:,electrodes.left,:,:) + dataLeft.powspctrm(:,electrodes.right,:,:)));
        
        ipsi = dataLeft;
        ipsi.powspctrm(:,electrodes.right,:,:) = 0.5 * ((dataRight.powspctrm(:,electrodes.right,:,:) + dataLeft.powspctrm(:,electrodes.left,:,:)));
        ipsi.powspctrm(:,electrodes.middle,:,:) = dataLeft.powspctrm(:,electrodes.middle,:,:) - dataLeft.powspctrm(:,electrodes.middle,:,:);
        ipsi.powspctrm(:,electrodes.left,:,:) = 0.5 * ((dataRight.powspctrm(:,electrodes.right,:,:) + dataLeft.powspctrm(:,electrodes.left,:,:)));
        
        %contra.powspctrm(:,33,:,:) = zeros(size(contra.powspctrm(:,33,:,:)));
        %ipsi.powspctrm(:,33,:,:) = zeros(size(ipsi.powspctrm(:,33,:,:)));
        
        selchan = ft_channelselection({'all' '-LM' '-LVEOG'},dataLeft.label); % '-RM' 
        contra = ft_selectdata(contra,'channel',selchan);
        ipsi = ft_selectdata(ipsi,'channel',selchan);
        
        dataContra = contra;
        dataIpsi = ipsi;
        
    case 'erp'
        
        contra = dataLeft;
        contra.individual(:,electrodes.right,:) = 0.5 * ((dataRight.individual(:,electrodes.left,:) + dataLeft.individual(:,electrodes.right,:)));
        contra.individual(:,electrodes.middle,:) = dataLeft.individual(:,electrodes.middle,:) - dataLeft.individual(:,electrodes.middle,:);
        contra.individual(:,electrodes.left,:) = 0.5 * ((dataRight.individual(:,electrodes.left,:) + dataLeft.individual(:,electrodes.right,:)));
        
        ipsi = dataLeft;
        ipsi.individual(:,electrodes.right,:) = 0.5 * ((dataRight.individual(:,electrodes.right,:) + dataLeft.individual(:,electrodes.left,:)));
        ipsi.individual(:,electrodes.middle,:) = dataLeft.individual(:,electrodes.middle,:) - dataLeft.individual(:,electrodes.middle,:);
        ipsi.individual(:,electrodes.left,:) = 0.5 * ((dataRight.individual(:,electrodes.right,:) + dataLeft.individual(:,electrodes.left,:)));
        
        selchan = ft_channelselection({'all' '-LM' '-LVEOG'},dataLeft.label); % '-RM' 
        contra = ft_selectdata(contra,'channel',selchan);
        ipsi = ft_selectdata(ipsi,'channel',selchan);
        
        dataContra = contra;
        dataIpsi = ipsi;
        
    otherwise
        disp('input power or ERP measure')
end

end
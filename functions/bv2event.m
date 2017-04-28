function EEG = bv2event(EEG)

% script for removing the prefixes in the event codes e.g. ('S 1' becomes '1')
% to use: in your matlab session add the directory of this script to 
% your path: 
% addpath(genpath('dir to bv2event.m'))
% load in the EEG data: pop_loadset('dir to your eeg data')
% then call the script by typing EEG = bv2event(EEG)
% if you have any questions please email me: berryv.dberg@gmail.com


	for i = 1:length(EEG.event)
		if(strcmp(EEG.urevent(1,4).code,'Stimulus'))
    			temp = strrep(EEG.event(i).type,'S','');
    			temp = strrep(temp,' ','');
    
    			EEG.urevent(i).type  = temp;
   		 	EEG.event(i).type = temp;
    		end
	end
end


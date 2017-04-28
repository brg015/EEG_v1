function EEG=eeg_reference(EEG)

global RUN;
%%%rereferencing
% Reref in EEGlab to average mastoid.
% first step is to add a ref channel to the channel location structure (its
% basically an empty "dummy" channel). channel locations are not needed at
% this point
EEG_1=pop_chanedit(EEG, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
    RUN.template.elp),'append',63, ...
    'changefield',{64 'labels' 'RM'},'setref',{'1:63' '63'});
% Changes
% Add an extra channel after 63              'append',63
% Changefield labels of added channel to RM  'changefield',{64 'labels' 'RM'}
% Set ref of all channels to be 63           'setref',{'1:63' '63'}
% - Has NO impact on the actual data

% rereference to average over all channels (in pop_reref the second input; [] does
% this), 'refloc' with label Ch33 retains the old reference, it basically
% fills in the dummy channel as the old reference channel....
EEG_2 = pop_reref( EEG_1, [],'refloc',struct('labels',{'RM'},'theta',{[]},'radius',{[]},'X',{[]},'Y',{[]},...
    'Z',{[]},'sph_theta',{[]},'sph_phi',{[]},'sph_radius',{[]},'type',{[]},'ref',{[]},'urchan',{[]}));
% Changes
% Label 33 is the old reference (channel 32)
% An additional channel has been added
% Ref is now averef and not common

% rereference to average mastoid NOT KEEPING the reference channels
% index channel 32 is the Right mastoid (channel 33 on the cap) 64 is the
% dummy channel we added in earlier and becomes the left mastoid...
%EEG = pop_reref( EEG, [32 64] );
% rereference to average mastoid KEEPING the reference channels

switch RUN.set.reref
    case 'average',  EEG_3 = pop_reref( EEG_2, 1:64 ,'keepref','on');
    case 'implicit', EEG_3 = pop_reref( EEG_2, 32, 'keepref','on');
    case 'mastoid',  EEG_3 = pop_reref( EEG_2, [32 64],'keepref','on');
    case 'r_mastoid',  EEG_3 = pop_reref( EEG_2, [64],'keepref','on');
    case 'l_mastoid',  EEG_3 = pop_reref( EEG_2, [32],'keepref','on');
    case 'Cz',         EEG_3 = pop_reref( EEG_2, [4] ,'keepref','on');
    otherwise
        error('Incorrect reference is set!');
end

EEG_4 = pop_chanedit(EEG_3, 'lookup',fullfile(fileparts(which('eeglab.m')), ...
    RUN.template.elp),'load',{RUN.template.ced 'filetype' 'autodetect'});

% Clear original and recast
clear EEG; EEG=EEG_4;
function EEG=eeg_SEFER(EEG,iSubj,scenario)
global RUN;

switch RUN.dir.sess
    case 'enc', 
        check_stim=1; Rresp=[5,6,10,9];
    case 'ret', 
        check_stim=[2,3]; Rresp=[4 5 6 7 8 9 10 11];
    case 'err'
        check_stim=[2,3]; Rresp=[4 5 6 7 8 9 10 11];
end

switch scenario
    case 'epoch'
%=========================================================================%
% Epoch data
%=========================================================================%
if (strcmp(RUN.dir.sess,'ret') && strcmp(RUN.dir.subjects{iSubj},'2735'))
    EEG=fix2735ret(EEG); % Special fix for our pilot subject.
end

ir_event='999';

% 1) Read in event data form behave
data=excel_reader(fullfile(RUN.dir.beh,['subject' RUN.dir.subjects{iSubj} '.csv']));
for ii=1:length(data), head{ii}=data{ii}.header{1}; end
% 2) Grab codes from here
codes=[data{strcmp(head,'code1')}.col; data{strcmp(head,'code2')}.col]'; 
for ii=1:size(codes,1)
    for jj=1:size(codes,2)
        codes_mat(ii,jj)=str2double(codes{ii,jj});
    end
end
% 3) Now define event epochs
events_idx=data{strcmp(head,'ID')}.col;

% 4) Modify code drops to better reflect events, we're looking
% at event and urevent for this

for ii=1:length(EEG.event), event_code(ii)=str2double(EEG.event(ii).type); end
% 5) Now find codes in EEG
I=zeros(1,length(event_code));
for ii=1:length(check_stim), I=I+(event_code==check_stim(ii)); end
I=find(I);

try
EEG_codes=[event_code(I+1); event_code(I+2)]';
catch err
    keyboard
end
event_count=1; inc_count=1; c=1; R=1;

for ii=1:size(EEG_codes,1) % Loop through EEG finds
    idx=find(ismember(codes_mat,EEG_codes(ii,:),'rows')); % Find code location
    if ~isempty(idx) % if found\
        if strcmp(EEG.event(I(ii)).type,'1')
            E_plus=10000;   E_trial=40000;
        elseif strcmp(EEG.event(I(ii)).type,'2')
            E_plus=20000;   E_trial=50000;
        elseif strcmp(EEG.event(I(ii)).type,'3')
            E_plus=30000;   E_trial=60000;
        else
            error('You messed up');
        end
        % Set start event to idx value + 10000
        EEG.event(I(ii)).otype=EEG.event(I(ii)).type;
        EEG.event(I(ii)).type=num2str(str2double(events_idx(idx))+ E_plus); 
        % End event is + 2 away, this is trial start
        EEG.event(I(ii)+2).otype=EEG.event(I(ii)+2).type;
        EEG.event(I(ii)+2).type=num2str(str2double(events_idx(idx))+ E_trial);

        inc_events(inc_count)=I(ii);    inc_count=inc_count+1;
        inc_events(inc_count)=I(ii)+2;  inc_count=inc_count+1;
        
        % If events are included then...
        if (length(event_code)>=(I(ii)+3))
            pos_resp=event_code(I(ii)+3);    
            if (sum(Rresp==pos_resp) && (EEG.event(I(ii)+3).latency-EEG.event(I(ii)).latency)<5*1000)
                RT(R)=EEG.event(I(ii)+3).latency-EEG.event(I(ii)).latency; R=R+1;
                inc_events(inc_count)=I(ii)+3; inc_count=inc_count+1;
            end
        end
        event_count=event_count+1;

        events_epoch{c}=num2str(str2double(events_idx(idx))+ E_plus);
        events_epoch{c+1}=num2str(str2double(events_idx(idx))+ E_trial);
        c=c+2;
    end
end

for jj=1:length(Rresp)
    events_epoch{end+1}=num2str(Rresp(jj));
end

% 6) Set other events to 'bad'
bad_events=setdiff(1:length(EEG.event),inc_events);
for ii=1:length(bad_events)
    % urevent is merely history i.e. not need updated
    % EEG.urevent(bad_events(ii)).otype=EEG.urevent(bad_events(ii)).type;
    % EEG.urevent(bad_events(ii)).type=ir_event;
    EEG.event(bad_events(ii)).otype=EEG.event(bad_events(ii)).type;
    EEG.event(bad_events(ii)).type=ir_event;
end   

EEG.events_epoch=events_epoch;

end
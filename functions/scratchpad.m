













%%

events = zeros(1,length(EEG.urevent));
times = zeros(1,length(EEG.urevent));


for i= 1:length(EEG.urevent)
    if(~strcmp(EEG.urevent(i).type,'boundary'))
        events(i) = str2double(EEG.urevent(i).type);
        times(i) = EEG.urevent(i).latency ;
    end
end

sum(ismember(events,100:132))
sum(ismember(events,1:4))


trialStruct.targetType = events(ismember(events,100:132))';
trialStruct.targetLatency = times(ismember(events,100:132))';
trialStruct.cueType = zeros(length(trialStruct.targetType),1);
trialStruct.cueSOA = zeros(length(trialStruct.targetType),1);
trialStruct.response = zeros(length(trialStruct.targetType),1);
trialStruct.responseLatency = zeros(length(trialStruct.targetType),1);
trialStruct.RT=zeros(length(trialStruct.targetType),1);

% How to calculate latency
% response Latency - target latency
%response = 1-4, target = 100-132


%%
for i = 1:length(trialStruct.targetType)
    trialStruct.cueType(i) = events(times < trialStruct.targetLatency(i)-500 & times > ...
        trialStruct.targetLatency(i)-1700 & ismember(events,40:44 ));
    
    trialStruct.cueSOA(i) = times(times == trialStruct.targetLatency(i))-times(times < ...
        trialStruct.targetLatency(i)-500 & times > trialStruct.targetLatency(i)-1700 ...
        & ismember(events,40:44 ));
    
   
    tmpEvents = times > trialStruct.targetLatency(i) & times < trialStruct.targetLatency(i)+600 ...
        & ismember(events,1:4);
    %between target and 1200ms posttarget, inserts a 1 if more than 1?
    tmpEventsIdx = find(tmpEvents);
    
    
    %if there is more than one response in time window: record first RT
    %and insert a 5 into responses -> should return 500 total responses
    if sum(tmpEvents) > 1
        trialStruct.responseLatency(i) = times(tmpEventsIdx(1));
        trialStruct.response(i) = 5;
    end
    
    %if there is a 'miss', record a zero
    %this is not succesfully recording misses.
    
    if sum(tmpEvents) == 0
        trialStruct.responseLatency(i) = trialStruct.targetLatency(i);
        %set these equal so that subtraction makes 0
        trialStruct.response(i) = 0;
    end
    
    
    %if there is only one response, record the response
    if sum(tmpEvents) == 1
        trialStruct.responseLatency(i) = times(tmpEventsIdx);
        trialStruct.response(i) = events(tmpEventsIdx);
    end
    
    %make sure latencies are multiplied correctly
    trialStruct.RT(i) = 2*(trialStruct.responseLatency(i)-trialStruct.targetLatency(i))
end

%%
%valid VS invalid
%when "cueType" port codes = trialType

trialStruct.valid = zeros((length(trialStruct.cueType)),1);
trialStruct.invalid = zeros((length(trialStruct.cueType)),1);

for i = 1:length(trialStruct.cueType)
    if trialStruct.cueType(i)==41 || trialStruct.cueType(i)==42
      trialStruct.valid(i) = trialStruct.RT(i);
    end
   
    if trialStruct.cueType(i) == 43 || trialStruct.cueType(i)==44
        trialStruct.invalid(i) = trialStruct.RT(i);
    end
    
end

%REMOVE ALL THE ZEROS: a(a==0) = [];
trialStruct.valid(trialStruct.valid==0) = [];
trialStruct.invalid(trialStruct.invalid==0) = [];
mean(trialStruct.valid) %+std(trialStruct.valid)
mean(trialStruct.invalid)

%%
%REWRITE in case "cueType" port codes = cueType(41 or 42)
trialStruct.valid = zeros((length(trialStruct.cueType)),1);
trialStruct.invalid = zeros((length(trialStruct.cueType)),1);


for i = 1:length(trialStruct.cueType)
    if trialStruct.targetType(i) < 117 ...
            && trialStruct.cueType(i) == 41
        %target on left and cue on left
           trialStruct.valid(i) = trialStruct.RT(i);
    end
    if trialStruct.targetType(i) > 116 && trialStruct.cueType(i) == 42
        trialStruct.valid(i) = trialStruct.RT(i);
    else
        trialStruct.invalid(i) = trialStruct.RT(i);
    
    end
    
end        

trialStruct.valid(trialStruct.valid==0) = [];
trialStruct.invalid(trialStruct.invalid==0) = [];
mean(trialStruct.valid)
mean(trialStruct.invalid)




%%
%congruent vs incongruent

trialStruct.cong=zeros((length(trialStruct.cueType)/2),1)
trialStruct.incong=zeros((length(trialStruct.cueType)/2),1)

for i = 1:length(trialStruct.response)
    if trialStruct.cueType(i)==41 | trialStruct.cueType(i)==43
      trialStruct.cong(i) = trialStruct.RT(i)
    end
   
    if trialStruct.cueType(i) == 42 | trialStruct.cueType(i)==44
        trialStruct.incong(i) = trialStruct.RT(i)
    end
    
end

trialStruct.cong(trialStruct.cong==0) = [];
trialStruct.incong(trialStruct.incong==0) = [];

mean(trialStruct.cong)
mean(trialStruct.incong)

function report_bool_trial(data,booleans,event)

fn=fieldnames(booleans);
idx=strsearch(event,fn);

bool_array=[];
for jj=1:length(idx)
    if jj==1 
        bool_array=booleans.(fn{idx(jj)}); 
    else
        bool_array=bool_array+booleans.(fn{idx(jj)});
    end
end
bool_array=bool_array>0;

trials=unique(data.trialStruct.trial(bool_array))
numel(trials)
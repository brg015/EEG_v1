%% ------------- redo positions ---------------%%
% helper file to redo the positions  in the 



%% first prepare the layout
% 
% cfg = [];cfg.rotate = 90; cfg.elec = data.elec;layoutmw64 =  ft_prepare_layout(cfg)
% cfg.layout = layoutmw64; ft_layoutplot(cfg);

%%
params.layout = layoutmw64;
params.height = 0.08;
params.width = 0.15;

params.setxparams = [-1 1];
tmp = (params.setxparams(1,2)-params.setxparams(1,1))/10;
params.setx = params.setxparams(1):tmp:params.setxparams(2);
params.setx = [params.setx(6) params.setx((-5:-1)*-1) params.setx(7:11)] ;

params.setyparams = [-0.58 0.58
    -0.6 0.6
    -0.45 0.45
    -0.4 0.4
    -0.3 0.3
    -0.4 0]



params.channels.middle = {'1' '2' '3' 'Cz' '38' '37' '36' '35' '34'};
params.channels.leftOne = {'5' '9' '11' '13' '49' '47' '45' '43' '41' '39'};
params.channels.leftTwo = {'7' '19' '59' '57' '55' '53' '51'};
params.channels.leftThree = {'15' '17' '25' '63' '61'};
params.channels.leftFour = {'21' '23' '27' 'LM'};
params.channels.leftFive = {'LVEOG' '29'};

params.channels.rightOne = {'6' '10' '12' '14' '50' '48' '46' '44' '42' '40'};
params.channels.rightTwo = {'8' '20' '60' '58' '56' '54' '52'};
params.channels.rightThree = {'16' '18' '26' '64' '62'};
params.channels.rightFour = {'22' '24' '28' 'RM'};
params.channels.rightFive = {'RVEOG' '30'};


params.sety.middle = params.setyparams(1,1):(params.setyparams(1,2) -params.setyparams(1,1))/8:params.setyparams(1,2);
params.sety.leftOne = params.setyparams(2,1):(params.setyparams(2,2) -params.setyparams(2,1))/9:params.setyparams(2,2);
params.sety.leftTwo = params.setyparams(3,1):(params.setyparams(3,2) -params.setyparams(3,1))/6:params.setyparams(3,2);
params.sety.leftThree = params.setyparams(4,1):(params.setyparams(4,2) -params.setyparams(4,1))/4:params.setyparams(4,2);
params.sety.leftFour = params.setyparams(5,1):(params.setyparams(5,2) -params.setyparams(5,1))/3:params.setyparams(5,2);
params.sety.leftFive = params.setyparams(6,1):(params.setyparams(6,2) -params.setyparams(6,1))/1:params.setyparams(6,2);
params.sety.rightOne = params.sety.leftOne;
params.sety.rightTwo = params.sety.leftTwo;
params.sety.rightThree = params.sety.leftThree;
params.sety.rightFour = params.sety.leftFour;
params.sety.rightFive = params.sety.leftFive;



fn = fieldnames(params.channels);
for i = 1:length(fn)
    params.idx.(fn{i}) =zeros(length(params.channels.(fn{i})),1);
    for j = 1:length(params.channels.(fn{i}))
        params.idx.(fn{i})(j) = find(strcmp(params.layout.label,params.channels.(fn{i})(j)));
    end
end
params.layout.width = params.layout.width + params.width;
params.layout.height = params.layout.height + params.height;

fn = fieldnames(params.idx);

for i = 1:length(fn);
    disp(fn{i})
    
    params.layout.pos(params.idx.(fn{i}),1) =  params.setx(i);% xposition
    params.layout.pos(params.idx.(fn{i}),2) =  -1*params.sety.(fn{i})';% yposition

    
end

cfg.layout = params.layout; ft_layoutplot(cfg);

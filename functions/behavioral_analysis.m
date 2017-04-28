%%plotting behavioral data
%RT - invalid vs valid


xlabel = ('Trial Type');
y = [mean_RT_Valid mean_RT_Invalid];
figure; bar(y);

%%
%validity effect
%RT - invalid vs valid + errorbar
y = [mean_RT_Valid mean_RT_Invalid];
e = [sd_RT_Valid sd_RT_Invalid];
figure; errorbar(y,e)

%multiple lines - how to change x axis to valid, invalid
x =[1 2];
y = [mean_RT_Valid mean_RT_Invalid];
z = [mean_RT_Cong mean_RT_Incong];
figure; plot(x, y, 'blue', x, z, 'red');
legend('valid', 'invalid');


%%
%congruency effect
%RT - incongruent vs congruent + errorbar

%insert title using GUI
%LATER: test for significance

figurepalette('show')


y = [mean_RT_Cong mean_RT_Incong];
e = [sd_RT_Cong sd_RT_Incong];
figure; errorbar(y,e)
xlabel('Trial Type');
ylabel('Reaction Time (ms)');

%%
%misc, examples


%BATRGRAPH with error bars - from google
%@Jaroslav Hlinka. Fixed this bug. Changed the plot errors loop to this. Set 'error_sides' to 3.
% Plot erros 

for i = 1:numbars 
x = get(get(handles.bars(i),'children'), 'xdata'); 
x = mean(x([1 3],:)); 

if error_sides == 3 
pValues = barvalues; 
pValues(pValues<0)=NaN; 
nValues = barvalues; 
nValues(nValues>0)=NaN; 
handles.posErrors(i) = errorbar(x, pValues(:,i), zeros(size(errors(:,i))), errors(:,i), 'k', 'linestyle', 'none', 'linewidth', 2); 
handles.negErrors(i) = errorbar(x, nValues(:,i), errors(:,i), zeros(size(errors(:,i))), 'k', 'linestyle', 'none', 'linewidth', 2); 
else 
handles.errors(i) = errorbar(x, barvalues(:,i), errors(:,i), 'k', 'linestyle', 'none', 'linewidth', 2); 
end 
ymax = max([ymax; barvalues(:,i)+errors(:,i)]); 
end





%Plot multiple vectors on one graph 
X = [3 9 27]; % dependent vectors of interest
Y = [10 8 6];
Z = [4 4 4];
t = [1 2 3]; % independent vector
figure
hold on % allow all vectors to be plotted in same
 % figure
plot(t, X, ‘blue’, t, Y, ‘red’, t, Z, ‘green’)
title(‘Plot of Distance over Time’) % title
ylabel(‘Distance (m)’) % label for y axis
xlabel(‘Time (s)’) % label for x axis
legend(‘Trial 1’, ‘Trial 2’, ‘Trial 3’)
legend(‘Location’,‘NorthWest’) % move legend to upper left





line('Valid',mean_RT_Valid, 'Invalid',mean_RT_Invalid )

x= 
y= trialStruct.reactionTimevalid
y(2,:)= trialStruct.reactionTimeinvalid
plot(x,y)


%los ejemplos
x = 0:.2:20;
y = sin(x)./sqrt(x+1);
y(2,:) = sin(x/2)./sqrt(x+1);
y(3,:) = sin(x/3)./sqrt(x+1);
plot(x,y)


load count.dat;
y = mean(count,2);
e = std(count,1,2);
figure
errorbar(y,e,'xr')
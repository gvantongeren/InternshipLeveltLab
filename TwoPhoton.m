%% Two photon data
% Made by Geke van Tongeren april 2024
% 123 orientation figure
% 113 orientation ground
% 122 figure
% 112 grey hold ground

na0 = Res.ax>=0 & Res.ax <=1; % Define the timeframe 

% Indexing the rows of the different combinations and then selecting them
rows_123 = info.Stim.log(:,2) == 2 & info.Stim.log(:,3) == 3;
data123 = Res.CaSigCorrected(na0,rows_123,:); % Put the data in a new data file of only combination 123
data123m = squeeze(mean(data123)); % Average over the time
data123a = mean(data123m); % Average over trial, for every neuron
data123sd = std(data123m); % Calculate the sd 

% Repeat for the other combinations
rows_113 = info.Stim.log(:,2) == 1 & info.Stim.log(:,3) == 3;
data113 = Res.CaSigCorrected(na0,rows_113,:);
data113m = squeeze(mean(data113));
data113a = mean(data113m);
data113sd = std(data113m);

rows_122 = info.Stim.log(:,2) == 2 & info.Stim.log(:,3) == 2;
data122 = Res.CaSigCorrected(na0,rows_122,:);
data122m = squeeze(mean(data122));
data122a = mean(data122m);
data122sd = std(data122m);

rows_112 = info.Stim.log(:,2) == 1 & info.Stim.log(:,3) == 2;
data112 = Res.CaSigCorrected(na0,rows_112,:);
data112m = squeeze(mean(data112));
data112a = mean(data112m);
data112sd = std(data112m);

% Calculate the d prime for the orientation
numerator = data123a - data113a;
figkwadraat = zeros(1,length(data123sd));
groundkwadraat = zeros(1,length(data123sd));
for i=1:length(data123sd)
figkwadraat(i) = data123sd(i)*data123sd(i);
groundkwadraat(i) = data113sd(i)*data113sd(i);
denumerator = sqrt(0.5* (figkwadraat +groundkwadraat));
dprime_orientation(i) = numerator(i) / denumerator(i);
end

% Calculate the d prime for the grey hole
numerator = data122a - data112a;
figkwadraat = zeros(1,length(data122sd));
groundkwadraat = zeros(1,length(data123sd));
for i=1:length(data122sd)
figkwadraat(i) = data122sd(i)*data122sd(i);
groundkwadraat(i) = data112sd(i)*data112sd(i);
denumerator = sqrt(0.5* (figkwadraat +groundkwadraat));
dprime_greyhole(i) = numerator(i) / denumerator(i);
end
%% Script, significant response. Over all pictures all neurons
% Made by Geke van Tongeren april 2024

% One mouse
%% Script

%Define variables
nrois = size(Res.CaSigCorrected, 3); %Number of neurons
chosenROIs = 1:nrois; %making a vector of all neurons
nImgs = length(unique(Log)); %Unique images, how many stimuli were shown
ImgsRepeats = sum (Log == Log(1)); %How many times was the first (and thus any picture) repeated in total
p_values = zeros(nImgs,nrois); % making a empty structure of zeros to put in the pvalues later
v = Res.ax>=-1 & Res.ax <=1; %Defining time frame of interest.
voor0 = Res.ax >= -1 & Res.ax<=0; %Before stimulus was shown
na0 = Res.ax>0 & Res.ax <=1; %After stimulus was shown

% Control if data is right
for i=1:nImgs %For every picture
    % figure;
    % for j=chosenROIs %For every neuron
    %     plot(Res.ax(v), Res.CaSigCorrected(v,Log==i ,j)); %Plot the raw
    %     % data
    %     hold on;
    %     plot(Res.ax(v), mean(Res.CaSigCorrected(v, Log==i, j), 2), 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Average Line'); %Plot the average line of all repeats
    % end
    % Splits de gegevens in twee delen: voor en na x=0
    before0 = Res.CaSigCorrected(voor0, Log==i, :);
    after0 = Res.CaSigCorrected(na0, Log==i, :);
    % Bereken het gemiddelde voor elk deel
    before0mean = mean(before0);
    after0mean = mean(after0);
    %ttest
    [h,p_values(i,:)]=ttest2(before0mean, after0mean);
end

%Significance information
significant = p_values<=0.05; %What p values are significant
significantperPic = sum(significant,2); %Overview of which picture works well
significantperNeuron = sum(significant); %Overview of which neuron reacts well
imagesc(significant) %Overview which neuron responds to which picture significantly

% %Histogram
h=bar(significantperPic)
xlabel('Picture number');
ylabel('Frequency');
title('Histogram of significant responses per picture');
xlim([0 41]);

%Define variables
predrows = 1:16;
mouserows = 17:32;
sc = 2:2:16;
uns = 1:2:15;

% Splice only from after 0, because that is what we are interested in
pred = after0(:, predrows);
mice = after0(:, mouserows);
% Calculate the average
predsc2keep = mean(pred(:, sc),2);
predunsc2keep = mean(pred(:, uns),2);
miceunsc2keep = mean(mice(:, sc),2);
micesc2keep = mean(mice(:, uns),2);

%ttest
[h,p1]=ttest2(predsc2keep, predunsc2keep); %significance between predator scrambled and unscrambled
[h,p2]=ttest2(miceunsc2keep, predunsc2keep); %significance between mice and predator unscrambled
[h,p3]=ttest2(miceunsc2keep, micesc2keep); %sign. between mice unscrambled and scrambled

%Divide mice picture scrambled and intact
faceSPics = 18:2:32
faceNPics = 17:2:31

%% For every neuron seperate
truena0 = sum (na0 ==1); %Define how many time points are measured
ABArightneurons = find(ABAroi.region ==2); %Define which regions you want to include
faceSPics = 18:2:32;
faceNPics = 17:2:31;

% Intact pictures
resultNP = zeros(length(faceNPics),length(ABArightneurons), truena0, ImgsRepeats);
for j = 1: length(ABArightneurons)
    for i = 1:length(faceNPics)
        [resultNP(i,j,:,:)] = Res.CaSigCorrected(na0, Log== faceNPics(i), ABArightneurons(1,j));
    end
end
resultsNPdims=size(resultNP);

% Scrambled pictures
resultSP = zeros(length(faceSPics),length(ABArightneurons), truena0, ImgsRepeats);

for j = 1: length(ABArightneurons)
    for i = 1:length(faceSPics)

        [resultSP(i,j,:,:)] = Res.CaSigCorrected(na0, Log== faceSPics(i), ABArightneurons(1,j));
    end
end
resultsSPdims=size(resultSP);

% Averaging over time
resultNPmean = squeeze(mean(resultNP,3));
% dimension = size(resultNPmean);
resultSPmean = squeeze(mean(resultSP,3));

%Average over picture
resultNPmeanP = squeeze(mean(resultNPmean))
resultSPmeanP = squeeze(mean(resultSPmean))

% Calculating the p values
pval = zeros(length(faceSPics), length(ABArightneurons)); % p value for every neuron selected for only the mouse pictures.
for j = 1:length(ABArightneurons)
    for i = 1:length(faceNPics)
        [h, pval(i,j)] = ttest2(resultSPmean(i,j,:), resultNPmean(i,j,:));
    end
end

%Significance
significantmousepics = pval<=0.05;
significantmousepicsperPicture = sum(significantmousepics,2);
significantmousepicsamountperNeuron = sum(significantmousepics);
imagesc(significantmousepics);

%% Plot lines to see if the scrambled or non-scrambled are higher
allmousePics = 17:32;
significantPictureNP = faceNPics(significantmousepicsperPicture ~= 0); % Picture numbers in raw data of which there is a significant effect of scrambled versus nonscrambled
significantPictureSP = faceSPics(significantmousepicsperPicture ~= 0);
significantNeuron = ABArightneurons(significantmousepicsamountperNeuron ~= 0); % Neurons which are significant
significantplotten = zeros(length(significantPictureSP), length(significantNeuron));
intactbigger = 0;
scatterdbigger = 0;

for j = 1:length(significantNeuron)
    for i = 1: length(significantPictureSP)
        whichpicture = find(faceNPics == significantPictureNP(i)); %Define which picture needs to be taken
        whichneuron = find(ABArightneurons == significantNeuron(j)); %Define which neuron needs to be looked at
        if pval(whichpicture, whichneuron) <= 0.05 %Control if in that spot, that neuron and that picture if the pvalue is below <0.05
            figure
            plot(mean(Res.CaSigCorrected(na0, Log == significantPictureSP(i), significantNeuron(j)),2))
            scatterd = mean(Res.CaSigCorrected(na0, Log == significantPictureSP(i), significantNeuron(j)),2);
            legend ({'Scrambled' ,'Unscrambled'})%Plot the mean data of the scatterd picture
            both = [significantNeuron(j),significantPictureNP(i)];
            title('Neuron and picture' , both) %Add in the title the neuron and picture number so it is clear when you only have the figure
            hold on
            plot(mean(Res.CaSigCorrected(na0, Log == significantPictureNP(i), significantNeuron(j)),2)) %Plot the mean data of the intact picture
            intact = mean(Res.CaSigCorrected(na0, Log == significantPictureNP(i), significantNeuron(j)),2);
            hold off
            if mean(intact) > mean(scatterd)
                intactbigger = intactbigger +1;
            elseif mean(scatterd) > mean(intact)
                scatterdbigger = scatterdbigger +1;
            end
            if intactbigger > scatterdbigger
                difference = intactbigger - scatterdbigger;
                allvar = [intactbigger scatterdbigger difference];
                fprintf('There were more intact pictures (%d) which evoked a bigger response than there were scatterd pictures (%d), with a difference of %d\n', allvar )
            else scatterdbigger > intactbigger
                difference = scatterdbigger - intactbigger;
                allvar = [intactbigger scatterdbigger difference];
                fprintf('There were less intact pictures (%d) which evoked a bigger response than there were scatterd pictures (%d), with a difference of %d\n', allvar )
            end
        end
    end
end

%% Baseline data
% Per picture, calculate the difference between baseline and reaction to
% picture for every repeat. Average this for the repeats so there is one
% value for every neuron every picture. Then treat the neurons as repeats


diffBRrep = zeros(nImgs, length(ABArightneurons), ImgsRepeats);

for i = 1:(nImgs) %For every picture
    iamgespot = find(Log == i); %Find the n-th picture
    for r = 1:ImgsRepeats
        for j = 1:ABArightneurons %chosenROIs
            baseline = Res.CaSigCorrected(voor0, iamgespot(r), j);
            baseline = mean(baseline); % Average the baseline
            response = Res.CaSigCorrected(na0, iamgespot(r), j);
            response = mean(response); % Average the response
            diffBRrep(i,j,r) = response - baseline; % Calculate difference between baseline and response and store for every picture, neuron and repeat
        end
    end
end

diffBR = mean(diffBRrep,3); % Average the difference value over the repeats

% %Divide mice picture scrambled and intact
% faceSPics = 18:2:32;
% faceNPics = 17:2:31;

diffNP = diffBR(faceNPics);
diffSP = diffBR(faceSPics);

% ttest
[h ,pdiff, ci ,stats] = ttest2(diffNP, diffSP) %Ttest for the mouse pictures, neurons as repeats.
if pdiff <0.05
    fprintf ('There is a significance difference between the scrambled and non scrambled pictures')
else 
    fprintf('There is no significance difference')
end


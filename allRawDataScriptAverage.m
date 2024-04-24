%% Load in the data separately
% Made by Geke van Tongeren april 2024

folder = 'X:\Huub\SarinahMaryamMatlabTest\rawData';
dirlist = dir(fullfile(folder, '*.mat'));
nFile = length(dirlist);
names = {dirlist.name};
logpart = 'log.mat';

%% Define the variables
load('Drie_20231228_006_log.mat'); %Loading one file is needed due to defining of the variables
load('Drie_20231228_006_normcorr_SPSIG_Res.mat');

nrois = size(Res.CaSigCorrected, 3); %Number of neurons
chosenROIs = 1:nrois; %making a vector of all neurons
nImgs = length(unique(Log)); %Unique images, how many stimuli were shown
ImgsRepeats = sum (Log == Log(1)); %How many times was the first (and thus any picture) repeated in total
p_values = zeros(nImgs,nrois); % making a empty structure of zeros to put in the pvalues later
v = Res.ax>=-1 & Res.ax <=1; %Defining time frame of interest.
voor0 = Res.ax >= -0.99 & Res.ax<=0; %Before stimulus was shown
na0 = Res.ax>0 & Res.ax <=1; %After stimulus was shown
truena0 = sum (na0 ==1); %Define how many time points are measured
faceSPics = 18:2:32; %Only mouse pics
faceNPics = 17:2:31; %Only mouse pics
% faceSPics = 2:2:32; %All face pictures, scrambled
% faceNPics = 1:2:31; %All face pictures, Nonscrambled/intact

%% Load in every Res file and log file which belong together
for i = 1:nFile %For the whole length of the folder see if the filename contains Res
    if  contains(names{i} , 'normcorr_SPSIG_Res') % If the file contains the later part of norm.. then load it.
        load(names{i})
        [filename, matches] = strsplit(names{i}, '_'); %Split the file name you just loaded into different parts
        logfile = append(filename(1), '_', filename(2) ,'_', filename(3), '_', logpart); %Attach the filename together and load the log file of the same trial
        load(logfile{1}) %Now both the mat and log file are loaded of the same trial
        trialname = append(filename(1), '_', filename(2) ,'_', filename(3));

        % Store every region of one file in different matrixes
        % For every trial this must be run
        whichregions = unique(ABAroi.region);%Which regions are within the trial present
        amountofregions = length(whichregions); %How many different regions are within the trial present
        for r = 1 :amountofregions
            ABArightneurons = find(ABAroi.region ==whichregions(r)); %Find the r-th region and select only the neurons that are in this region
            amofthisregion = sum(ABAroi.region == whichregions(r)); %How many neurons are in this r-th region

            % Intact pictures
            resultNP = zeros(length(faceNPics),length(ABArightneurons), truena0, ImgsRepeats);
            for j = 1: length(ABArightneurons)
                for i = 1:length(faceNPics)
                    [resultNP(i,j,:,:)] = Res.CaSigCorrected(na0, Log== faceNPics(i), ABArightneurons(1,j));
                end
            end

            % Scrambled pictures
            resultSP = zeros(length(faceSPics),length(ABArightneurons), truena0, ImgsRepeats);

            for j = 1: length(ABArightneurons)
                for i = 1:length(faceSPics)

                    [resultSP(i,j,:,:)] = Res.CaSigCorrected(na0, Log== faceSPics(i), ABArightneurons(1,j));
                end
            end

            resultNPmean = squeeze(mean(resultNP,3));
            resultNPmean = squeeze(mean(resultNPmean,1));
            resultNPmeanM = mean(resultNPmean,2);

            resultSPmean = squeeze(mean(resultSP,3));
            resultSPmean = squeeze(mean(resultSPmean,1));
            resultSPmeanM = mean(resultSPmean,2);

            mask = resultNPmeanM > resultSPmeanM;

            % Calculating the p values
            pval = zeros(1, length(ABArightneurons)); % p value for every neuron selected for only the mouse pictures.
            for j = 1:length(ABArightneurons)
                [h, pval(j)] = ttest2(resultSPmean(j,:), resultNPmean(j,:));
                number = whichregions(r);
            end

            pval(mask) = NaN;

            nameoftrial = trialname{1}; % Defining the trial as a name
            pvalfilename = sprintf('pval_%d_%s', number, nameoftrial); %Make an individual file based on number (based on r the region) and on the nameof the trial
            eval([pvalfilename '= pval;']);
        end
    end
end

% Make a pval file per regionregion =[];
region =[];
wholist = who;
listofpv = wholist(startsWith(wholist, 'pval_'));
pvnumber = sum(contains(who, 'pval_'));

for contain=0:7 %Define the region possibilities
    containment = sprintf('%d_',contain); %Define the name of the file based on the region
    loop = sum(contains(listofpv, containment)); %How many files are there that include this region
    list = wholist(contains(wholist, containment)); %What are the names of the files that include the region
    for l = 1:loop
        listpva = list{l};
        listpva = eval(listpva);
        region = [region, listpva];
    end
    regionname = sprintf('region_%s', containment); %Rename file based on region
    eval([regionname '= region;']);

    sig = region <= 0.05; %Make a figure of all significant results per region
    sNeuR = sum(sig);
    figure, imagesc(sig)
    title('Region', contain) %Figure per region

    amountofneuronstotal = length(region);
    spPicR = sum(sig,2);
    % spPicR = sum(spPicR)

    sNeuR_region = sprintf('TotaalSigNeu_%d', contain);
    eval([sNeuR_region '=spPicR']);

    NeuTot = sprintf('NeuTot_%d',contain);
    eval([NeuTot '=amountofneuronstotal']);

    Fraction = (spPicR / amountofneuronstotal)*100;
    FractionSig = sprintf('FractionSig_%d', contain);
    eval([FractionSig '=Fraction']);

    region =[] ;
end

if FractionSig_1 > FractionSig_7
    fprintf('The fraction of significant neurons in region 1 (%.3f) is bigger then in region 7 (%.3f)\n', FractionSig_1, FractionSig_7);
elseif FractionSig_1 == FractionSig_7
    fprintf('The fraction of significant neurons in region 1 (%.3f) is equal then in region 7 (%.3f)\n', FractionSig_1, FractionSig_7);
else
    fprintf('The fraction of significant neurons in region 1 (%.3f) is smaller then in region 7 (%.3f) \n', FractionSig_1, FractionSig_7);
end

%% Baseline the data
for F = 1:nFile %For the whole length of the folder see if the filename contains Res
    if  contains(names{F} , 'normcorr_SPSIG_Res') % If the file contains the later part of norm.. then load it.
        load(names{F})
        [filename, matches] = strsplit(names{F}, '_'); %Split the file name you just loaded into different parts
        logfile = append(filename(1), '_', filename(2) ,'_', filename(3), '_', logpart); %Attach the filename together and load the log file of the same trial
        load(logfile{1}) %Now both the mat and log file are loaded of the same trial
        trialname = append(filename(1), '_', filename(2) ,'_', filename(3));

        % Store every region of one file in different matrixes
        % For every trial this must be run
        whichregions = unique(ABAroi.region);%Which regions are within the trial present
        amountofregions = length(whichregions); %How many different regions are within the trial present
        for g = 1 :amountofregions
            ABArightneurons = find(ABAroi.region ==whichregions(g)); %Find the r-th region and select only the neurons that are in this region
            amofthisregion = sum(ABAroi.region == whichregions(g)); %How many neurons are in this r-th region

            diffBRrepN = zeros(length(faceNPics), length(ABArightneurons), ImgsRepeats);
            diffBRrepS = zeros(length(faceSPics), length(ABArightneurons), ImgsRepeats);

            % Intact pictures
            for i = 1:length(faceNPics) %For every picture
                iamgespot = find(Log == faceNPics(i)); %Find the n-th picture
                for r = 1:ImgsRepeats
                    for j = 1:length(ABArightneurons) %chosenROIs
                        baseline = Res.CaSigCorrected(voor0, iamgespot(r), j);
                        baseline = mean(baseline); % Average the baseline
                        response = Res.CaSigCorrected(na0, iamgespot(r), j);
                        response = max(response); % Find maximum of the response = amplitude
                        value = response - baseline;
                        diffBRrepN(i,j,r) = value; % Calculate difference between baseline and response and store for every mouse picture, neuron and repeat
                    end
                end
            end

            %Scrambled pictures
            for i = 1:length(faceSPics) %For every picture
                iamgespot = find(Log == faceSPics(i)); %Find the n-th picture
                for r = 1:ImgsRepeats
                    for j = 1:length(ABArightneurons) %chosenROIs
                        baseline = Res.CaSigCorrected(voor0, iamgespot(r), j);
                        baseline = mean(baseline); % Average the baseline
                        response = Res.CaSigCorrected(na0, iamgespot(r), j);
                        response = max(response); % Find maximum of the response = amplitude
                        value = response - baseline;
                        diffBRrepS(i,j,r) = value; % Calculate difference between baseline and response and store for every mouse picture, neuron and repeat
                    end
                end
            end

            diffBRrepN = squeeze(mean(diffBRrepN));
            diffBRrepS = squeeze(mean(diffBRrepS));

            % Significant difference
            pva = zeros(length(ABArightneurons),1); % p value for every neuron selected for only the mouse pictures.
            for lrn = 1:length(ABArightneurons)
                [h, pva(lrn)] = ttest2(diffBRrepN(lrn,:), diffBRrepS(lrn,:));
            end

            % Deleting all in which the mean of amplitude of the trials is in the
            % scrambled higher than intact
            meanN =mean(diffBRrepN,2);
            meanS = mean(diffBRrepS,2);
            mask = meanN < meanS;
            pva(mask) = NaN;
            % pva = pva.';
            regionnumber =  whichregions(g);
            
            %Name specificly to region and trial
            % Fractiontemp = sprintf('Fraction_Sig_%d', regionnumber);
            % eval([Fractiontemp '=Fraction']);
            pvaluestemp = sprintf(['Pvalues_%d_%s'], regionnumber,trialname{1});
            eval([pvaluestemp '=pva'])
        end
    end
end

region =[];
wholist = who;
listofpv = wholist(startsWith(wholist, 'Pvalues_'));
pvnumber = sum(contains(who, 'Pvalues_'));

for contain=0:7 %Define the region possibilities
    containment = sprintf('_%d_',contain); %Define the name of the file based on the region
    loop = sum(contains(listofpv, containment)); %How many files are there that include this region
    list = wholist(contains(wholist, containment)); %What are the names of the files that include the region
    for l = 1:loop
        listpva = list{l};
        listpva = eval(listpva);
        region = vertcat(region, listpva) ;
    end
    regionname = sprintf('region%s', containment); %Rename file based on region
    eval([regionname '= region;']);
  
    sig = region <= 0.05; %Make a figure of all significant results per region
    sNeuR = sum(sig);

    amountofneuronstotal = length(region);
    % spPicR = sum(sig,2);
    % spPicR = sum(spPicR)

    sNeuR_region = sprintf('TotaalSigNeu%d', contain);
    eval([sNeuR_region '=sNeuR']);

    NeuTot = sprintf('NeuTot%d',contain);
    eval([NeuTot '=amountofneuronstotal']);

    Fraction = (sNeuR/ amountofneuronstotal)*100;
    FractionSig = sprintf('FractionSig%d', contain);
    eval([FractionSig '=Fraction']);
    sig= sig.'
    figure, imagesc(sig)
    titlename =[contain, Fraction]
    title('Region', titlename) %Figure per region

    region =[]
end

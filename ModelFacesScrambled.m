%% Load in the data separately
% Made by Geke van Tongeren april 2024

folder = 'X:\Huub\SarinahMaryamMatlabTest\rawData';
dirlist = dir(fullfile(folder, '*.mat'));
nFile = length(dirlist);
names = {dirlist.name};
logpart = 'log.mat';
% Load in one file to be able to define the variables
load('Drie_20231228_006_log.mat'); %Loading one file is needed due to defining of the variables
load('Drie_20231228_006_normcorr_SPSIG_Res.mat');
%Define variables
nrois = size(Res.CaSigCorrected, 3); %Number of neurons
chosenROIs = 1:nrois; %making a vector of all neurons
nImgs = length(unique(Log)); %Unique images, how many stimuli were shown
ImgsRepeats = sum (Log == Log(1)); %How many times was the first (and thus any picture) repeated in total
na0 = Res.ax>0 & Res.ax <=1; %After stimulus was shown
voor0 =Res.ax>-1 & Res.ax <=0;
truena0 = sum (na0 ==1); %Define how many time points are measured
faceSPics = 18:2:32; %Only mouse pics
faceNPics = 17:2:31; %Only mouse pics
facePics = 17:32; %All mouse pics

%% Loading every file

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
            whatregion =whichregions(g);

            diffBRrepN = zeros(length(faceNPics), length(ABArightneurons), ImgsRepeats);
            diffBRrepS = zeros(length(faceSPics), length(ABArightneurons), ImgsRepeats);

            % For all mouse Npics
            for j = 1:length(ABArightneurons) %chosenROIs
                for i = 1:length(faceNPics) %For every picture
                    iamgespot = find(Log == faceNPics(i)); %Find the n-th picture
                    for r = 1:ImgsRepeats
                        baseline = Res.CaSigCorrected(voor0, iamgespot(r), j);
                        baseline = mean(baseline); % Average the baseline
                        response = Res.CaSigCorrected(na0, iamgespot(r), j);
                        response = max(response); % Find maximum of the response = amplitude
                        value = response - baseline;
                        diffBRrepN(i,j,r) = value; % Calculate difference between baseline and response and store for every mouse picture, neuron and repeat
                    end
                end
            end

            diffBRrepN = permute(diffBRrepN, [3,2,1]);
            diffBRrepN = reshape(diffBRrepN,240,length(ABArightneurons));

            % For all mouse Spics
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

            diffBRrepS = permute(diffBRrepS, [3,2,1]);
            diffBRrepS = reshape(diffBRrepS,240,length(ABArightneurons));

            % diffBRrep = [diffBRrepN; diffBRrepS];
            % difftemp = sprintf(['diffBRrep_%d_%s'], whatregion, trialname{1});
            % eval([difftemp '=diffBRrep']);

            difftemp = sprintf(['diffBRrepN_%d_%s'], whatregion, trialname{1});
            eval([difftemp '=diffBRrepN']);
            difftemp = sprintf(['diffBRrepS_%d_%s'], whatregion, trialname{1});
            eval([difftemp '=diffBRrepS']);
        end
    end
end

% Patching all the files together
region =[];
wholist = who;
listofdiff = wholist(startsWith(wholist, 'diffBRrep'));
listN = wholist(startsWith(wholist, 'diffBRrepN'));
listS = wholist(startsWith(wholist, 'diffBRrepS'));

for contain=0:7 %Define the region possibilities
    containment = sprintf('_%d_',contain); %Define the name of the file based on the region
    loopN = sum(contains(listN, containment)); %How many files are there that include this region
    loopS = sum(contains(listS, containment)); %How many files are there that include this region
    list = wholist(contains(wholist, containment)); %What are the names of the files that include the region
    for l = 1:loopN
        listpva = list{l};
        listpva = eval(listpva);
        region = [region, listpva] ;
    end
    regionname = sprintf('useregionN%s', containment); %Rename file based on region
    eval([regionname '= region;']);

    region=[];

    for l = 1:loopS
        listpva = list{l};
        listpva = eval(listpva);
        region = [region, listpva] ;
    end
    regionname = sprintf('useregionS%s', containment); %Rename file based on region
    eval([regionname '= region;'])

    region =[];
end

wholist= who;
listuse = wholist(startsWith(wholist, 'useregion'));
theregion =[];

for contain =0:7
    containment = sprintf('_%d_',contain); %Define the name of the file based on the region
    loop = sum(contains(listuse, containment)); %How many files are there that include this region
    list = listuse(contains(listuse, containment)); %What are the names of the files that include the region
    for y=1:loop
        listU = list{y};
        listU = eval(listU);
        theregion = [theregion; listU]; %Stick them together, below each other
    end

    theregion = zscore(theregion);

    tempregion = sprintf('theregion%d', contain); %Rename the structure for every individual region
    eval([tempregion '=theregion'])
    theregion =[];
end

%% 
clear theregion
wholist= who;
listofregions =wholist(startsWith(wholist, 'theregion'));

% Making a model and training it
for contain=0:7
    containt = sprintf('%d', contain);
    listofcertainregion = listofregions(contains(listofregions, containt));
    regionX = listofcertainregion{1};
    regionX = eval(regionX);
    if isempty(regionX)
        sprintf('This region %d is empty, no average accuracy can be calculated', contain)
    else

        % %Second possibility for the model
        % accuracies=zeros(10,1)
        % 
        % for a=1:10
        % 
        %     hpartition = cvpartition(n, 'Holdout', 0.2)
        %     idxTrain = training(hpartition)
        % 
        %     tb1Train = regionX(idxTrain,:)
        %     tb1Labels = Labels(idxTrain,:)
        %     LabelsTest = Labels(idxNew)
        %     idxNew = test(hpartition)
        %     tb1New = regionX(idxNew,:)
        % 
        %     Md = fitcsvm(tb1Train, tb1Labels)
        %     cvMdl = crossval(Md); % Performs stratified 10-fold cross-validation
        %     cvtrainError = kfoldLoss(cvMdl)
        %     Accuracy = 1 - cvtrainError
        % 
        %     predictedLabels = predict(Md, tb1New);
        %     accuracies(a) = sum(predictedLabels == LabelsTest) / numel(LabelsTest)
        %     accuraciesm = mean(accuracies)
        %     AvgAccuTemp = sprintf('AvgAccu_Region_%d', contain); %Make an individual variable for every region and the average accuracy for the model
        %     eval([AvgAccuTemp '= accuraciesm'])
        % end

        % %Old model
        % % Split matrix into two
        % NShalf = regionX(1:240, :);
        % NShalf = NShalf
        % Shalf = regionX(241:end, :);
        % % regionX = [NShalf; Shalf]

        ncv = 100;
        accuracies = zeros(ncv,1);

        for a = 1:ncv % How many times do you want to train the model
            % Randomly select 20% of rows from each part for testing
            halfSize = size(NShalf,1);
            NumberOfRows = round(0.2 * halfSize); % 20% rows
            NShalf_test = randperm(halfSize, NumberOfRows);
            LabelNStest = ones(1,length(NShalf_test)); %Defining labels
            Shalf_test = randperm(halfSize, NumberOfRows);
            LabelStest = ones(1,length(Shalf_test))*2;

            % Randomly select the remaining 80% of rows for training
            NShalf_train = setdiff(1:size(NShalf, 1), NShalf_test); %Remaining amount for Training
            LabelNSTrain = ones(1,length(NShalf_train)); %Defining labels
            Shalf_train = setdiff(1:size(Shalf, 1), Shalf_test);
            LabelSTrain = ones(1,length(Shalf_train))*2;

            % Combine Labels
            LabelsTest = [LabelNStest, LabelStest]';
            % LabelsTest = LabelsTest';
            LabelsTrain = [LabelNSTrain, LabelSTrain]';
            % LabelsTrain = LabelsTrain';
            % LabelN = ones(1, length(regionX)/2)'
            % LabelS = (ones(1, length(regionX)/2)*2)'
            % Labels = [LabelN; LabelS]

            % Combine training and test
            Train = [NShalf_train, (Shalf_train + 240)]; %index the rows
            Test = [NShalf_test, (Shalf_test + 240)]; %index the rows

            % Produce train and test data set
            TrainData = regionX(Train, :);
            TestData = regionX(Test, :);

            % Model
            classifier = fitcdiscr(TrainData, LabelsTrain, 'DiscrimType','pseudoQuadratic');
            predictedLabels = predict(classifier, TestData);
            amountof1(a) = sum(ismember(predictedLabels, 1));
            accuracies(a) = sum(predictedLabels == LabelsTest) / numel(LabelsTest);
        end

        AvgAccu = mean(accuracies);
        AvgAccuTemp = sprintf('AvgAccu_R%d', contain); %Make an individual variable for every region and the average accuracy for the model
        eval([AvgAccuTemp '= AvgAccu';])

        amountof1 =mean(amountof1);
        amountof1temp = sprintf('Amount1_%d', contain);
        eval([amountof1temp '=amountof1'])
    end
end
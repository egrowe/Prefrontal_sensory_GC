%% Visualise and decode NP coherence, power and Granger Causality results

function supp_s08_decode_Oz(arrayInput)

% This script takes the scalp-EEG waveform from electrode Oz and runs
% classificaton analysis using 70/30 (train/test) split and 10-fold cross
% validation (CV).

% The output decoding accuracy represents the ability of the classifying to
% distinguish between a FACE and RANDOM trial

% The script runs for each participant and experimental phase

% This script requires libSVM toolbox to run

%Add the libSVM decoding library to the path
cd(['../decode_libSVM/libsvm/matlab'])
make

addpath(genpath(['../decode_libSVM/libsvm/matlab']))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countA = 1; countB = 1; countC = 1;

sprintf(['This is the current array input number:' num2str(arrayInput)])

%Define parameters
load('PLIST.mat'); %participant list
PName = PList(arrayInput,:);
sprintf(['This P selected:' num2str(PName)])
thisPhase = 'Phase1';
PhIdx = 1; %PhIdx 3 = Ph3, 2 = Ph 2 and 1 = Ph 1

trainProp = 0.7; %proportion of trials to use for training
timeWindow = [0,500]; %timewindow of interest

useThisName = ['_AT_OZ_']; %use savename relevant to above coords
maxFreq = 55; %maximum freq bin to examine (from 0 to maxFreq)
nCVReps = 10; %number of CV repetitions
zVal = 1; %if set to one then zscore the input test/train data

%Setup some counters (so we only take power values ONCE)
count = 1; cA = 1; cB = 1; cA_PFC = 1; cB_PFC = 1;
cC = 1; cD = 1; cC_PFC = 1; cD_PFC = 1;

for thisCVRep = 1:nCVReps

    % Setup the load names and save names
    loadnameFaces = ['NEW_ALL_OZ_faces_TrialData_' num2str(PName) '_All_Phases.mat']; %all data from Oz for this participant
    loadnameRandom = ['NEW_ALL_OZ_random_TrialData_' num2str(PName) '_All_Phases.mat']

    savename_Results = ['results_OZ_decode_' (PName) '_' (thisPhase) '_CVmethod_10REPS_' ...
        num2str(useThisName) '_0to500ms_thisRep' num2str(thisCVRep) '.mat'];

    % ASSIGN DECODING PARAMETERS (libSVM)
    nReps = 1; % number of repetitions of decoding steps
    costRange = [-10:5:10]; %cost range = 1 value (%-20:5:15)
    zscore = 1; % zscore decoding data (1 = 'yes', 0 = 'no')

    %% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load(loadnameFaces,'allTrial_Data_faces')
    load(loadnameRandom,'allTrial_Data_random')

    %Extract data from structure
    extractFaceData = allTrial_Data_faces(PhIdx,:);
    extractRandomData = allTrial_Data_random(PhIdx,:);

    %Place into matrix
    for hh = 1:length(extractFaceData)
        if isempty(extractFaceData{1,hh})
            continue
        end
        faceData_thisPhase(hh,:) = double(extractFaceData{1,hh});
    end

    for gg = 1:length(extractRandomData)
        if isempty(extractRandomData{1,gg})
            continue
        end
        randomData_thisPhase(gg,:) = double(extractRandomData{1,gg});
    end

    %Randomise the random number generator each loop
    rng('shuffle')

    %Randomly pick nTrials for RANDOM that equal nTrials for FACE
    if size(faceData_thisPhase,1) ~= size(randomData_thisPhase,1)
        randomData_thisPhase_thisRep = randomData_thisPhase(randperm(size(randomData_thisPhase,1),size(faceData_thisPhase,1)),:);
    end

    %Set trials in test/train
    trainSplit = ceil(size(faceData_thisPhase,1)*trainProp); %number of trials for 70%

    %Partition the FACE data at NODE X and Y randomly (randomise trial order)
    randomiseTrials_Face = randperm(size(faceData_thisPhase,1));  %select trial order randomly
    Face_train = faceData_thisPhase(randomiseTrials_Face(1:trainSplit),:); % select 70% train trials: Obs x Features
    Face_test = faceData_thisPhase(randomiseTrials_Face(trainSplit+1:end),:); % select 30% test trials: Obs x Features

    %Partition the RANDOM data NODE X and Y randomly (randmomise trial order)
    randomiseTrials_Rand = randperm(size(randomData_thisPhase_thisRep,1));  %select trial order randomly
    Random_train = randomData_thisPhase_thisRep(randomiseTrials_Rand(1:trainSplit),:); % select 70% train trials: Obs x Features
    Random_test = randomData_thisPhase_thisRep(randomiseTrials_Rand(trainSplit+1:end),:); % select 30% test trials: Obs x Features

    % NORMALISE THE TRAINING AND TEST DATA SEPARATELY (demean/divide stdev)
    % APPLY DEMEAN PROCESS (NEW VERSION)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TRAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_ERP_train = [Face_train; Random_train];
    all_ERP_train = all_ERP_train - mean(all_ERP_train); %demean
    ERPforGC_train = all_ERP_train./ std(all_ERP_train); % std across trials

    Face_train = ERPforGC_train(1:size(Face_train,1),:);
    Random_train = ERPforGC_train(size(Face_train,1)+1:end,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_ERP_test = [Face_test; Random_test];
    all_ERP_test = all_ERP_test - mean(all_ERP_test); %demean
    ERPforGC_test = all_ERP_test./ std(all_ERP_test); % std across trials

    Face_test = ERPforGC_test(1:size(Face_test,1),:);
    Random_test = ERPforGC_test(size(Face_test,1)+1:end,:);

    %% NOW GATHER FOR DECODING
    % %Gather GRANGER for BOTH hemispheres
    all_Train_faces = Face_train;
    all_Test_faces = Face_test;

    %% RANDOM: SQUEEZE DATA INTO REQUIRED FORMAT
    all_Train_rand = Random_train;
    all_Test_rand = Random_test;

    all_FACES_allReps(thisCVRep,:,:,:) = [all_Train_faces;all_Test_faces];
    all_RANDOM_allReps(thisCVRep,:,:,:) = [all_Train_rand;all_Test_rand];

    allTrain_data_FINAL = [all_Train_faces; all_Train_rand];
    allTest_data_FINAL = [all_Test_faces; all_Test_rand];

    %% NORMALISE ONCE ALL DATA FOR THIS REP IS EXTRACTED
    % Gather all LABELS for TRAIN and TEST sets
    classLabel_train = [repmat(1, size(allTrain_data_FINAL,1)/2,1); repmat(2, size(allTrain_data_FINAL,1)/2,1)];
    classLabel_test = [repmat(1, size(allTest_data_FINAL,1)/2,1); repmat(2, size(allTest_data_FINAL,1)/2,1)];

    if zVal == 1
        % Normalise training data
        means = mean(allTrain_data_FINAL, 1);
        stds = std(allTrain_data_FINAL, [], 1);
        means_mat = repmat(means, [size(allTrain_data_FINAL, 1), 1]);
        stds_mat = repmat(stds, [size(allTrain_data_FINAL, 1), 1]);
        allTrain_data_FINAL = (allTrain_data_FINAL - means_mat) ./ stds_mat;

        % Normalise validation data (using same parameters as for
        % training data)
        means_mat = repmat(means, [size(allTest_data_FINAL, 1), 1]);
        stds_mat = repmat(stds, [size(allTest_data_FINAL, 1), 1]);
        allTest_data_FINAL = (allTest_data_FINAL - means_mat) ./ stds_mat;
    end


    %% Optimise the cost parameter through grid search
    % bestcv = [];
    thisChanCV = 0;
    for costCounter = 1:length(costRange)
        log2c = 2^costRange(costCounter);

        %Setup the model parameters and define data + labels
        cmdCV = ['-s 0 -t 0 -c ' num2str(log2c) ' -v 10']; %
        labels = classLabel_train;
        data = sparse(allTrain_data_FINAL);

        % 10-fold CROSS VALIDATION TO TUNE COST PARAMETER
        cv = svmtrain(labels, data , cmdCV); %CV accuracy
        cv_acc(costCounter) = cv; %save the CV accuracy

        % Determine the best cost parameter (i.e. highest CV accuracy)
        if (cv > thisChanCV)
            thisChanCV = cv; bestcv = costRange(costCounter);
        end
    end
    allBestCostParams = bestcv; %Save each of the best cost parameters (per repetition of outer loop)
    useCostFinal = 2^allBestCostParams; %choose best C parameter

    %% TRAIN THE CLASSIFIER
    cmdModel = ['-s 0 -t 0 -c ' num2str(useCostFinal) ' -q']; %

    %TRAIN THE CLASSIFIER
    labels = classLabel_train;
    data = sparse(allTrain_data_FINAL);
    model = svmtrain(labels, data , cmdModel);

    %% TEST THE CLASSIFIER
    [predicted_label, accuracy, decision_values] = svmpredict(classLabel_test,sparse(allTest_data_FINAL),model);
    decodability(thisCVRep) = accuracy(1,:)

    clear classLabel_test allTest_data_FINAL model allTrain_data_FINAL all_Train_faces ...
        all_Test_faces all_Train_random all_Test_random

end

%MEAN ACCURACY
mean_accuracy =  mean(decodability)
std_accuracy = std(decodability)

%% SAVE DATA AND RESULTS
save(savename_Results)

end

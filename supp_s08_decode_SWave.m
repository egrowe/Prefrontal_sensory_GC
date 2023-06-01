%% Visualise and decode NP coherence, power and Granger Causality results

function supp_s08_decode_SWave(arrayInput)

% This script uses SVM classification (10-fold cross-validation) to
% classify whether the source-level waveforms contain information about the
% current stimulus (i.e., is it a face or random trial?)

%The ROIs are all PFC and sensory ROIs


%Add the libSVM decoding library to the path
cd(['../decode_libSVM/libsvm/matlab'])
make

addpath(genpath(['../decode_libSVM/libsvm/matlab']))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT: Settings and filenames for extracting ROI data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countA = 1; countB = 1; countC = 1;

sprintf(['This is the current array input number:' num2str(arrayInput)])

% SELECT THE TWO NODES TO DECODE
possCoords = {'L_SP_PFC','R_SP_PFC','L_Inf_PFC','R_Inf_PFC',...
    'L_FFA','R_FFA','L_Occ','R_Occ'}; % examine all of the ROIs individually

%Define parameters
load('PLIST.mat')
PName = PList(arrayInput,:);
sprintf(['This P selected:' num2str(PName)])
thisPhase = 'Phase1'; %examine this phase

trainProp = 0.7; %this training proportion of trials
timeWindow = [0,500]; %this timewindow

maxFreq = 55; %maximum freq bin to examine (from 0 to maxFreq)
nCVReps = 10; %number of CV repetitions
zVal = 1; %set to 1 to zscore the data

%Setup some counters (so we only take power values ONCE)
count = 1; cA = 1; cB = 1; cA_PFC = 1; cB_PFC = 1;
cC = 1; cD = 1; cC_PFC = 1; cD_PFC = 1;

for coordIdx = 1:length(possCoords)

        ROIcount = 1;

    for thisCVRep = 1:nCVReps


        useCoord = possCoords{coordIdx};

        loadname = ['CoordSorted_rmline_' (PName) '_allPhases_by_TopCoordPh3_BC_NoHP_EBRem_' ...
            num2str(useCoord) '_by_TstatOverTime_0_to_500ms.mat'];


        savename_Results = ['results_SWAVE_rmline_' (PName) '_' (thisPhase) '_CVmethod_10REPS_MultiROI_' ...
            num2str(useCoord) '_0to500ms_thisRep' num2str(thisCVRep) '.mat'];

        % ASSIGN DECODING PARAMETERS (libSVM)
        nReps = 1; % number of repetitions of decoding steps
        costRange = [-10:5:10]; %cost range = 1 value (%-20:5:15)
        zscore = 1; % zscore decoding data (1 = 'yes', 0 = 'no')

        %% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        load([loadname]); %load face data
        % SPLIT THE DATA INTO 70/30 for TRAIN/TEST
        allFACE_data_Node = faceData_thisCoord_Phase1;
        allRANDOM_data_Node = randData_thisCoord_Phase1; % CHANGE THIS DEPENDING ON PHASE!!!

        %Randomly pick nTrials for RANDOM that equal nTrials for FACE
        if size(allFACE_data_Node,1) ~= size(allRANDOM_data_Node,1)
            allRANDOM_data_Node_thisRep = allRANDOM_data_Node(randperm(size(allRANDOM_data_Node,1),size(allFACE_data_Node,1)),:);
        end

        %Set trials in test/train
        trainSplit = ceil(size(allFACE_data_Node,1)*trainProp); %number of trials for 70%

        %Partition the FACE data at NODE X and Y randomly (randomise trial order)
        randomiseTrials_Face = randperm(size(allFACE_data_Node,1));  %select trial order randomly
        Face_train = allFACE_data_Node(randomiseTrials_Face(1:trainSplit),:); % select 70% train trials: Obs x Features
        Face_test = allFACE_data_Node(randomiseTrials_Face(trainSplit+1:end),:); % select 30% test trials: Obs x Features

        %Partition the RANDOM data NODE X and Y randomly (randmomise trial order)
        randomiseTrials_Rand = randperm(size(allRANDOM_data_Node_thisRep,1));  %select trial order randomly
        Random_train = allRANDOM_data_Node_thisRep(randomiseTrials_Rand(1:trainSplit),:); % select 70% train trials: Obs x Features
        Random_test = allRANDOM_data_Node_thisRep(randomiseTrials_Rand(trainSplit+1:end),:); % select 30% test trials: Obs x Features

        % %Gather GRANGER for BOTH hemispheres
        all_Train_SWave_faces(:,:,ROIcount) = Face_train;
        all_Test_SWave_faces(:,:,ROIcount) = Face_test;

        %% RANDOM: SQUEEZE DATA INTO REQUIRED FORMAT
        all_Train_SWave_rand(:,:,ROIcount) = Random_train;
        all_Test_SWave_rand(:,:,ROIcount) = Random_test;

        all_FACES_Sensory_to_PFC_allReps(thisCVRep,:,:,:) = [all_Train_SWave_faces;all_Test_SWave_faces];
        all_RANDOM_Sensory_to_PFC_allReps(thisCVRep,:,:,:) = [all_Train_SWave_rand;all_Test_SWave_rand];


        %% ORGANISE THE DATA ONCE ALL GATHERED
        all_Train_Face = reshape(all_Train_SWave_faces,[size(all_Train_SWave_faces,1),...
            size(all_Train_SWave_faces,2)*size(all_Train_SWave_faces,3)]);
        all_Test_Face = reshape(all_Test_SWave_faces,[size(all_Test_SWave_faces,1),...
            size(all_Test_SWave_faces,2)*size(all_Test_SWave_faces,3)]);

        all_Train_Rand = reshape(all_Train_SWave_rand,[size(all_Train_SWave_rand,1),...
            size(all_Train_SWave_rand,2)*size(all_Train_SWave_rand,3)]);
        all_Test_Rand = reshape(all_Test_SWave_rand,[size(all_Test_SWave_rand,1),...
            size(all_Test_SWave_rand,2)*size(all_Test_SWave_rand,3)]);

        allTrain_Sensory_to_PFC = [all_Train_Face; all_Train_Rand];
        allTest_Sensory_to_PFC = [all_Test_Face; all_Test_Rand];

        %% NORMALISE ONCE ALL DATA FOR THIS REP IS EXTRACTED
        % Gather all LABELS for TRAIN and TEST sets
        classLabel_toPFC_train = [repmat(1, size(allTrain_Sensory_to_PFC,1)/2,1); repmat(2, size(allTrain_Sensory_to_PFC,1)/2,1)];
        classLabel_toPFC_test = [repmat(1, size(allTest_Sensory_to_PFC,1)/2,1); repmat(2, size(allTest_Sensory_to_PFC,1)/2,1)];

        if zVal == 1
            % Normalise training data
            means = mean(allTrain_Sensory_to_PFC, 1);
            stds = std(allTrain_Sensory_to_PFC, [], 1);
            means_mat = repmat(means, [size(allTrain_Sensory_to_PFC, 1), 1]);
            stds_mat = repmat(stds, [size(allTrain_Sensory_to_PFC, 1), 1]);
            allTrain_Sensory_to_PFC = (allTrain_Sensory_to_PFC - means_mat) ./ stds_mat;

            % Normalise validation data (using same parameters as for
            % training data)
            means_mat = repmat(means, [size(allTest_Sensory_to_PFC, 1), 1]);
            stds_mat = repmat(stds, [size(allTest_Sensory_to_PFC, 1), 1]);
            allTest_Sensory_to_PFC = (allTest_Sensory_to_PFC - means_mat) ./ stds_mat;
        end


        %% Optimise the cost parameter through grid search
        % bestcv = [];
        thisChanCV = 0;
        for costCounter = 1:length(costRange)
            log2c = 2^costRange(costCounter);

            %Setup the model parameters and define data + labels
            cmdCV = ['-s 0 -t 0 -c ' num2str(log2c) ' -v 10']; %
            labels = classLabel_toPFC_train;
            data = sparse(allTrain_Sensory_to_PFC);

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
        labels = classLabel_toPFC_train;
        data = sparse(allTrain_Sensory_to_PFC);
        model = svmtrain(labels, data , cmdModel);

        %% TEST THE CLASSIFIER
        [predicted_label, accuracy, decision_values] = svmpredict(classLabel_toPFC_test,sparse(allTest_Sensory_to_PFC),model);
        decodability(thisCVRep) = accuracy(1,:)

        clear classLabel_toPFC_test allTest_Sensory_to_PFC model allTrain_Sensory_to_PFC all_Train_SWave_faces ...
            all_Test_SWave_faces all_Train_GC_TO_PFC_random all_Test_GC_TO_PFC_random ...
            all_Train_SWave_faces all_Test_SWave_faces all_Test_SWave_rand all_Train_SWave_rand ...
            allFACE_data_Node allRANDOM_data_Node faceData_thisCoord_Phase1 randData_thisCoord_Phase1

    end

    %MEAN ACCURACY
    mean_accuracy =  mean(decodability)
    std_accuracy = std(decodability)

    %% SAVE DATA AND RESULTS
    save(savename_Results)

end


% script to loop over the construction of LIF networks and get
% corresponding statistics
clear all
Ratio_All = 1.0:0.2:5.0;       % clustering parameters
repetitions = 10;              % number of repetitions to draw statistics

% Set ExcInhNet = 1 if want to simulate EI networks - MUST CHANGE 
% NetworkLIF.m weight matrix section also
ExcInhNet = 0;

% Creating matrices to store different spike rate variability metrics
SlamScoreRV1 = zeros(repetitions,length(Ratio_All));
SlamScoreRV2 = SlamScoreRV1;
SlamScoreRV3 = SlamScoreRV1;
SlamScoreRV4 = SlamScoreRV1;

% These matrcies are for the variability across time
SlamScore2_sorted = zeros(size(SlamScoreRV1));
SlamScore2_shuff = zeros(size(SlamScoreRV1));

% Creating matrices to store eigenvalues and eigenspectral gaps
gap = zeros(repetitions, length(Ratio_All));
lambda_all = cell(repetitions,length(Ratio_All));

% Creating matrices to store the angles between schurs vectors and
% eigenvecotrs
PCA_EIG_Ang = zeros(repetitions,length(Ratio_All));
PCA_SCH_Ang = zeros(repetitions,length(Ratio_All));
PCA_EIG_Ang_All = cell(repetitions,length(Ratio_All));
PCA_SCH_Ang_All = cell(repetitions,length(Ratio_All));

% Matrix to store the random seed values of the random number generator if
% interested in repeating any particular cases
random_seeds = cell(repetitions,length(Ratio_All));



%%%% If want to run any network EXCEPT the EI scripts %%%%
if (ExcInhNet ~= 1) 
    
    % Loop over all values of the ratio being studied. Also perform the desired
    % number of repetitions. Store all the values
    for kkk = 1:repetitions
        for ss = 1:length(Ratio_All);
            Ratio = Ratio_All(ss);
            random_seeds{kkk,ss} = rng;
            
            %%% IT IS IMPORTANT TO ADJUST THIS FILE IN THE WEIGHT MATRIX
            %%% SECTION TO ENSURE CREATING RIGHT WEIGHT MATRIX. IT SIMPLY
            %%% INVOLVES COMMENTING/UNCOMMENTING THE WRITE LINE
            NetworkLIF
            conversion_spike_data;
            rast = rast(1:EneuronNum,:);
            
            % Calculate spike rate variability metric across groups
            [SlamScoreRV1(kkk,ss), SlamScoreRV2(kkk,ss), SlamScoreRV3(kkk,ss), SlamScoreRV4(kkk,ss)] = ...
                measure_rate_variability_groups(rast,numClusters);
            
            % Calculate spike rate variability metric across time
            [SlamScore2_sorted(kkk,ss), SlamScore2_shuff(kkk,ss)] = ...
                measure_rate_variability_groups_across_time(rast,numClusters);
            
            % Determine eigenvalues and the eigengap
            [lambdas, WW] = get_eigenvalues_LIF(weightsEE,weightsIE,weightsEI,weightsII);
            lambda_all{kkk,ss} = lambdas;
            gap(kkk,ss)  = extract_eigenvalue_gap_real(lambdas);
            
            % Find the PCs of the raster plots and compare their alignment with
            % the Schur vectors by calculating the angle span of the subspace
            calculateFiringRate;
            [AngEig, AngSchur] = calculateSubspaceAngle(WW,numClusters,firingRate);
            PCA_EIG_Ang_All{kkk,ss} = AngEig;
            PCA_SCH_Ang_All{kkk,ss} = AngSchur;
            PCA_EIG_Ang(kkk,ss) = min(AngEig);
            PCA_SCH_Ang(kkk,ss) = min(AngSchur);
            
        end
    end

 %%%%   If want to run the EI scripts   %%%%
else          

    % Loop over all values of the ratio being studied. Also perform the desired
    % number of repetitions. Store all the values
    for kkk = 1:repetitions
        for ss = 1:length(Ratio_All);
            
            Ratio = Ratio_All(ss);
            random_seeds{kkk,ss} = rng;
            
            %%% IT IS IMPORTANT TO ADJUST THIS FILE IN THE WEIGHT MATRIX
            %%% SECTION TO ENSURE CREATING RIGHT WEIGHT MATRIX. IT SIMPLY
            %%% INVOLVES COMMENTING/UNCOMMENTING THE WRITE LINE
            NetworkLIF            
            conversion_spike_data;

            % Calculate spike rate variability metric across groups
            [SlamScoreRV1(kkk,ss), SlamScoreRV2(kkk,ss), SlamScoreRV3(kkk,ss), SlamScoreRV4(kkk,ss)] = ...
                measure_rate_variability_groups_EI(rast,numClusters);
            
            % Calculate spike rate variability metric across time
            [SlamScore2_sorted(kkk,ss), SlamScore2_shuff(kkk,ss)] = ...
                measure_rate_variability_groups_across_time_EI(rast,numClusters);

            % Determine eigenvalues and the eigengap
            [lambdas, WW] = get_eigenvalues_LIF(weightsEE,weightsIE,weightsEI,weightsII);
            lambda_all{kkk,ss} = lambdas;
            gap(kkk,ss)  = extract_eigenvalue_gap_EI_real(lambdas);
            
            % Find the PCs of the raster plots and compare their alignment with
            % the Schur vectors by calculating the angle span of the subspace
            calculateFiringRate;
            [AngEig, AngSchur] = calculateSubspaceAngle(WW,numClusters,firingRate);
            PCA_EIG_Ang_All{kkk,ss} = AngEig;
            PCA_SCH_Ang_All{kkk,ss} = AngSchur;
            PCA_EIG_Ang(kkk,ss) = min(AngEig);
            PCA_SCH_Ang(kkk,ss) = min(AngSchur);
            
            
            
        end
    end
    
end

% Store all the ratios used
Ratio = Ratio_All;

% Save the data variable of interest. Can save everything of course!
% However this is done to save space. Storing the random seeds allows us to
% regenerate anything if needed.
save('stat_clustered', 'EneuronNum', 'Etm', 'Etr', 'IneuronNum', 'Itm', 'Itr', 'Ratio', ...
    'SlamScoreRV1', 'SlamScoreRV2', 'SlamScoreRV3',...
    'SlamScoreRV4', 'gap', 'neuronNum', 'numClusters', 'tEnd',...
    'tStart', 'tStep', 't_EE', 't_EI', 't_IE', 't_II', 'wEE', 'wEI', 'wIE', ...
    'wII','lambda_all','random_seeds',  'PCA_SCH_Ang', ...
    'PCA_SCH_Ang_All', 'SlamScore2_sorted','SlamScore2_shuff')

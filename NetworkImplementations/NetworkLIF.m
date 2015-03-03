%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart = 0;
tEnd = 2000;                      %Simulation in milli-seconds
tStep = 0.1;                      %0.1 millisecond time step

time = [tStart:tStep:tEnd];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting Parameters For Neuron Number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mult = 2;

EneuronNum  = mult*800;                 %Number of excitatory neurons in the network
numClusters = mult*10;                  %Number of clusters
IneuronNum  = round(0.25*EneuronNum);   %Number of inhibitory neurons in the network
neuronNum   = EneuronNum + IneuronNum;  %Total number of neurons


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT MATRIX - uncomment/comment required lines for desired network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%                 Create REE or WRatio weight matrix            %%%%
createWeightMatrix_REE_or_WRatio

%%%%                 Crease small world weight matrix              %%%%
% createWeightMatrix_SmallWorld

%%%%                 Create hierarchcal weight matrix              %%%%
% createWeightMatrix_Hierarchical

%%%%                 Create scale free weight matrix               %%%%
% createWeightMatrix_ScaleFree

%%%%       Create excitatory-inhibitory grouping weight matrix     %%%%
% createWeightMatrix_ExcInh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the LIF Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vthres = 1;      %Threshold for both exc and inh neurons
Vreset = 0;      %Reset for both exc and inh neurons

%Excitatory neuron params
Etm = 15;        %Membrane Time Constant
Etr = 5;         %Refractory period

%Inhibitory neuron params
Itm = 10;        %Membrane Time Constant
Itr = 5;         %Refractory period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the synapses and neuronal conductances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_EE = 3;                                       %time constant for excitatory to excitatory synapses
t_IE = 3;                                       %time constant for excitatory to inhibitory synapses
t_EI = 2;                                       %time constant for inhibitory to excitatory synapses
t_II = 2;                                       %time constant for inhibitory to inhibitory synapses

gEE = zeros(1,EneuronNum);                      %conductance for excitatory to excitatory synapses
gIE = zeros(1,IneuronNum);                      %conductance for excitatory to inhibitory synapses
gEI = zeros(1,EneuronNum);                      %conductance for inhibitory to excitatory synapses
gII = zeros(1,IneuronNum);                      %conductance for inhibitory to inhibitory synapses

gBE = 1.1 + .1*rand(1,EneuronNum);
gBI = 1 + .05* rand(1,IneuronNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solving Network with LIF neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rast = zeros(neuronNum,(tEnd - tStart)/tStep + 1);     %Matrix storing spike times for raster plots
lastAP  = -50 * ones(1,neuronNum);                     %last action potential for refractor period calculation (just big number negative put)

memVol = rand(neuronNum,(tEnd - tStart)/tStep + 1);

for i =2:(tEnd - tStart)/tStep
    for j = 1:neuronNum
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CONNCECTIVITY CALCULATIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (j <= EneuronNum)
            gEE(j) = gEE(j) - gEE(j)*tStep/t_EE;
            gEI(j) = gEI(j) - gEI(j)*tStep/t_EI;
        else
            gIE(j-EneuronNum) = gIE(j-EneuronNum) - gIE(j-EneuronNum)*tStep/t_IE;
            gII(j-EneuronNum) = gII(j-EneuronNum) - gII(j-EneuronNum)*tStep/t_II;
        end
        
        if (rast(j,i-1) ~= 0 && j <= EneuronNum)       %If excitatory neuron fired
            gEE = gEE + weightsEE(:,j)';
            gIE = gIE + weightsIE(:,j)';
        end
        
        if (rast(j,i-1) ~= 0 && j > EneuronNum)
            gEI = gEI + weightsEI(:,j-EneuronNum)';
            gII = gII + weightsII(:,j-EneuronNum)';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %EXCITATORY NEURONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(j <= EneuronNum)
            
            gE= gEE(j);
            gI= gEI(j);
            
            
            v = memVol(j,i-1) + (tStep/Etm)*(-memVol(j,i-1) + gBE(j)) + tStep*(gE -gI);
            
            if ((lastAP(j) + Etr/tStep)>=i)   %Refractory Period
                v = Vreset;
            end
            
            if (v > Vthres)           %Fire if exceed threshold
                v = Vthres;
                lastAP(j) = i;
                rast(j,i) = j;
            end
            
            memVol(j,i) = v;
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %INHIBITORY NEURONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (j > EneuronNum)
            gE = gIE(j - EneuronNum);
            gI = gII(j - EneuronNum);                          %If want to add Inh to Inh
            
            
            v = memVol(j,i-1) + (tStep/Itm)*(-memVol(j,i-1) + gBI(j-EneuronNum)) + tStep*(gE -gI);
            
            if ((lastAP(j) + Itr/tStep)>=i)     %Refractory Period
                v = Vreset;
            end
            
            if (v > Vthres)                  %Fire if exceed threshold
                v = Vthres;
                lastAP(j) = i;
                rast(j,i) = j;
            end
            
            memVol(j,i) = v;
        end
        
    end
end
   
[lambdas, W]= get_eigenvalues_LIF(weightsEE,weightsIE,weightsEI,weightsII);
% figure; plot(real(lambdas),imag(lambdas),'.','color',[0 0.3 0], 'MarkerSize',5)
% set(gca,'FontSize',20)
% set(gcf,'Color',[1 1 1])
% plotRASTER

% save('network_type_parameterValue')


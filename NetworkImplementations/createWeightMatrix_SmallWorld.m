%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE WEIGHT MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REE = 3.5;                   %Ratio of pin/pout for small world network

%If run_sim_loop original file being excuted, it will have a Ratio
%parameter and the value is updated according the current value of Ratio.
%It is upto the user of this code to change the line inside the if
%condition to REE or WRatio depending on which case they are studying
if exist('Ratio')
    REE = Ratio;
end

f = 1/sqrt(mult);          %Factor to scale by synaptic weight parameters by network size

wEI     = f*.042;          %Average weight of inhibitroy to excitatory cells
wIE     = f*0.0105;        %Average weight of excitatory to inhibitory cells
wEE     = f*.022;          %Average weight of excitatory to excitatory cells
wII     = f*0.042;         %Average weight of inhibitory to inhibitory cells

weightsEI = random('binom',1,0.5,[EneuronNum,IneuronNum]);      %Weight matrix of inhibioty to excitatory LIF cells
weightsEI = wEI.* weightsEI;

weightsIE = random('binom',1,0.5,[IneuronNum, EneuronNum]);     %Weight matrix of excitatory to inhibitory cells
weightsIE = wIE.* weightsIE;

weightsII = random('binom',1,0.5,[IneuronNum, IneuronNum]);     %Weight matrix of inhibitory to inhibitory cells
weightsII = wII.* weightsII;

p_in = 0.2/(1/numClusters+(1-1/numClusters)/REE); %internal probability
p_out = p_in/REE;


% create k-nearest neighbour graph
neurons_per_group = floor(EneuronNum/numClusters);
EE_ring = zeros(EneuronNum);
for jj = 1:neurons_per_group/2;
    EE_ring = EE_ring + diag(ones(EneuronNum-jj,1),jj);
end
EE_ring(1,end) =1;
EE_ring = EE_ring + EE_ring';

% weight matrix outside k-nearest neighbour ring
weightsEE = random('binom',1,p_out,[EneuronNum, EneuronNum]);  
weightsEE(weightsEE == EE_ring) = 0; % remove connections inside ring

% make probabilistic connections inside k-nearset neighbours ring
weightsEE_ring = (rand(EneuronNum)<=p_in).*EE_ring;
weightsEE = wEE.* (weightsEE+weightsEE_ring);


%Ensure the diagonals are zero
weightsEE = weightsEE -diag(diag(weightsEE));
weightsII = weightsII -diag(diag(weightsII));
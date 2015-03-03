
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WEIGHT MATRIX HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%Weight matrix of excitatory to excitatory cells
weightsEE = full(scalefree(EneuronNum, 64)).*wEE.*random('binom',1,1,[EneuronNum, EneuronNum]);

%Ensure the diagonals are zero
weightsII = weightsII - diag(diag(weightsII));
weightsEE = weightsEE - diag(diag(weightsEE));


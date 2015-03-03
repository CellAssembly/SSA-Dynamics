%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE WEIGHT MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factorEI = 3;
factorIE = 3;

if exist('Ratio')
    factorEI = Ratio;
    factorIE = Ratio;
end

pfactorEI = 2;             %Change the connection ratio also
pfactorIE = 2;             %Change the connection ratio also

f = 1/sqrt(mult);          %Factor to scale by synaptic weight parameters by network size

wEI     =f*0.042;          %Average weight of inhibitroy to inhibitory excitatory LIF units (single compartment)
wEIsub = wEI/(1/numClusters+(1-1/numClusters)*factorEI);
wEI = wEIsub*factorEI;
pEI = .5;
pEIsub = pEI/(1/numClusters+(1-1/numClusters)*pfactorEI);
pEI = pEIsub*pfactorEI;


wIE     = f*0.0105;        %Average weight of excitatory to inhibitory cells
wIEsub = wIE/(1/numClusters+(1-1/numClusters)/factorIE);
wIE = wIEsub/factorIE;
pIE = .5;
pIEsub = pIE/(1/numClusters+(1-1/numClusters)/pfactorIE);
pIE = pIEsub/pfactorIE;


wEE     = f*0.022;         %Average weight of excitatory to excitatory cells (all of which go to soma)
pEE     = .2;
wII     = f*0.042;         %Average weight of inhibitory to inhibitory cells
pII =.5;

weightsEI = random('binom',1,pEI,[EneuronNum,IneuronNum]);          %Weight matrix of inhibioty to single compartment excitatory LIF units
weightsEI = wEI.* weightsEI;

weightsIE = random('binom',1,pIE,[IneuronNum, EneuronNum]);         %Weight matrix of excitatory to inhibitory cells
weightsIE = wIE.* weightsIE;

weightsII = random('binom',1,pII,[IneuronNum, IneuronNum]);         %Weight matrix of inhibitory to inhibitory cells
weightsII = wII.* weightsII;

weightsEE = random('binom',1,pEE,[EneuronNum, EneuronNum]);         %Weight matrix of excitatory to excitatory cells
weightsEE = wEE.* weightsEE;   

%Create the group weight matrices for Exc to Inh and update the total weight matrix
for i = 1:numClusters
    weightsIEsub = random('binom',1,pIEsub,[IneuronNum/numClusters, EneuronNum/numClusters]);
    weightsIEsub = wIEsub.* weightsIEsub;              
    kk = (i-1)*IneuronNum/numClusters+1:i*IneuronNum/numClusters;
    jj = (i-1)*EneuronNum/numClusters+1:i*EneuronNum/numClusters;
    weightsIE(kk,jj) = weightsIEsub;
end

%Create the group weight matrices for Inh to Exc and update the total weight matrix
for i = 1:numClusters
    weightsEIsub = random('binom',1,pEIsub,[EneuronNum/numClusters, IneuronNum/numClusters]);
    weightsEIsub = wEIsub.* weightsEIsub;              
    kk = (i-1)*IneuronNum/numClusters+1:i*IneuronNum/numClusters;
    jj = (i-1)*EneuronNum/numClusters+1:i*EneuronNum/numClusters;
    weightsEI(jj,kk) = weightsEIsub;
end

%Ensure the diagonals are zero
weightsEE = weightsEE -diag(diag(weightsEE));
weightsII = weightsII -diag(diag(weightsII));
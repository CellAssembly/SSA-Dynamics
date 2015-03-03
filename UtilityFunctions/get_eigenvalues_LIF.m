function [ lambdas, W] = get_eigenvalues_LIF(weightsEE,weightsIE,weightsEI,weightsII)
% compute the eigenvalues of the LIF network

W = [weightsEE, -weightsEI; weightsIE, -weightsII];

lambdas = eig(W);

end


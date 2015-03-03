function [AngEig, AngSchur] = calculateSubspaceAngle(weightsAll,numClusters, firingRate)
% Caluclate subspace angles of firing rates with schur decomposition and eigenvector
[V,D] = eig(weightsAll);        %Find eigenvalues of weight matrix - stored in V
[U,T] = schur(weightsAll,'complex');      %Find Schur vectors - stored in U

% sort vectors according to real value
dd = diag(D);
[~,index] = sort(real(dd),'descend');
V = V(:,index);

dd = diag(T);
[~,index] = sort(real(dd),'descend');
U = U(:,index);

%If you use an older version of matlab, use the function princomp
princComps = princomp(firingRate');   

AngEig   = mPrinAngles(V(:,1:numClusters-1),princComps(:,1:numClusters-1));
AngSchur = mPrinAngles(U(:,1:numClusters-1),princComps(:,1:numClusters-1));

end





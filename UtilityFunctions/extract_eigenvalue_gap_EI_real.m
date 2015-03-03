function [ gap ] = extract_eigenvalue_gap_EI_real(lambdas)
%extract_eigenvalue_gap -- get the eigenvalue gap from the spectrum of LIF

% remove lambdas with negative real part to compute gap only on positive
% side..
lambdas(real(lambdas)<0) = [];

% now sort according to real value
lambda = sort(real(lambdas));

%  differences of lambda
d_lambda = diff(lambda);
[~, index] = sort(d_lambda);

% find real values of the eigenvalues that are at the boundaries of 
% the (now) largest gap
l1 = lambda(index(end));
l2 = lambda(index(end)+1);

% actual gap
gap = l2 -l1;

end


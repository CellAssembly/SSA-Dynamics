function [ gap ] = extract_eigenvalue_gap_real(lambdas)
%extract_eigenvalue_gap -- get the eigenvalue gap from the spectrum of LIF

% sort lambdas according to magnitude
lambda = sort(lambdas);

% NOTE: case where there is no outlier is not covered
% remove largest "outlier" eigenvalues
if conj(lambda(end)) == lambda(end-1) % if complex conjugate pair
    lambdas = lambda(1:end-2);
elseif isreal(lambda(end))  %single real outlier eigenvalue
    lambdas = lambda(1:end-1); 
end

% now that outliers are removed, sort according to real value
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


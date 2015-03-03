function [Sig, Sig_shuff] = measure_rate_variability_groups_across_time(rast,num_groups)
% compute covariance scores between spike trains
if nargin < 2 
    num_groups = 10;
end
[num_neurons, rast_length] =   size(rast);


% 0.1 ms is stepsize, we want to have 100ms windows...
bin100ms = 100/.1;

% new number of bins and length of last bin
n_new_bins = floor(rast_length/bin100ms);
last_bin_length = -n_new_bins*bin100ms + rast_length;

if last_bin_length > 1
    % if it is just one then this just comes from thes simulation time..
    warning('time of recording does not match with 100ms bins')
end

% partition vector for coarse graining, last bin is .1 ms longer unless
% warning appears
cg_vec = [kron(1:n_new_bins,ones(1,bin100ms)), (n_new_bins)*ones(1,last_bin_length) ];    
Hcg = transformPartitionVectorToHMatrix(cg_vec);
% note rast encodes neuron numbers not num of spikes, afterwards Z encodes
% number of spikes in each bin
Z = (rast>0)*Hcg;


% size of groups and partition vector + matrix
grpsize = num_neurons/num_groups;
pvec = kron(1:num_groups,ones(1,grpsize));
H= transformPartitionVectorToHMatrix(pvec);

Zfr = Z/100e-3; % Zfr is now the firing rate of each neuron in a 100ms
Zgrouped = H'*Zfr/grpsize; % same at group level

Sig = mean(std(Zgrouped')); % std deviation along columns (groups), mean along rows (over time)



num_rep = 10;
Sig_shuff = zeros(1,num_rep);
for i = 1:num_rep
    Z = Z(randperm(num_neurons),:); % permute neuron labels

    Zfr = Z/100e-3; % Zfr is now the firing rate of each neuron in a 100ms
    Zgrouped = H'*Zfr/grpsize; % same at group level

    Sig_shuff(i) = mean(std(Zgrouped'));
end
Sig_shuff = mean(Sig_shuff);




end

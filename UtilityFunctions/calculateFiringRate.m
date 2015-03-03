
%% Find the firing rate

spksAnalyse = spksAll;  

%in seconds
tWindow = 0.25;      
%Total number of time data points
totalPoints = round(max(spksAnalyse(:,2)))/tWindow;           
%Store firing rate for each neuron in N x T matrix
firingRate = zeros(EneuronNum+IneuronNum,floor(totalPoints));        

for i = 1:neuronNum
    %Find spikes and spike times for neuron i
    neuronSpks = spksAnalyse(spksAnalyse(:,1) == i,:);            
    for k = 1: floor(totalPoints)
        temp = neuronSpks(:,2);
        temp(temp > k*tWindow)     = [];                  
        temp(temp < (k-1)*tWindow) = [];
        %Number of spikes in time window
        numSpikes = length(temp);                         
        firingRate(i,k) = numSpikes./tWindow;            
    end
end



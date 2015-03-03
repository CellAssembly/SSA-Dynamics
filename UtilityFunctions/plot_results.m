% script to extract and plot results

% load data
load('stat_clustered.mat')
color = 'k';
% load('Wmatrices.mat')

% %% do plots vs REE/WRatio
% dataset = 'W';
% dirname = 'figures_results/data/';
% 
% switch dataset
%     case '1000'
%         load([dirname 'EEclusteredREE_1000/stat_clustered_new.mat'])
%         color = 'b';
%         REE = 2:.2:5; % define REE range to produce figures
%     case '2000'
%         load([dirname 'EEclusteredREE_2000/stat_clustered_new.mat'])
%         color = 'k';
%         REE = 2:.2:5; % define REE range to produce figures
%     case '3000'
%         load([dirname 'EEclusteredREE_3000/stat_clustered_new.mat'])
%         color = 'r';
%         REE = 2:.2:5; % define REE range to produce figures
%     case '2000_32'
%         load([dirname 'EEclusteredREE_2000_32groups/stat_clustered_new.mat'])
%         color = 'g';
%         REE = 2:.2:5; % define REE range to produce figures
%     case '2000_40'
%         load([dirname 'EEclusteredREE_2000_40groups/stat_clustered_new.mat'])
%         color = 'm';
%         REE = 2:.2:5; % define REE range to produce figures
%     case 'W'
%         load([dirname 'EEclusteredWratio_2000/stat_clustered_new.mat'])
%         color = 'k';
%         REE = 1:.2:4; % define range for WRatio, but store in REE variable
%     case 'WEI'
%         load([dirname 'EIclusteredW_2000_full/stat_clustered_new.mat'])
%         color = 'k';
%         REE = 1:.2:5; % define range for WRatio, but store in REE variable
% end


% plot REE vs coefficient of variation plot
figure(3); hold all
plot(Ratio,SlamScoreRV1-SlamScoreRV3,[color '.'])
m =mean(SlamScoreRV1-SlamScoreRV3);
s =std(SlamScoreRV1-SlamScoreRV3);
jbfill(Ratio, m+s, m-s,color); hold all
plot(Ratio,m,color)
set(gca,'LineWidth',1.5)

%% second part: plots vs eigenvalue gap

figure(9); hold all
plot(gap,SlamScoreRV1-SlamScoreRV3,[color '.'])

%% gap vs REE
figure(13); hold all
plot(Ratio,gap, [color '.'])

%% For plotting subspace angle figures

load('stat_clustered.mat')

figure; jbfill(Ratio,mean(cos(PCA_SCH_Ang))+std(cos(PCA_SCH_Ang)),mean(cos(PCA_SCH_Ang))-std(cos(PCA_SCH_Ang)),'k','k');
hold on; plot(Ratio,(cos(PCA_SCH_Ang)),'k.');
hold on; plot(Ratio,mean(cos(PCA_SCH_Ang)),'k');
set(gca,'FontSize',15)

% Remember to change the x-axis label if you need to.
xlabel('Ratio'); ylabel('PCA-SCHUR Angle');
set(gcf,'Color',[1 1 1]); 
set(gca,'FontSize',15)
set(gca,'LineWidth',2.0)

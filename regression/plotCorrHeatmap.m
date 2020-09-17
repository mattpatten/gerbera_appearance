function plotCorrHeatmap(allcorr, header, saveName)

%% plot correlation matrix for attributes

header = rename_headers(header);

f1 = figure;
img = imagesc(allcorr,[-1 1]); 
set(img,'AlphaData',~isnan(allcorr));
%colormap('parula'); 
colormap(colormapFromRSA);
colorbar; 
set(gca,'XTick',1:length(header),'XTickLabels',header,'FontSize',8,'TickLabelInterpreter','None');
set(gca,'YTick',1:length(header),'YTickLabels',header,'FontSize',8,'TickLabelInterpreter','None');
title('Attribute correlations','FontSize',12,'Interpreter','None');
xtickangle(45);
set(gca,'linewidth',1);
axis square; 
box on;
set(gcf, 'Position', [1, 41, 1128, 917]); %set size on screen
saveas(f1,[saveName '.png']);
close;

end
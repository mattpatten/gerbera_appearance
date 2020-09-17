function plotRegression(mean_relatedness, errorbar_lower, errorbar_higher, sort_idx, header, title_text, saveName)

%% Draw final figure

nAttributes = length(header);

f = figure;
hold all;

bar(1:nAttributes,mean_relatedness(sort_idx),'FaceColor',[68/225 131/225 149/225],'EdgeColor','None');
errorbar(1:nAttributes,mean_relatedness(sort_idx),errorbar_lower(sort_idx), errorbar_higher(sort_idx), ...
    'LineStyle','None','LineWidth',1,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','None','Marker','None','MarkerSize',3,'CapSize',3); %parameters

%Make it look nice
title(title_text,'Interpreter','None');
set(gca,'XLim',[0 nAttributes+1],'XTick',1:nAttributes,'XTickLabels',header(sort_idx),'TickLabelInterpreter','None');
%set(gca,'YLim',[0 71],'YTick',[1 10:10:70],'YDir','reverse');
xtickangle(45);
%xlabel('Flower ID');
set(gca,'ylim',[-0.25 0.25],'ytick',-0.2:0.1:0.2);
%ylabel('Relatedness (Spearman corr)');
ylabel('Beta values');
set(gca,'linewidth',1);
box off; %axis square;

set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
saveas(f,[saveName '.png']);
close;

end
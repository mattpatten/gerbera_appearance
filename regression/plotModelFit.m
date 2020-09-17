function plotModelFit(observed_data,predicted_median,predicted_lower,predicted_upper,rsquare_val,analysis_string,imgDir,outputDir,labels,sort_by_pred,session)

nFlowers = size(predicted_median,1);

%observed_data = eval(['lm.Variables.' lm.VariableNames{end}]);
%predicted_median = lm.Fitted;
diff_data = abs(observed_data - predicted_median);

if sort_by_pred
    [~, sort_idx] = sort(predicted_median,'descend');
    sortMethod = 'pred';
else
    [~, sort_idx] = sort(observed_data,'descend');
    sortMethod = 'obs';
end

%% Find influential observations
if nFlowers==70
    influential_flowers = get_influential_observers(session, 1);
    idx_regular = ~ismember(1:size(observed_data,1),influential_flowers); %cut off necessary positions/ranks
    idx_influential = ~idx_regular;
end

%ceilfix = @(x)ceil(abs(x)).*sign(x);
%ylim_vals = ceilfix([min([observed_data; predicted_median]) max(([observed_data; predicted_median]))]);
ylim_vals = [-1.65 1.65];

f = figure;
hold all;

%horizontal line at zero
line([0 nFlowers+1],[0 0],'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',2);

%Bootstrapped error bar
errorbar(1:nFlowers,predicted_median(sort_idx),predicted_lower(sort_idx),predicted_upper(sort_idx),...
    'LineStyle','None','LineWidth',1.5,'Color',rgb('Salmon'),'Marker','None','Capsize',3);

%observed values
if nFlowers==70
    sorted_obs = observed_data(sort_idx);
    %regular flowers
    plot(find(idx_regular(sort_idx)),sorted_obs(idx_regular(sort_idx)),'LineStyle','None','LineWidth',1,'Color',rgb('blue'),'MarkerFaceColor',rgb('blue'),'MarkerEdgeColor',rgb('blue'),'Marker','o','MarkerSize',6);
    %influential
    plot(find(idx_influential(sort_idx)),sorted_obs(idx_influential(sort_idx)),'LineStyle','None','LineWidth',1,'Color',rgb('blue'),'MarkerFaceColor','None','MarkerEdgeColor',rgb('blue'),'Marker','o','MarkerSize',6);
    legend({'Mean rating','95% CI of model predictions','Observed values - Regular flowers','Observed values - Influential flowers'},'Location','NorthEast');
else
    plot(1:nFlowers,observed_data(sort_idx),'LineStyle','None','LineWidth',1,'Color',rgb('blue'),'MarkerFaceColor',rgb('blue'),'MarkerEdgeColor',rgb('blue'),'Marker','o','MarkerSize',6);
    legend({'Mean rating','95% CI of model predictions','Observed values'},'Location','NorthEast');
end
%plot(1:nFlowers,predicted_median(sort_idx),'LineStyle','None','LineWidth',1,'Color',rgb('red'),'Marker','x','MarkerSize',8);

% Put lines indicating magnitude of residual going up from the x-axis
%line([(1:nFlowers); (1:nFlowers)],[repmat(ylim_vals(1),1,nFlowers); repmat(ylim_vals(1),1,nFlowers) + diff_data(sort_idx)'], ...
%        'Color',rgb('Chocolate'),'LineStyle','-','LineWidth',4);


    
%Make it look nice
title(sprintf('Model evaluation - %s  r²=%.2f',analysis_string,rsquare_val),'Interpreter','None');
ylabel('Rating (z-score)');

set(gca,'YLim',ylim_vals); 
set(gca,'XLim',[0 nFlowers+1],'XTick',1:nFlowers,'XTickLabels',[]);
set(gca,'linewidth',1);

box off; %axis square;
set(gcf, 'Position', [1, 41, 1828, 917]); %set size on screen
flwr_width=0.035; %nFlowers/200
for i=1:nFlowers
    axes('pos',[0.13-flwr_width/2+i/(nFlowers+1)*0.775 0.11-flwr_width-flwr_width*mod(i+1,2) flwr_width flwr_width]) %[0.1300 0.1100 0.7750 0.8150]
    imshow(sprintf('%s%s.png',imgDir,labels{sort_idx(i)}));
end

saveas(f,[outputDir 'modelEval_' sortMethod '_' analysis_string '.png']);
close;

end

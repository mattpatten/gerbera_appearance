function cluster_pref_data

%load data
rootDir = get_dir([],'root');
load([rootDir 'labels_pref.mat']); %load preference data (reference RDMs)

%get data and parameters
vals = sortedFlowerData;
n = 1:50;

for i=n

    [clusters, ~, centroids, centDist] = perform_kmeans(vals, i, 0); %perform k-means clustering
    centDist = min(centDist,[],2); %the minimum distance will be the distance to the closest cluster (compared to the distance to all the other clusters)
    sumSquares(i) = sum(centDist);
end

figure;
f = subplot(1,2,1);
plot(sumSquares);
xlabel('Number of clusters');
ylabel('Sum of distance to nearest cluster');
title('Participant clusters');

subplot(1,2,2);
plot(abs(diff(sumSquares)));
xlabel('Number of clusters');
ylabel('Sum of distance reduction by adding i-th cluster');
title('Difference in clusters');

end
function totSum = kmeansWSS(nClusters, data)
for i = 1:nClusters
  [~, ~, sumd] = kmeans(data, i);
  totSum(i) = sum(sumd);
end
plot(totSum) % plot of totals versus number (same as index)
hold on
plot(totSum, '.')
xlabel('Number of clusters')
ylabel('Within-cluster sum-of-squares')
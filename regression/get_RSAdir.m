function RSAdir = get_RSAdir(outputDir, analysis_name, new_session, remove_influential_observers, flowers_to_remove, shuffling, bootstrapping, rankTransform, use_zscore, nClusters, cluster_num)

shuffling_text   = {'','shuffled_'};
bootstrap_text   = {'','btstrp_'};
rank_text        = {'','rank_'};
zscore_text      = {'','zscore_'};
influential_text = {'',['_not' sprintf('%d',flowers_to_remove)]};

RSAdir = [outputDir 'RSA' filesep analysis_name '_' new_session influential_text{1+remove_influential_observers} '_' ...
          shuffling_text{1+shuffling} bootstrap_text{1+bootstrapping} rank_text{1+rankTransform} zscore_text{1+use_zscore} ...
          num2str(nClusters) 'cl_' num2str(cluster_num) filesep];

end

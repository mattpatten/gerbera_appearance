function [ranks] = convert_to_ranks(data,order)

% input a vector or matrix and have each row turned into a rank
% order is 'ascend' or 'descend' as set by Matlab sort function
%
% Ties are sorted as one would expect - duplicating tie and leaving space thereafter
% data = [-2 2 0 -1 0 0 0 0 0 0];
% order = 'ascend';
%
% rank = [1 10 3 2 3 3 3 3 3 3];
%
% solution pinched off this matlab central query: 
% https://au.mathworks.com/matlabcentral/answers/33296-ranking-ordering-values-with-repeats
%
% Created by Matt Patten
% Created on 21/5/2020

%convert values to rank, and allow for ties (e.g., [0.5 4 5 5 7] converts to [1 2 3 3 5])
ranks = NaN(size(data));
for row=1:size(data,1)
    row_data = data(row,:);
    row_sorted = sort(row_data,order);
    [~, row_rank] = ismember(row_data,row_sorted);
    ranks(row,:) = row_rank;
    clear row_rank;
end

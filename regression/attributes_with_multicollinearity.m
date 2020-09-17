function attribute_idxs = attributes_with_multicollinearity(session,influential_flwrs_removed)

% A centralized function so we don't have to write this more than once, across scripts.
%
% Created by Matt Patten
% Created on 8th July, 2020

if ~influential_flwrs_removed
    % ALL FLOWERS
    %v2
    attribute_idxs = [              ...
        23:24                       ... %remove gradient threshold as flowers without gradient don't have a reasonable continuous value
        ,28,38,34,3,7,37,6,22,21,18 ... %VIF<10
        ,14,29,33,17,1,11           ... %VIF<5
        ];

%     %v1
%     attribute_idxs = [              ...
%         23:24                       ... %remove gradient threshold as flowers without gradient don't have a reasonable continuous value
%         ,28,38,34,3,7,37,6,22,21,18 ... %VIF<10
%         ,14,39,17,35,1,33           ... %VIF<5
%         ];

else
    % AFTER INFLUENTIAL OBSERVERS REMOVED
    if ismember(session,{'both_appeal','appeal'})
        %v2 - but with new influential observers
         attribute_idxs = [                 ...
             23:24                          ... %remove gradient threshold as flowers without gradient don't have a reasonable continuous value
            ,28,38,34,3,7,37,6,18,21,22,29  ... %VIF<10
            ,14,33,1,17,11                  ... %VIF<5
             ];
        
        %v2
%         attribute_idxs = [                 ...
%             23:24                          ... %remove gradient threshold as flowers without gradient don't have a reasonable continuous value
%             ,28,38,34,3,7,37,6,21,18,29,22 ... %VIF<10
%             ,33,14,1,17,11                 ...
%             ];

%        %v1
%        attribute_idxs = [                 ...
%            23:24                          ... %remove gradient threshold as flowers without gradient don't have a reasonable continuous value
%            ,28,38,34,3,7,37,6,21,18,22,39 ... %VIF<10
%            ,33,14,1,17,36,29              ... %VIF<5
%            ];
        
    elseif strcmpi(session,'both_interest')
         %v2 - but with new influential observers
         attribute_idxs = [              ...
             23:24                       ... %remove gradient threshold as flowers without gradient don't have a reasonable continuous value
             ,28,38,34,3,7,37,6,22,18,21 ... %VIF<10
             ,14,33,29,11,1,17,30           ... %VIF<5
             ];

%        %v1
%         attribute_idxs = [              ...
%             23:24                       ... %remove gradient threshold as flowers without gradient don't have a reasonable continuous value
%             ,28,38,34,3,7,37,6,22,18,21 ... %VIF<10
%             ,14,33,29,11,1,17              ... %VIF<5
%             ];

%    elseif strcmpi(session,'perceptual')
%        attribute_idxs = 1:35;
%        
%    elseif strcmpi(session,'objective')
%        attribute_idxs = 36:41;
        
    else
        error(['VIF analysis for attribute ' session ' has not been previously computed. Must be appeal or interest']);
    end

end          

end
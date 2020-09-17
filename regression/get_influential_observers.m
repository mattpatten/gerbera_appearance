function flowers_to_remove = get_influential_observers(session, remove_influential_observers)

% A centralized function so we don't have to write this more than once, across scripts.
%
% Created by Matt Patten
% Created on 10th July, 2020


rng('shuffle');

if remove_influential_observers
    if ismember(session,{'both_appeal','appeal'})
        %v2
        flowers_to_remove = [3 35 62 63];  %ascending by image number - appeal/both_appeal
        %v1
        %flowers_to_remove = [3 42 53 62 63];  %ascending by image number - appeal/both_appeal
    elseif strcmpi(session,'both_interest')
        %v2
        flowers_to_remove = [43 55 63];  %ascending by image number - appeal/both_appeal
        %v1
        %flowers_to_remove = [24 43 48 53 63];  %ascending by image number - both_interest
    end
else
    flowers_to_remove = [];
    %flowers_to_remove = 10+randi(50,1,3); %choose 3 numbers between 10 and 60 - this is their position/ranking, not the image number
end

end
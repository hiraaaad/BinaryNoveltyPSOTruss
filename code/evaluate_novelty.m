function novelty_score = evaluate_novelty(pos,archive,radius_novelty,pnovel)
pop = [pos;archive];

if size(archive,1) <= radius_novelty
    novelty_score =  repmat(1e15,size(pos,1),1);
else

if nargin == 3
    ps =  size(pos,1);
    novelty_score = zeros(1,ps);
    for i = 1 : ps
        x = pos(i,:);
        dist = sum(xor(x',pop'))';
        dist(:,2) = [1:1:size(pop,1)]';
        sorted_dist = sortrows(dist);
        delete_idx = find(sorted_dist(:,1)==0);
        sorted_dist(delete_idx,:) = []; % same as itself should be eliminated
        % lets see how many closest are in population
        % 
%         radius_novelty = sum(sorted_dist(:,1)==1);
%         radius_novelty = 10;
        nearest_neighbors = sorted_dist(1:radius_novelty,1);
        novelty_score(i) = mean(nearest_neighbors);
    end
    
elseif nargin == 4
    % reclacluate pnovelval wrt. pop
    dist = sum(xor(pnovel',pop'))';
    sorted_dist = sortrows(dist);
    sorted_dist(1,:) = []; % same as itself should be eliminated
    nearest_neighbors = sorted_dist(1:radius_novelty,1);
    novelty_score = mean(nearest_neighbors);
end
     
    
novelty_score = 1./novelty_score;
end

end

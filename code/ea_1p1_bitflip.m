function [x,kx] = ea_1p1_bitflip(x)
% function [x,kx] = ea_1p1_bitflip(x,infeasible_upper,feasible_upper)
    change = false;
%     while pass == false
        % bit flip
        for k = 1:length(x)
            if rand < 1/length(x)
                x(k) = 1-x(k);
                change = true;
            end
        end
        if change == false
            idx = randi([1,numel(x)]);
            x(idx) = 1-x(k);
        end
        kx = key_p(x);
%         if ~isKey(infeasible_upper,kx) && ~isKey(feasible_upper,kx)
%             pass = true; % just produce new unseen offspring 
%         end
end
% end
function y = sigmoid(x,phi)

if nargin < 2
    y = 1./(1 + exp(-x));
else
    y = 1./(1 + exp(-x/phi));
end
    

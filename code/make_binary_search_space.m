function binary_search_space = make_binary_search_space(benchmark)

switch benchmark
    case 10
        total = 2^10;
        initpop = dec2bin(0:2^10-1)' - '0';
        initpop = initpop';
        M = zeros(sqrt(total),sqrt(total));

end

%
hamming_distance_10bar_wrt_upperbound = b;
a=squareform(pdist(initpop,'minkowski',1));
b=a(end,:);
size(b)
b
hamming_distance_10bar_wrt_upperbound = b;

[val,rnk]=sort(hamming_distance_10bar_wrt_upperbound)
% val : hammind distance value
% rnk: real index in <
%

all_designs = [1:1:total];


M = ones(32,32)
X=[];
Y=[];
ctr=1
for i = 1: 32
    for j = 1:32 
        X(end+1)=i;
        Y(end+1)=j;
        M(i,j) = hamming_distance_10bar_wrt_upperbound(ctr);
        ctr=ctr+1;
    end
end
% M = 

k=0;
for i = 1: 32
    for j = 1:32
    k=k+1;
    v = make_rectangle(i,j);
    color = patch_color(k);
    f= 1:1:4;
    patch('Faces',f,'Vertices',v,'FaceColor',color)
    hold on
    end
end

C= zeros(54,3);
% FF=find(F>1e7);
f=f';
C(f,:) = repmat([255,255,0],numel(f),1);
scatter(X,Y,[],C,'filled')




row = floor(idx/32);
col = mod(idx,32)


    
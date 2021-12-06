function javab=find_largest_element_ind(A)
[b,ind]=max(A);
[c,ind2]=max(b);
javab=[ind(ind2) ind2];

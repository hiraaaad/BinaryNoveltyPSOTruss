function javab=calc_var(a,b,m,s) % calculates the variance of truncated normal distribution
%clear;clc;a=-rand(1,10);b=rand(1,10);m=rand(1,10),s=rand(1,10);

ta=(a-m)./s;
tb=(b-m)./s;
tapdf=normpdf(ta);
tbpdf=normpdf(tb);
denom=normcdf(tb)-normcdf(ta);
term1=(ta.*tapdf-tb.*tbpdf)./denom;
term2=(tapdf-tbpdf)./denom;
javab=s.^2.*(1+term1-term2.^2);

%r_d=normcdf(a,m,s);
%r_u=normcdf(b,m,s);
%for k=1:10000
%    r(k,:)=norminv(r_d+(r_u-r_d).*rand(1,10),m,s);
%end
%jjavab=[javab ;var(r)]


function [shortFvals,shortdist]=simplify_dist_metric(Fvals,disterr)

if size(Fvals,2)~=length(disterr),
    error('Fvals must have same number of columns as disterr')
end;

if disterr(1)~=0,
    error('expecting 0 dist as first element')
end;
Norig=size(Fvals,2);
Nnew=length(disterr)/2+1; %% new columns in Fvals
Nrows=size(Fvals,1);
shortFvals=zeros(Nrows,Nnew);
for i=1:Nrows,
    shortFvals(i,1)=Fvals(i,1); %% 0 error keeps same value
    for j=1:Nnew-1,
        pairvals=[Fvals(i,1+j) Fvals(i,Norig+1-j)];
        shortFvals(i,Nnew+1-j)=max(pairvals);
        
    end;
end;

for j=1:Nnew-1,
        pairvals=[disterr(1+j) disterr(Norig+1-j)];
        shortdist(Nnew+1-j)=max(pairvals);
    end;
function outputpara=qpSolver(fpara,spara,UPvalue,tervalue)
% min 0.5a'Qa-e'a
% s.t. 0<=a<=C
%t is between 0 and 2
[m,n]=size(fpara);
inipara=zeros(m,1);
trfpara=tril(fpara);
diagpara=diag(fpara);
interpara=inipara;
for j=1:n
    interpara(j,1)=inipara(j,1)-(spara/diagpara(j,1))*(fpara(j,:)*interpara(:,1)-1+trfpara(j,:)*(interpara(:,1)-inipara));
    if interpara(j,1)<0
        interpara(j,1)=0;
    elseif interpara(j,1)>UPvalue
        interpara(j,1)=UPvalue;
    else
        ;
    end
end
toutpara=[inipara,interpara];
iter=0;
while norm(toutpara(:,2)-toutpara(:,1))>tervalue && iter<100
    iter=iter+1;
    for j=1:n
        interpara(j,1)=toutpara(j,2)-(spara/diagpara(j,1))*(fpara(j,:)*interpara(:,1)-1+trfpara(j,:)*(interpara(:,1)-toutpara(:,2)));
        if interpara(j,1)<0
            interpara(j,1)=0;
        elseif interpara(j,1)>UPvalue
            interpara(j,1)=UPvalue;
        else
            ;
        end
    end
    toutpara(:,1)=[];
    toutpara=[toutpara,interpara];
end
outputpara=toutpara(:,2);
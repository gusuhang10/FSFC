function colwj=initialnizeW(dataj)
    [numsampA,~]=size(dataj);
    cj=mean(dataj);
    repdataj=repmat(cj,numsampA,1);
    Djp=(dataj-repdataj)'*(dataj-repdataj);
    [V,D]=eig(Djp);
    [~,pos]=min(diag(D));
    colwj=V(:,pos);
end
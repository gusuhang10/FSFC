function  [ttYpre,Wj]= FSFC(data,iniY,b1,b2,theta,lamda,col_W,U,m)
    % data:the given dataset with N*d, N is the number of samples and d is the
    % data dimension
    % iniY:the initial label set Y with N*1, where the element is between 1 and k
    % U and m are the initial fuzzy membership matrix and the fuzzy index, respectively
    %col_W: the column of W
    numsamp=size(data,1);
    numdimension=size(data,2);
    unitE=ones(numdimension,col_W);
    unitE(1:numdimension,1:numdimension)=eye(numdimension);
    termination=0;
    numiter=0;
    tYpre=iniY;
    while termination==0 && numiter<30 && (length(unique(tYpre))==length(unique(iniY)))
        numiter=numiter+1;
        ttYpre=tYpre;
        lenYpre=unique(ttYpre);
        numY=length(lenYpre);
        dis=zeros(numsamp,numY);
        Wj=zeros(numdimension,col_W,numY);
        centerj=zeros(numY,numdimension);
        for i=1:numY
            dataj=data(ttYpre==lenYpre(i),:);
            odataj=data(ttYpre~=lenYpre(i),:);
            centerj(i,:)=mean(dataj);
            for p=1:col_W
                if p==1
                    colWj=zeros(numdimension,1);
                elseif norm(colWj)~=0
                    colWj=colWj/norm(colWj);
                end
                dataj=dataj-dataj*(colWj*colWj');
                odataj=odataj-odataj*(colWj*colWj');
                iniW=initialnizeW(dataj);
                colWj=funSolver(dataj,odataj,iniW,b1,b2,theta,lamda,unitE(:,p));
                Wj(:,p,i)=colWj;
            end
        end
        for i=1:numY
            eachS(:,:,:,i)=data*Wj(:,:,i)-repmat(centerj(i,:)*Wj(:,:,i),numsamp,1);
            com_eachS(:,i)=sum(abs(eachS(:,:,:,i)),2);
        end
        U=ones(size(U,1),size(U,2))./((com_eachS./repmat(sum(com_eachS,2),1,numY)).^m);
        U=U./repmat(sum(U,2),1,numY);
        for i=1:numY
            eS=eachS(:,:,:,i);
            for j=1:numsamp
                dis(j,i)=U(j,i)^m*norm(eS(j,:));
            end
        end
        [~,tYpre]=min(dis,[],2);
        if accuracy(ttYpre, tYpre)>=1
            termination=1;
        end
    end
end
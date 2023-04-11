function wjp=funSolver(Cdata,Otdata,w0,g1,g2,sigma,lamda,unitEcol)
    data=[Cdata;Otdata];%the given dataset
    numCdata=size(Cdata,1);
    dim=size(Cdata,2);
    numOtdata=size(Otdata,1);
    e1=ones(numCdata,1);
    e2=ones(numOtdata,1);
    %eunitE=ones(size(w0,1),1);
    diff=10^-3;
    wdiff=1;    
    iterH=30;
    h=0;
    while wdiff>diff && h<iterH
        h=h+1;
        A1=sign(w0'*(data'*data)*w0-1);
        A2=diag(sign((Otdata-1/numCdata*e2*e1'*Cdata-lamda*e2*unitEcol')*w0));
        A3=(Cdata-repmat(mean(Cdata),numCdata,1))'*(Cdata-repmat(mean(Cdata),numCdata,1));
        A4=((1+lamda)*eye(dim)+g1*A3+sigma*A1*(data'*data))\(Otdata'-1/numCdata*Cdata'*e1*e2'-lamda*unitEcol*e2')*A2;
        H=A2*(Otdata-1/numCdata*e2*e1'*Cdata-lamda*e2*unitEcol')*A4;
        H=(H+H')/2;
        gamma=qpSOR(H,0.7,g2,0.05);
        wjp=A4*gamma;
        wdiff=norm(wjp-w0);
        w0=wjp;
    end
end
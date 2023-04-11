function bestalpha=qpSOR(Q,t,C,smallvalue)
%This refers to Prof. Shao's work at http://www.optimal-group.org/Resource.html
% min 0.5a'Qa-e'a
% s.t. 0<=a<=C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Q���������
%    t�ǣ�0��2���Ĳ�����
%    C���Ͻ磻
%    smallvalue����ֹ������
%    bestalpha�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ʼ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(Q);
alpha0=zeros(m,1);
% i=0;
L=tril(Q);
E=diag(Q);
twinalpha=alpha0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�ж�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:n
%     i=i+1;
    twinalpha(j,1)=alpha0(j,1)-(t/E(j,1))*(Q(j,:)*twinalpha(:,1)-1+L(j,:)*(twinalpha(:,1)-alpha0));
    if twinalpha(j,1)<0
        twinalpha(j,1)=0;
    elseif twinalpha(j,1)>C
        twinalpha(j,1)=C;
    else
        ;
    end
%     twinalpha(j,1)=tmp_alpha;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ѭ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=[alpha0,twinalpha];
ite=0;
while norm(alpha(:,2)-alpha(:,1))>smallvalue && ite<100
    ite=ite+1;
    for j=1:n
        twinalpha(j,1)=alpha(j,2)-(t/E(j,1))*(Q(j,:)*twinalpha(:,1)-1+L(j,:)*(twinalpha(:,1)-alpha(:,2)));
        if twinalpha(j,1)<0
            twinalpha(j,1)=0;
        elseif twinalpha(j,1)>C
            twinalpha(j,1)=C;
        else
            ;
        end
%         twinalpha(j,1)=tmp_alpha;
    end
    alpha(:,1)=[];
    alpha=[alpha,twinalpha];
%     twinalpha=alpha(:,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bestalpha=alpha(:,2);



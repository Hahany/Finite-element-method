function [connect_mat,A] = connect_mat( numnodx,numnody,nel)
%�����������Ľڵ���Ŀ���͵�Ԫ���ɶ�
%������Ӿ���ÿ����Ԫ�漰�Ľڵ�ı��
xn = 1:(numnodx*numnody);%����һ�����
A = reshape(xn,numnody,numnodx);%ͬ��״���
B = zeros(numnody,numnodx);
for i = 1:numnody
B(i,:)=A(numnody+1-i,:);%��������һֱ�����Һ�����Ϊ����
end
A=B;
for i = 1:(numnodx-1)*(numnody-1)
    yg = rem(i,numnody-1);%xg��ʾ��ԪΪ�±߽�����ڼ���
    if yg == 0
        yg = numnody-1;
    end
    xg = ceil(i/(numnody-1));%��߽������ڼ���
    a = A(numnody-yg:numnody-yg+1,xg:xg+1);%���С������ֱ�˾������Ӿ���
    a_veg = a(:);
    connect_mat(i,1:nel) =a_veg([2 4 3 1]);   %��ʱ����룬���׼�Ľڵ㵥Ԫ����һ�£�����Ϊ1
end
end
function [ estress,estrain ] = stress_strain(e, tota,connect,D,B)
%�ú�����������Ӧ��ֵ��Ӧ��ֵ
% ��������λ�õ�Ӧ����Ӧ�䣬xi��eta�������⣬��B��ֵ������
nodes = connect(e,:);
unodes = zeros(2*length(nodes),1);
for i = 1:length(nodes)
    unodes(2*i-1) = 2*nodes(i)-1;
    unodes(2*i) = 2*nodes(i);
end
ae = tota(unodes);
estress = D*B*ae;
estrain = B*ae;
end


function [ K ] = totalstiff( D,connect,coord,numel,nel,numnod )
%���㵥Ԫ�նȾ���Ke���������ܸն�K
% connect �ǵ�Ԫ�ڵ��� coord�ǵ�Ԫ�ڵ��Ӧ������   D �Ǳ�������ϵ������ numel�ǵ�Ԫ������ nel��Ԫ�ڵ���
% numnod���ܽ����
  gaussxi=1/sqrt(3)*[-1,1,1,-1];  %��˹��xi����
  gausseta=1/sqrt(3)*[-1,-1,1,1];  %��˹��eta����
  K = sparse(2*numnod,2*numnod); % �նȾ���[K]����ʼ��Ϊ0��ʹ��ϡ�����洢
for e = 1:numel %ͬһά���������Ȼ����Ԫ��ɨ��
    Ke = zeros(8,8);
    for i_gauss = 1:4 %���ø�˹���ַ���
        [Bgauss,J] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%�������Ԫ��˹�㴦��Ӧ�����
        Ke = Ke+ Bgauss'*D*Bgauss*det(J); %��Ԫ�նȾ���
    end
  totn = connect(e,:); %��Ԫ�ڵ���
  totn2 = zeros(2*length(totn),1);  %��ˮƽ����ֱλ�ƣ���������һ��
  for i =1:nel
      totn2(2*i-1:2*i,1) =  [2*totn(i)-1;2*totn(i)];  %�����Ԫ�ڵ��Ӧ���к�
  end
  K(totn2,totn2) = K(totn2,totn2)+Ke;  %��װ����ն�
end

end


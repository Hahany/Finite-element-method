function [B,J] = B(e,coord,connect,gaussxi,gausseta)
%���������ڼ���gauss�㴦��Ӧ�����B��ֵ��
nodes = connect(e,:);%����κ������ڵ㣩���
xe    = coord(nodes,:); %��������ϵ�£���Ԫ�ڵ������
% �����ֲ�����ϵ�µ��κ���N{i}
  J_par=zeros(2,4);
    %�����Ÿ��Ⱦ��� J  %%����һ�в���Ԫ���������м���ķ�����N{i} =@(xi,eta) 1/4*(1+Ixi(i)*xi)*(1+Ieta(i)*eta); %�κ���   ע�⣺���õ�ʱ��Ҫ��{}
  B=zeros(3,8); %�����ʼ��Ӧ�����B
  Ixi=[-1,1,1,-1];
  Ieta=[-1,-1,1,1]; %���ſ���
              for k = 1:4
                  eta =gausseta ;%gauss���Ӧ�ľֲ����� eta
                  xi = gaussxi;%gauss���Ӧ�ľֲ����� xi
                  J_par(1,k) = Ixi(k)*1/4*(1+Ieta(k)*eta);  %�Ÿ��Ⱦ���ƫ��ϵ����һ�У���i����˹���ֵ
                  J_par(2,k) = Ieta(k)*1/4*(1+Ixi(k)*xi);  %�Ÿ������ƫ��ϵ���ڶ��У���i����˹���ֵ
              end
%               Jgauss = permute(J_par(i,:,:),[2 3 1]);  %��i����˹�㴦���Ÿ��Ⱦ���
                  J =  J_par*xe;  
              for i = 1:4%��i����˹����Ÿ��Ⱦ���J
                  N_ipar(:,i) = inv(J)*J_par(:,i);
                  B(1,2*i-1)=N_ipar(1,i);
                  B(2,2*i)=N_ipar(2,i);
                  B(3,2*i-1:2*i)=[N_ipar(2,i),N_ipar(1,i)];
              end
return
end

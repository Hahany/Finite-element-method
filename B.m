function [B,J] = B(e,coord,connect,gaussxi,gausseta)
%本函数用于计算gauss点处的应变矩阵B的值。
nodes = connect(e,:);%相关形函数（节点）编号
xe    = coord(nodes,:); %整体坐标系下，单元节点的坐标
% 构建局部坐标系下的形函数N{i}
  J_par=zeros(2,4);
    %计算雅各比矩阵 J  %%另外一中采用元胞方法进行计算的方法：N{i} =@(xi,eta) 1/4*(1+Ixi(i)*xi)*(1+Ieta(i)*eta); %形函数   注意：调用的时候要用{}
  B=zeros(3,8); %定义初始的应变矩阵B
  Ixi=[-1,1,1,-1];
  Ieta=[-1,-1,1,1]; %符号控制
              for k = 1:4
                  eta =gausseta ;%gauss点对应的局部坐标 eta
                  xi = gaussxi;%gauss点对应的局部坐标 xi
                  J_par(1,k) = Ixi(k)*1/4*(1+Ieta(k)*eta);  %雅各比矩阵偏导系数第一行，第i个高斯点的值
                  J_par(2,k) = Ieta(k)*1/4*(1+Ixi(k)*xi);  %雅各矩阵比偏导系数第二行，第i个高斯点的值
              end
%               Jgauss = permute(J_par(i,:,:),[2 3 1]);  %第i个高斯点处的雅各比矩阵
                  J =  J_par*xe;  
              for i = 1:4%第i个高斯点的雅各比矩阵J
                  N_ipar(:,i) = inv(J)*J_par(:,i);
                  B(1,2*i-1)=N_ipar(1,i);
                  B(2,2*i)=N_ipar(2,i);
                  B(3,2*i-1:2*i)=[N_ipar(2,i),N_ipar(1,i)];
              end
return
end

function [ K ] = totalstiff( D,connect,coord,numel,nel,numnod )
%计算单元刚度矩阵Ke，并集成总刚度K
% connect 是单元节点编号 coord是单元节点对应的坐标   D 是本构方程系数矩阵 numel是单元总数量 nel单元节点数
% numnod是总结点数
  gaussxi=1/sqrt(3)*[-1,1,1,-1];  %高斯点xi坐标
  gausseta=1/sqrt(3)*[-1,-1,1,1];  %高斯点eta坐标
  K = sparse(2*numnod,2*numnod); % 刚度矩阵[K]，初始化为0，使用稀疏矩阵存储
for e = 1:numel %同一维的情况，依然按单元来扫描
    Ke = zeros(8,8);
    for i_gauss = 1:4 %采用高斯积分方案
        [Bgauss,J] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%计算各单元高斯点处的应变矩阵
        Ke = Ke+ Bgauss'*D*Bgauss*det(J); %单元刚度矩阵
    end
  totn = connect(e,:); %单元节点编号
  totn2 = zeros(2*length(totn),1);  %有水平和竖直位移，矩阵扩大一倍
  for i =1:nel
      totn2(2*i-1:2*i,1) =  [2*totn(i)-1;2*totn(i)];  %扩大后单元节点对应的行号
  end
  K(totn2,totn2) = K(totn2,totn2)+Ke;  %组装整体刚度
end

end


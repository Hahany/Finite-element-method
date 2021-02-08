function [ sigma0_max, sigma0,tota0,smnode0 ] = solver( K,P0,D,numnod,numel,nel,bnoderow,outp,Ixi,Ieta,coord,connect,gaussxi,gausseta )
%该函数用来求解给定总刚度，荷载和本构关系时的各单元最大等效应力
%   除了以上的变量以外，其余均为网格划分时的常数。
%该方程的缺点：不是变刚度的方程。
a0=K\P0;%计算P0作用时弹性解
tota0 = zeros(2*numnod,1);
rowk0=setdiff(1:2*numnod,bnoderow);
tota0(rowk0)=a0;
u0=zeros(numnod,1); %P0荷载对应的水平位移
v0=zeros(numnod,1); %P0荷载对应的纵向位移
for i = 1:numnod
    u0(i) = tota0(2*i-1);%P0荷载对应取出水平位移
    v0(i) = tota0(2*i); %P0荷载对应取出纵向位移
end
sigma0=zeros(numel,nel*3);
strain0=zeros(numel,nel*3);
for e = 1:numel
    %计算gauss点处的应力和应变
    for i_gauss = 1:4 %逐个取高斯点
        [Bgauss,~] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%计算各单元高斯点处的应变矩阵
        [sigma0(e,3*i_gauss-2:3*i_gauss),~]=stress_strain(e,tota0,connect,D,Bgauss);
    end
end
smnode0=sigma0*outp; %得到节点处磨平后的应力值
sigma0_max = zeros(numel,1);
for e =1:numel
    sigma0_i = zeros(nel,1);
    for i =1:nel
        [~,sigma0_i(i)] = sigmainner(smnode0, e,Ixi(i),Ieta(i));
    end
    sigma0_max(e) =max(sigma0_i);  %取得单元的最大等效应力，行数对应单元数
end

end


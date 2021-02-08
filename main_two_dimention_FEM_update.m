%% 这是一个平面应力弹塑性有限元程序
%% 四边形网格
clc;  % 清空命令行窗口
clear; %清除工作空间
close all; %关闭所有图像
%% 参数设置及网格划分
error_L2s = [];
E=1e5; %弹性模量
Et=0.2*E;
nu = 0.25; %泊松比
sigmas=25; %mises屈服应力强度
Lx = 8; %定义单元右边界（左边界为0）
Ly = 1;%定义单元上边界;
numelx = 160;%定义分割的x方向单元数目
numely = 20;%定义分割的y方向单元数目
hx = Lx/numelx;%x方向上的单元长度
hy = Ly/numely;%y方向上的单元长度
numel = numelx*numely;%单元的数目
u_b = zeros(2*(numely+1),1); %定义第一类边界条件，x=0，处位移为0
numnodx = numelx + 1; % x方向节点个数比单元个数多1
numnody = numely + 1; % y方向节点个数比单元个数多1
numnod = numnodx*numnody;%总的节点个数
nel = 4;%每个单元的节点数目,即每个单元上有几个形函数参与作用
coordx = linspace(0,Lx,numnodx)'; %等分节点的坐标（为了方便，我这里采用等分的方式，事实上单元长度可以不一致，非均匀网格）
coordy = linspace(0,Ly,numnody)'; %等分节点的坐标（为了方便，我这里采用等分的方式，事实上单元长度可以不一致）
[X, Y] = meshgrid(coordx,coordy);%张成网格，X和Y分别表示对应位置的横纵坐标
coord = [X(:) Y(:)];%把网格一列一列扯开，coord的每一行是对应节点的坐标，按顺序排列
[connect,A] = connect_mat(numnodx,numnody,nel);%连接矩阵，表示每个单元周围的节点编号，也就是涉及的形函数编号
nodeb = 1:numnody; % 强制性边界点的编号
nodeval = u_b; %假设边界值都为u_b
P=zeros(2*numnod,1);
P(end) = -1;  %节点荷载

%%  计算系数矩阵K和右端荷载矩阵，单刚组装总刚，以及边界处理
D=elasm(E,nu); %广义胡克定律系数矩阵
K= totalstiff(D,connect,coord,numel,nel,numnod );
bnoderow = zeros(2*length(nodeb),1);
for i =1:length(nodeb)
    bnoderow(2*i-1:2*i,1)=[2*nodeb(i)-1;2*nodeb(i)];
end
K(bnoderow,:)=[];
K(:,bnoderow)=[]; %删除掉位移为0的点对应，刚度矩阵的行和列
P(bnoderow)=[];  %删除固定边界对应的荷载分量
a = K\P;  %计算弹性解
tota = zeros(2*numnod,1);
rowke = setdiff(1:numnod,nodeb);
rowk=setdiff(1:2*numnod,bnoderow);
tota(rowk)=a;
u=zeros(numnod,1); %水平位移
v=zeros(numnod,1); %纵向位移
wholea=zeros(numnod,1); %总位移
for i = 1:numnod
    u(i) = tota(2*i-1)*1000;%取出水平位移
    v(i) = tota(2*i)*1000; %取出纵向位移
    wholea(i) = sqrt(u(i)^2+v(i)^2)*1000; %计算总位移，注意和tota的区别
end

%% 绘图
u = reshape(u,numnody,numnodx);   %重构横向位移矩阵
v = reshape(v,numnody,numnodx); %重构纵向位移矩阵
subplot(2,2,1)
surface(coordx,coordy,v);colorbar;  %绘制位移表面
shading interp
ylim([0,1])    %控制轴坐标范围
xlim([0,8])
axis equal   %轴单位相同
subplot(2,2,3)
mesh(coordx,coordy,v)
ylim([0,1])    %控制轴坐标范围
xlim([0,8])
% axis equal   %轴单位相同

%% 结果后处理
%应力计算
%应变计算
gaussxi=1/sqrt(3)*[-1,1,1,-1];  %高斯点xi坐标
gausseta=1/sqrt(3)*[-1,-1,1,1];
stress=zeros(numel,nel*3);
strain=zeros(numel,nel*3);
for e = 1:numel
    %计算gauss点处的应力和应变
    for i_gauss = 1:4 %逐个取高斯点
        [Bgauss,~] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%计算各单元高斯点处的应变矩阵
        %1-3列为第一个高斯点的值，依此类推
        [stress(e,3*i_gauss-2:3*i_gauss),strain(e,3*i_gauss-2:3*i_gauss)]=stress_strain(e,tota,connect,D,Bgauss);
    end
end
           
%% 磨平
%利用高斯点外推节点应力 sigma每一行是每一个单元，高斯点处的应力值，分别为sigmax，sigmay,tauxy
%注意到sigma行向量的特点，构建外推矩阵
apara=1+sqrt(3)/2;
bpara=-1/2;
cpara=1-sqrt(3)/2;
aeye = eye(3)*apara;
beye = eye(3)*bpara;
ceye = eye(3)*cpara;
eoutp=eye(4)*apara+diag([1,1,1]*bpara,1)+diag([1,1,1]*bpara,-1)+diag([1,1]*cpara,2)+diag([1,1]*cpara,-2)+diag(bpara,3)+diag(bpara,-3);
outp=zeros(12,12);
for i =1:4  %利用for循环构造外推矩阵
    outp(3*i-2:3*i,3*i-2:3*i)=aeye;
    if i < 4
    outp(3*i-2:3*i,3*i+1:3*i+3)=beye;
    outp(3*i+1:3*i+3,3*i-2:3*i)=beye;
    if i < 3
    outp(3*i-2:3*i,3*i+4:3*i+6)=ceye;
    outp(3*i+4:3*i+6,3*i-2:3*i)=ceye;
    end
    else
        outp(1:3,10:12)=beye;
        outp(10:12,1:3)=beye;
    end
end
smnode=stress*outp; %得到节点处磨平后的应力值

%% 计算每一个单元的最大等效应力
%单元的最大等效应力在单元节点处，见test.m中的绘图结果。
sigma_max = zeros(32,1);
Ixi=[-1,1,1,-1];
Ieta=[-1,-1,1,1];
sigplotf=zeros(nel,numel);
for e =1:numel
    sigma_i = zeros(4,1);
    for i =1:4
        [~,sigma_i(i)] = sigmainner(smnode, e,Ixi(i),Ieta(i));
    end
    sigma_max(e) =max(sigma_i);  %取得单元的最大等效应力，行数对应单元数
    sigplotf(:,e)=sigma_i;
end

%% 变刚度切线法求弹塑性问题
%首先计算出最大弹性荷载P0对应的弹性解
totsigma = max(sigma_max); %判断sigma_max 与sigmas之间的关系
L = totsigma/sigmas;
P0=1/L*P;
n = 10; %设定剩余荷载的增量数量
deltP=1/n*(1-1/L)*P;
[sigma0_max,sigma0,tota0,smnode0] = solver( K,P0,D,numnod,numel,nel,bnoderow,outp ,Ixi,Ieta,coord,connect,gaussxi,gausseta);

%% 施加固定荷载增量deltP
%这里的关键在于变单元刚度
% sigma1_max = sigma0_max; %初始化最大等效应力为最大弹性荷载时的的值
elsigma_max = sigma0_max;
gausstress=sigma0;
Kepll=sparse(2*numnod,2*numnod);
me=zeros(numel,1);
for e = 1:numel
    if (elsigma_max(e)-sigmas)*(abs(elsigma_max(e)-sigmas)>1e-10)<0   %设定精度为1e-10.
        me(e)=1;
    else
        me(e)=0;
    end
    Ke = zeros(8,8);
    for i_gauss = 1:4 %采用高斯积分方案
        [Bgauss,J] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%计算各单元高斯点处的应变矩阵
        Dep = e_plasm( sigmas,E,Et,nu,gausstress(e,(3*i_gauss-2:3*i_gauss)));
        Dhat = me(e)*D+(1-me(e))*Dep;
        Ke = Ke+ Bgauss'*Dhat*Bgauss*det(J); %单元刚度矩阵
    end
    totn = connect(e,:); %单元节点编号
    totn2 = zeros(2*length(totn),1);  %有水平和竖直位移，矩阵扩大一倍
    for i =1:nel
        totn2(2*i-1:2*i,1) =  [2*totn(i)-1;2*totn(i)];  %扩大后单元节点对应的行号
    end
    Kepll(totn2,totn2) = Kepll(totn2,totn2)+Ke;  %组装整体刚度
    
end
Kepll(bnoderow,:) = [];
Kepll(:,bnoderow) = [];
[kx,ky]=size(Kepll);
Kep=zeros(kx,ky,n+1);
Kep(:,:,1) = full(K);
delu=zeros(n,length(P));
deltotu = zeros(2*numnod,n);
delsigma = zeros(numel,3*nel,n);
delsigma_max =zeros(e,n);
delsigma_max(:,1)=sigma0_max;
m=zeros(numel,n);
m(:,1)=me;
sigmatot(:,:,1)=sigma0;
smdelsigma(:,:,1)=smnode0;
totu=tota0;
j = 0;
sigplot=zeros(nel,numel);
for ip = 1:n;
    jwhile=1;
    kwhile=0;
%     while 1
    kwhile=kwhile+1;
    delu(ip,:)=Kep(:,:,ip)\deltP;  %deltaU_2
    deltotu(:,ip) = zeros(2*numnod,1);
    deltotu(rowk,ip)=delu(ip,:);
    for e = 1:numel
    %计算gauss点处的应力和应变
        for i_gauss = 1:4 %逐个取高斯点
            [Bgauss,~] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%计算各单元高斯点处的应变矩阵
            Dep = e_plasm( sigmas,E,Et,nu,sigmatot(e,(3*i_gauss-2:3*i_gauss),ip));
            Dhat = m(e,ip)*D+(1-m(e,ip))*Dep;
            [delsigma(e,3*i_gauss-2:3*i_gauss,ip),~]=stress_strain(e,deltotu(:,ip),connect,Dhat,Bgauss);
        end
        s = zeros(ip,3*nel);
        sm = zeros(ip,3*nel);
        for is=1:ip
            s(is,:) =delsigma(e,:,is);
            sm(is,:)=delsigma(e,:,is)*outp;
        end
        sigmatot(e,:,ip+1)=sigmatot(e,:,1)+sum(s);
        smdelsigma(e,:,ip+1)=smdelsigma(e,:,1)+sum(sm); %每个单元应力磨平后的节点值
        delsigma_i = zeros(nel,1);
        for i =1:nel
            [~,delsigma_i(i)] = sigmainner(smdelsigma(:,:,ip+1), e,Ixi(i),Ieta(i));
        end
        sigplot(:,e)=delsigma_i;
        delsigma_max(e,ip+1) =max(delsigma_i);  %取得单元的最大等效应力，行数对应单元数
        if (delsigma_max(e,ip)-sigmas)*(abs(delsigma_max(e)-sigmas)>1e-10)<0 
            if (delsigma_max(e,ip+1)-sigmas)*(abs(delsigma_max(e,ip+1)-sigmas)>1e-10)<0
                m(e,ip+1) = 1;
            else
                m(e,ip+1)=(sigmas-delsigma_max(e,ip))/(delsigma_max(e,ip+1)-delsigma_max(e,ip));
                jwhile=jwhile+1;
            end
        else
            m(e,ip+1)=0;
        end
        Ke = zeros(8,8);
        for i_gauss = 1:4 %采用高斯积分方案
            [Bgauss,J] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%计算各单元高斯点处的应变矩阵
            Dep = e_plasm( sigmas,E,Et,nu,sigmatot(e,(3*i_gauss-2:3*i_gauss),ip+1));
            Dhat = m(e,ip+1)*D+(1-m(e,ip+1))*Dep;
            Ke = Ke+ Bgauss'*Dhat*Bgauss*det(J); %单元刚度矩阵
        end
        totnin = connect(e,:); %单元节点编号
        [~,In] = setdiff(totnin,nodeb);
        totn = totnin(sort(In));
        [~,kkill,~] = intersect(totnin,nodeb);
        if isempty(kkill)==0
                kl = [2*kkill-1;2*kkill];
                Ke(kl,:) = [];
                Ke(:,kl) = [];
        end        
        totn2 = zeros(2*length(totn),1);  %有水平和竖直位移，矩阵扩大一倍
        for i =1:length(totn)
            nr = find(rowke==totn(i));
          totn2(2*i-1:2*i,1) =  [2*nr-1;2*nr];  %扩大后单元节点对应的行号
        end
        Kep(totn2,totn2,ip+1) = Kep(totn2,totn2,ip+1)+Ke;  %组装整体刚度
    end
    epm(:,kwhile)=m(:,ip+1);
%     if jwhile >1
%         Kep(:,:,ip) = Kep(:,:,ip+1);
%         continue
%     else
%         break
% %         if norm(epm(:,jwhile+1)-epm(:,ip))<1e-3
% %             j=j+1
% %             break
% %         else
% %             Kep(:,:,ip) = Kep(:,:,ip+1);
% %         end
%     end
%     end
    totu = totu + deltotu(:,ip);
end
uep=zeros(numnod,1); %水平位移
vep=zeros(numnod,1); %纵向位移
wholea=zeros(numnod,1); %总位移
for i = 1:numnod
    uep(i) = totu(2*i-1)*1000;%取出水平位移
    vep(i) = totu(2*i)*1000; %取出纵向位移
end

%% 绘图
uep = reshape(uep,numnody,numnodx);   %重构横向位移矩阵
vep = reshape(vep,numnody,numnodx); %重构纵向位移矩阵
subplot(2,2,2)
surface(coordx,coordy,vep);colorbar;  %绘制位移表面
shading interp
ylim([0,1])    %控制轴坐标范围
xlim([0,8])
axis equal   %轴单位相同
subplot(2,2,4)
mesh(coordx,coordy,vep)
ylim([0,1])    %控制轴坐标范围
xlim([0,8])
% axis equal   %轴单位相同



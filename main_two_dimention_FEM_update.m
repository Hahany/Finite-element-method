%% ����һ��ƽ��Ӧ������������Ԫ����
%% �ı�������
clc;  % ��������д���
clear; %��������ռ�
close all; %�ر�����ͼ��
%% �������ü����񻮷�
error_L2s = [];
E=1e5; %����ģ��
Et=0.2*E;
nu = 0.25; %���ɱ�
sigmas=25; %mises����Ӧ��ǿ��
Lx = 8; %���嵥Ԫ�ұ߽磨��߽�Ϊ0��
Ly = 1;%���嵥Ԫ�ϱ߽�;
numelx = 160;%����ָ��x����Ԫ��Ŀ
numely = 20;%����ָ��y����Ԫ��Ŀ
hx = Lx/numelx;%x�����ϵĵ�Ԫ����
hy = Ly/numely;%y�����ϵĵ�Ԫ����
numel = numelx*numely;%��Ԫ����Ŀ
u_b = zeros(2*(numely+1),1); %�����һ��߽�������x=0����λ��Ϊ0
numnodx = numelx + 1; % x����ڵ�����ȵ�Ԫ������1
numnody = numely + 1; % y����ڵ�����ȵ�Ԫ������1
numnod = numnodx*numnody;%�ܵĽڵ����
nel = 4;%ÿ����Ԫ�Ľڵ���Ŀ,��ÿ����Ԫ���м����κ�����������
coordx = linspace(0,Lx,numnodx)'; %�ȷֽڵ�����꣨Ϊ�˷��㣬��������õȷֵķ�ʽ����ʵ�ϵ�Ԫ���ȿ��Բ�һ�£��Ǿ�������
coordy = linspace(0,Ly,numnody)'; %�ȷֽڵ�����꣨Ϊ�˷��㣬��������õȷֵķ�ʽ����ʵ�ϵ�Ԫ���ȿ��Բ�һ�£�
[X, Y] = meshgrid(coordx,coordy);%�ų�����X��Y�ֱ��ʾ��Ӧλ�õĺ�������
coord = [X(:) Y(:)];%������һ��һ�г�����coord��ÿһ���Ƕ�Ӧ�ڵ�����꣬��˳������
[connect,A] = connect_mat(numnodx,numnody,nel);%���Ӿ��󣬱�ʾÿ����Ԫ��Χ�Ľڵ��ţ�Ҳ�����漰���κ������
nodeb = 1:numnody; % ǿ���Ա߽��ı��
nodeval = u_b; %����߽�ֵ��Ϊu_b
P=zeros(2*numnod,1);
P(end) = -1;  %�ڵ����

%%  ����ϵ������K���Ҷ˺��ؾ��󣬵�����װ�ܸգ��Լ��߽紦��
D=elasm(E,nu); %������˶���ϵ������
K= totalstiff(D,connect,coord,numel,nel,numnod );
bnoderow = zeros(2*length(nodeb),1);
for i =1:length(nodeb)
    bnoderow(2*i-1:2*i,1)=[2*nodeb(i)-1;2*nodeb(i)];
end
K(bnoderow,:)=[];
K(:,bnoderow)=[]; %ɾ����λ��Ϊ0�ĵ��Ӧ���նȾ�����к���
P(bnoderow)=[];  %ɾ���̶��߽��Ӧ�ĺ��ط���
a = K\P;  %���㵯�Խ�
tota = zeros(2*numnod,1);
rowke = setdiff(1:numnod,nodeb);
rowk=setdiff(1:2*numnod,bnoderow);
tota(rowk)=a;
u=zeros(numnod,1); %ˮƽλ��
v=zeros(numnod,1); %����λ��
wholea=zeros(numnod,1); %��λ��
for i = 1:numnod
    u(i) = tota(2*i-1)*1000;%ȡ��ˮƽλ��
    v(i) = tota(2*i)*1000; %ȡ������λ��
    wholea(i) = sqrt(u(i)^2+v(i)^2)*1000; %������λ�ƣ�ע���tota������
end

%% ��ͼ
u = reshape(u,numnody,numnodx);   %�ع�����λ�ƾ���
v = reshape(v,numnody,numnodx); %�ع�����λ�ƾ���
subplot(2,2,1)
surface(coordx,coordy,v);colorbar;  %����λ�Ʊ���
shading interp
ylim([0,1])    %���������귶Χ
xlim([0,8])
axis equal   %�ᵥλ��ͬ
subplot(2,2,3)
mesh(coordx,coordy,v)
ylim([0,1])    %���������귶Χ
xlim([0,8])
% axis equal   %�ᵥλ��ͬ

%% �������
%Ӧ������
%Ӧ�����
gaussxi=1/sqrt(3)*[-1,1,1,-1];  %��˹��xi����
gausseta=1/sqrt(3)*[-1,-1,1,1];
stress=zeros(numel,nel*3);
strain=zeros(numel,nel*3);
for e = 1:numel
    %����gauss�㴦��Ӧ����Ӧ��
    for i_gauss = 1:4 %���ȡ��˹��
        [Bgauss,~] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%�������Ԫ��˹�㴦��Ӧ�����
        %1-3��Ϊ��һ����˹���ֵ����������
        [stress(e,3*i_gauss-2:3*i_gauss),strain(e,3*i_gauss-2:3*i_gauss)]=stress_strain(e,tota,connect,D,Bgauss);
    end
end
           
%% ĥƽ
%���ø�˹�����ƽڵ�Ӧ�� sigmaÿһ����ÿһ����Ԫ����˹�㴦��Ӧ��ֵ���ֱ�Ϊsigmax��sigmay,tauxy
%ע�⵽sigma���������ص㣬�������ƾ���
apara=1+sqrt(3)/2;
bpara=-1/2;
cpara=1-sqrt(3)/2;
aeye = eye(3)*apara;
beye = eye(3)*bpara;
ceye = eye(3)*cpara;
eoutp=eye(4)*apara+diag([1,1,1]*bpara,1)+diag([1,1,1]*bpara,-1)+diag([1,1]*cpara,2)+diag([1,1]*cpara,-2)+diag(bpara,3)+diag(bpara,-3);
outp=zeros(12,12);
for i =1:4  %����forѭ���������ƾ���
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
smnode=stress*outp; %�õ��ڵ㴦ĥƽ���Ӧ��ֵ

%% ����ÿһ����Ԫ������ЧӦ��
%��Ԫ������ЧӦ���ڵ�Ԫ�ڵ㴦����test.m�еĻ�ͼ�����
sigma_max = zeros(32,1);
Ixi=[-1,1,1,-1];
Ieta=[-1,-1,1,1];
sigplotf=zeros(nel,numel);
for e =1:numel
    sigma_i = zeros(4,1);
    for i =1:4
        [~,sigma_i(i)] = sigmainner(smnode, e,Ixi(i),Ieta(i));
    end
    sigma_max(e) =max(sigma_i);  %ȡ�õ�Ԫ������ЧӦ����������Ӧ��Ԫ��
    sigplotf(:,e)=sigma_i;
end

%% ��ն����߷�����������
%���ȼ��������Ժ���P0��Ӧ�ĵ��Խ�
totsigma = max(sigma_max); %�ж�sigma_max ��sigmas֮��Ĺ�ϵ
L = totsigma/sigmas;
P0=1/L*P;
n = 10; %�趨ʣ����ص���������
deltP=1/n*(1-1/L)*P;
[sigma0_max,sigma0,tota0,smnode0] = solver( K,P0,D,numnod,numel,nel,bnoderow,outp ,Ixi,Ieta,coord,connect,gaussxi,gausseta);

%% ʩ�ӹ̶���������deltP
%����Ĺؼ����ڱ䵥Ԫ�ն�
% sigma1_max = sigma0_max; %��ʼ������ЧӦ��Ϊ����Ժ���ʱ�ĵ�ֵ
elsigma_max = sigma0_max;
gausstress=sigma0;
Kepll=sparse(2*numnod,2*numnod);
me=zeros(numel,1);
for e = 1:numel
    if (elsigma_max(e)-sigmas)*(abs(elsigma_max(e)-sigmas)>1e-10)<0   %�趨����Ϊ1e-10.
        me(e)=1;
    else
        me(e)=0;
    end
    Ke = zeros(8,8);
    for i_gauss = 1:4 %���ø�˹���ַ���
        [Bgauss,J] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%�������Ԫ��˹�㴦��Ӧ�����
        Dep = e_plasm( sigmas,E,Et,nu,gausstress(e,(3*i_gauss-2:3*i_gauss)));
        Dhat = me(e)*D+(1-me(e))*Dep;
        Ke = Ke+ Bgauss'*Dhat*Bgauss*det(J); %��Ԫ�նȾ���
    end
    totn = connect(e,:); %��Ԫ�ڵ���
    totn2 = zeros(2*length(totn),1);  %��ˮƽ����ֱλ�ƣ���������һ��
    for i =1:nel
        totn2(2*i-1:2*i,1) =  [2*totn(i)-1;2*totn(i)];  %�����Ԫ�ڵ��Ӧ���к�
    end
    Kepll(totn2,totn2) = Kepll(totn2,totn2)+Ke;  %��װ����ն�
    
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
    %����gauss�㴦��Ӧ����Ӧ��
        for i_gauss = 1:4 %���ȡ��˹��
            [Bgauss,~] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%�������Ԫ��˹�㴦��Ӧ�����
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
        smdelsigma(e,:,ip+1)=smdelsigma(e,:,1)+sum(sm); %ÿ����ԪӦ��ĥƽ��Ľڵ�ֵ
        delsigma_i = zeros(nel,1);
        for i =1:nel
            [~,delsigma_i(i)] = sigmainner(smdelsigma(:,:,ip+1), e,Ixi(i),Ieta(i));
        end
        sigplot(:,e)=delsigma_i;
        delsigma_max(e,ip+1) =max(delsigma_i);  %ȡ�õ�Ԫ������ЧӦ����������Ӧ��Ԫ��
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
        for i_gauss = 1:4 %���ø�˹���ַ���
            [Bgauss,J] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%�������Ԫ��˹�㴦��Ӧ�����
            Dep = e_plasm( sigmas,E,Et,nu,sigmatot(e,(3*i_gauss-2:3*i_gauss),ip+1));
            Dhat = m(e,ip+1)*D+(1-m(e,ip+1))*Dep;
            Ke = Ke+ Bgauss'*Dhat*Bgauss*det(J); %��Ԫ�նȾ���
        end
        totnin = connect(e,:); %��Ԫ�ڵ���
        [~,In] = setdiff(totnin,nodeb);
        totn = totnin(sort(In));
        [~,kkill,~] = intersect(totnin,nodeb);
        if isempty(kkill)==0
                kl = [2*kkill-1;2*kkill];
                Ke(kl,:) = [];
                Ke(:,kl) = [];
        end        
        totn2 = zeros(2*length(totn),1);  %��ˮƽ����ֱλ�ƣ���������һ��
        for i =1:length(totn)
            nr = find(rowke==totn(i));
          totn2(2*i-1:2*i,1) =  [2*nr-1;2*nr];  %�����Ԫ�ڵ��Ӧ���к�
        end
        Kep(totn2,totn2,ip+1) = Kep(totn2,totn2,ip+1)+Ke;  %��װ����ն�
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
uep=zeros(numnod,1); %ˮƽλ��
vep=zeros(numnod,1); %����λ��
wholea=zeros(numnod,1); %��λ��
for i = 1:numnod
    uep(i) = totu(2*i-1)*1000;%ȡ��ˮƽλ��
    vep(i) = totu(2*i)*1000; %ȡ������λ��
end

%% ��ͼ
uep = reshape(uep,numnody,numnodx);   %�ع�����λ�ƾ���
vep = reshape(vep,numnody,numnodx); %�ع�����λ�ƾ���
subplot(2,2,2)
surface(coordx,coordy,vep);colorbar;  %����λ�Ʊ���
shading interp
ylim([0,1])    %���������귶Χ
xlim([0,8])
axis equal   %�ᵥλ��ͬ
subplot(2,2,4)
mesh(coordx,coordy,vep)
ylim([0,1])    %���������귶Χ
xlim([0,8])
% axis equal   %�ᵥλ��ͬ



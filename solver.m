function [ sigma0_max, sigma0,tota0,smnode0 ] = solver( K,P0,D,numnod,numel,nel,bnoderow,outp,Ixi,Ieta,coord,connect,gaussxi,gausseta )
%�ú��������������ܸնȣ����غͱ�����ϵʱ�ĸ���Ԫ����ЧӦ��
%   �������ϵı������⣬�����Ϊ���񻮷�ʱ�ĳ�����
%�÷��̵�ȱ�㣺���Ǳ�նȵķ��̡�
a0=K\P0;%����P0����ʱ���Խ�
tota0 = zeros(2*numnod,1);
rowk0=setdiff(1:2*numnod,bnoderow);
tota0(rowk0)=a0;
u0=zeros(numnod,1); %P0���ض�Ӧ��ˮƽλ��
v0=zeros(numnod,1); %P0���ض�Ӧ������λ��
for i = 1:numnod
    u0(i) = tota0(2*i-1);%P0���ض�Ӧȡ��ˮƽλ��
    v0(i) = tota0(2*i); %P0���ض�Ӧȡ������λ��
end
sigma0=zeros(numel,nel*3);
strain0=zeros(numel,nel*3);
for e = 1:numel
    %����gauss�㴦��Ӧ����Ӧ��
    for i_gauss = 1:4 %���ȡ��˹��
        [Bgauss,~] = B(e,coord,connect,gaussxi(i_gauss),gausseta(i_gauss));%�������Ԫ��˹�㴦��Ӧ�����
        [sigma0(e,3*i_gauss-2:3*i_gauss),~]=stress_strain(e,tota0,connect,D,Bgauss);
    end
end
smnode0=sigma0*outp; %�õ��ڵ㴦ĥƽ���Ӧ��ֵ
sigma0_max = zeros(numel,1);
for e =1:numel
    sigma0_i = zeros(nel,1);
    for i =1:nel
        [~,sigma0_i(i)] = sigmainner(smnode0, e,Ixi(i),Ieta(i));
    end
    sigma0_max(e) =max(sigma0_i);  %ȡ�õ�Ԫ������ЧӦ����������Ӧ��Ԫ��
end

end


% clc
% %���Կ�����ЧӦ������ֵ�����ڽڵ㴦
% for e = 1:32
%     for i = 1:10
%         xi = -1+2/10*i;
%         for j = 1:10
%             eta = -1+2/10*i;
%             [~,x]=sigmainner(smnode, e,xi,eta);
%             xx(i,j) = x;
%         end
%     end
%    mesh(1:10,1:10,xx)
%    hold on
% end
% 
pp=zeros(numnod,1);
n=1;
for e =1:numel
    aa=connect(e,:);
    for j=1:nel
        m=aa(j);
        pp(m,n)=sigplotf(j,e);
    end
     n=n+1;
end
for i = 1:numnod
    py(i,1) = sum(pp(i,:))/length(find(pp(i,:)));
end           
sig = reshape(py,numnody,numnodx);   %�ع�����λ�ƾ���
surface(coordx,coordy,sig);colorbar;  %����λ�Ʊ���
shading interp
ylim([0,1])    %���������귶Χ
xlim([0,8])
axis equal


for ani=1:4
    xint=fix(numelx/5*ani+1);
    yint=fix(numely/2+1);
    rateright=numelx/5*ani+1-xint;
    vpoint(ani)=(1-rateright)*v(yint,xint)+rateright*v(yint,xint+1);
    veppoint(ani)=(1-rateright)*vep(yint,xint)+rateright*vep(yint,xint+1);
end
    
    
% clc
% for e = 1:32
%     for i = 1:10
%         xi = -1+2/10*i;
%         for j = 1:10
%             eta = -1+2/10*i;
%             [~,x]=sigmainner(smnode, e,xi,eta);
%             xx(i,j) = x;
%         end
%     end
%     row = (rem(e,2)+1)*10-9:1:(rem(e,2)+1)*10;
%     colum = ceil(e/2)*10-9:1:ceil(e/2)*10;
%     ck(row,colum) = xx
% end
% surface(1:160,1:20,ck)
% axis equal
% atest = [1,2,3;1,2,3];
% btest=[1,1]';
% s=zeros(2,3)
% for is=1:3
%     s(:,is) = atest(:,is);
%             
% end
% % c =sum(s,2)+btest
% 
% atest=[1;1];
% btest=[1,2;1,2];
% for i =1: 2
%    atest = atest + btest(:,i)
% end
pp=zeros(numnod,1);
n=1;
for e =1:numel
    aa=connect(e,:);
    for j=1:nel
        m=aa(j);
        pp(m,n)=sigplot(j,e);
    end
     n=n+1;
end
for i = 1:numnod
    py(i,1) = sum(pp(i,:))/length(find(pp(i,:)));
end       
sig = reshape(py,numnody,numnodx);   %重构横向位移矩阵
surface(coordx,coordy,sig);colorbar;  %绘制位移表面
shading interp
ylim([0,1])    %控制轴坐标范围
xlim([0,8])
axis equal

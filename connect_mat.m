function [connect_mat,A] = connect_mat( numnodx,numnody,nel)
%输入横纵坐标的节点数目，和单元自由度
%输出连接矩阵，每个单元涉及的节点的编号
xn = 1:(numnodx*numnody);%拉成一条编号
A = reshape(xn,numnody,numnodx);%同形状编号
B = zeros(numnody,numnodx);
for i = 1:numnody
B(i,:)=A(numnody+1-i,:);%和坐标轴一直，向右和向上为正。
end
A=B;
for i = 1:(numnodx-1)*(numnody-1)
    yg = rem(i,numnody-1);%xg表示单元为下边界数起第几个
    if yg == 0
        yg = numnody-1;
    end
    xg = ceil(i/(numnody-1));%左边界其数第几个
    a = A(numnody-yg:numnody-yg+1,xg:xg+1);%这个小矩阵，拉直了就是连接矩阵
    a_veg = a(:);
    connect_mat(i,1:nel) =a_veg([2 4 3 1]);   %逆时针编码，与标准四节点单元编码一致，左下为1
end
end
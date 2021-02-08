function [ sigmainner,sigmahat ] = sigmainner( smnode, e,xi,eta )
%单元内部得到任意位置的应力值和该点的等效应力值
%计算单元内部的应力，用节点处磨平后的应力值内插   smnode 是磨平以后的单元节点值；e 是单元数，xi和eta是单元等参坐标
% sigmainner 是一个一行三列的矩阵。分别对应：sigmainnerx sigmainnery tauxy
%sigmahat 是等效应力值 判断是否进入塑性，及判断sigmahat 与sigmas 之间的关系
N1 = 1/4*(1-xi)*(1-eta);
N2 = 1/4*(1+xi)*(1-eta);
N3 = 1/4*(1+xi)*(1+eta);
N4 = 1/4*(1-xi)*(1+eta);
N = [eye(3)*N1;eye(3)*N2;eye(3)*N3;eye(3)*N4];
sigmainner = smnode(e,:)*N;
sigmam=1/3*(sigmainner(1)+sigmainner(2));
Sx = sigmainner(1)-sigmam;
Sy = sigmainner(2)-sigmam;
Sz = 0-sigmam;
tauxy=sigmainner(3);
sigmahat = sqrt(3)*sqrt(1/2*(Sx^2+Sy^2+Sz^2)+tauxy^2);
end


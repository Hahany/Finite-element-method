function [ D ] = e_plasm( sigmas,E,Et,nu,sigma )
%该函数用来生成弹塑性矩阵函数，用于增量法的计算
% 第一步：构建应力偏量,函数需要输入的是应力向量sigma=（sigmax,sigmay,tauxy）
sigmam=1/3*(sigma(1)+sigma(2));
Sx = sigma(1)-sigmam;
Sy = sigma(2)-sigmam;
tauxy=sigma(3);
G=E/(2*(1+nu));
A=E*Et/(E-Et);
B = Sx^2+Sy^2+2*nu*Sx*Sy+2*(1-nu)*tauxy^2+2*(1-nu)*A*sigmas^2/(9*G);
Dp = E/(B*(1-nu^2))*[(Sx+nu*Sy)^2,(Sx+nu*Sy)*(Sy+nu*Sx),(1-nu)*(Sx+nu*Sy)*tauxy;
    (Sx+nu*Sy)*(Sy+nu*Sx),(Sy+nu*Sx)^2,(1-nu)*(Sy+nu*Sx)*tauxy;
    (1-nu)*(Sx+nu*Sy)*tauxy,(1-nu)*(Sy+nu*Sx)*tauxy,(1-nu)^2*tauxy^2];
D = elasm(E,nu) - Dp ;
return
end


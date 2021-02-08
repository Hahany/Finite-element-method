function [ sigmainner,sigmahat ] = sigmainner( smnode, e,xi,eta )
%��Ԫ�ڲ��õ�����λ�õ�Ӧ��ֵ�͸õ�ĵ�ЧӦ��ֵ
%���㵥Ԫ�ڲ���Ӧ�����ýڵ㴦ĥƽ���Ӧ��ֵ�ڲ�   smnode ��ĥƽ�Ժ�ĵ�Ԫ�ڵ�ֵ��e �ǵ�Ԫ����xi��eta�ǵ�Ԫ�Ȳ�����
% sigmainner ��һ��һ�����еľ��󡣷ֱ��Ӧ��sigmainnerx sigmainnery tauxy
%sigmahat �ǵ�ЧӦ��ֵ �ж��Ƿ�������ԣ����ж�sigmahat ��sigmas ֮��Ĺ�ϵ
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


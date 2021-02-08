function [ estress,estrain ] = stress_strain(e, tota,connect,D,B)
%该函数用来计算应力值和应变值
% 计算任意位置的应力和应变，xi和eta坐标任意，由B的值来体现
nodes = connect(e,:);
unodes = zeros(2*length(nodes),1);
for i = 1:length(nodes)
    unodes(2*i-1) = 2*nodes(i)-1;
    unodes(2*i) = 2*nodes(i);
end
ae = tota(unodes);
estress = D*B*ae;
estrain = B*ae;
end


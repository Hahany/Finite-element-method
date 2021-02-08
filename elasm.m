function [ De ] = elasm(E,nu  )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
De=E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
return
end


function [ De ] = elasm(E,nu  )
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
De=E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
return
end


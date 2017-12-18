function [c, ceq]=g_func_pinhole_mb(x)
%Wrapper function to calculate the beam width after pinhole and perform
%paramter scaling.
B=x(1)*1e-5;
P=x(2)*1e5;
d=x(3)*1e-6;
c=[];
ceq = 1e9*g_func_pinhole(B,P,d);

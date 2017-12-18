function [c, ceq]=g_func_mb(x)
%Wrapper function to calculate the beam width after a zone plate and perform
%paramter scaling.
B=x(1)*1e-5;
P=x(2)*1e5;
c=[];
ceq = 1e9*g_func(B,P);

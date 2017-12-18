function fun=I_func_mb(x)
%Wrapper function to calculate the intensity after a zone plate and perform
%paramter scaling.
B=x(1)*1e-5;
P=x(2)*1e5;
fun = I_func(B,P);
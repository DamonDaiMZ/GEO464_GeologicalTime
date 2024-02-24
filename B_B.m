clear; clc;
a = 0.89123;
sig_a = 0.00012;
b = 132.34;
sig_b = b*0.01;
c = 0.70721;
lambda87 = 1.42e-11;
sig_lambda87 = lambda87.*2e-2;
t = log((a-c)./b + 1)./lambda87
disp(t/1e6)
dtda = 1./lambda87./(a-c+b);
dtdb = (c-a)./lambda87./((a-c).*b + b.^2);
dtdlamb = -log((a-c)./b + 1)./lambda87^2;
sig_t = sqrt( (dtda.*sig_a).^2 + (dtdb.*sig_b).^2 )
disp(sig_t/1e6)


sig_tt = sqrt( (dtda.*sig_a).^2 + (dtdb.*sig_b).^2 + (dtdlamb.*sig_lambda87).^2)
disp(sig_tt/1e6)

d = 0.015495;
sig_d = d.*0.02e-2;
lambda238 = 1.55125e-10;
sig_lambda238 = lambda238.*0.107e-2;

t2 = log(d+1)./lambda238
disp(t2/1e6)
dt2dd = 1./(lambda238*(d+1));
dt2dlam = -log(d+1)./lambda238.^2;
sig_t2 = dt2dd.*sig_d;
disp(sig_t2/1e6)

sig_tt2 = sqrt((dt2dd.*sig_d).^2 + (dt2dlam.*sig_lambda238).^2)
disp(sig_tt2/1e6);
data=importdata('/data/kabir/output/ppf/old_stuff/RSnostim_750_leak_noPMCA/s_00001/dat/pmca&leak_ca_flux.dat');
%data=importdata('/data/kabir/output/ppf/RSnostim_750_leak_noPMCA_correction/s_00001/dat/pmca&leak_ca_flux.dat');

A=data.data;
Y=A(:,4);
X=A(:,1);
clf
plot(X,Y)
hold on
P=polyfit(X(1:100),Y(1:100),1);
x=[0:0.01:0.1];
plot(x,P(1)*x+P(2));

leak_rate=P(1)
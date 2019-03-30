
data=[importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.2/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.4/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.6/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.8/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.2/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.4/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.6/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.8/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long2/s_00001/dat/ca.dat')];


P=zeros(9,2);
for k=1:9
    A=data(k).data;
    Y=A(:,4);
    X=A(:,1);
    P(k,:)=polyfit(X,Y,1);
end
P
A=data(9).data;
Y=A(:,4);
X=A(:,1);
p=polyfit(X,Y,1);
clf
scatter(X,Y)
hold on
x=[0:0.01:0.5];
plot(x,p(1)*x+p(2))
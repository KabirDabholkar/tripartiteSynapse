%range analysis

data=[importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/0.2/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/0.4/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/0.6/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/0.8/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/1/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/1.5/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/1.8/s_00001/dat/ca.dat');
importdata('/data/kabir/output/ppf/two_way_leak_range/RSnostim_750_leak_noPMCA/2/s_00001/dat/ca.dat');];

% data=[importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.2/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.4/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.6/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.8/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.2/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.4/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.6/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long1.8/s_00001/dat/serca_ca_flux.dat');
% importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long2/s_00001/dat/serca_ca_flux.dat')];



P=zeros(8,2);
end_point=zeros(8,1);
for k=1:8
    A=data(k).data;
    Y=A(:,4);%A(:,3)-A(:,2);
    X=A(:,1);
    P(k,:)=polyfit(X,Y,1);
    end_point(k)=Y(end);
end
P
end_point
% A=data(5).data;
% Y=A(:,3)-A(:,2);
% X=A(:,1);
% p=polyfit(X,Y,1);
% clf
% scatter(X,Y)
% hold on
% x=[0:0.01:0.5];
% plot(x,p(1)*x+p(2))
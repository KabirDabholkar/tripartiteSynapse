%plotting
d=importdata('/data/kabir/output/ppf/old_stuff/RSnostim_750_long (copy)/s_00001/dat/ca.dat');
d=importdata('/data/kabir/output/ppf/range/RSnostim_750_leak_noPMCA_correction_long0.2/s_00001/dat/ca.dat');
A=d.data;

N_avo=6.0221409e23;
vol_er = 3.9*0.1*0.1;
plot(A(:,1)*1e3,A(:,4)*1e15/N_avo/vol_er*1e6,'LineWidth',3)
ylabel("Ca2+ concentration in ER (uM)")
xlabel("Time (ms)")
axis([0 500 0 1000])
title("SERCA* + Leak")
saveas(gcf,'correction_0.2','epsc')
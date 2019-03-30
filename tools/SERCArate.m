sm=[1 2 3 5 10];

% d1=importdata('/data/kabir/output/ppf/RSnostim_750_emptyER_sm/sm1/s_00001/dat/ca.dat').data;
% d2=importdata('/data/kabir/output/ppf/RSnostim_750_emptyER_sm/sm2/s_00001/dat/ca.dat').data;
% d3=importdata('/data/kabir/output/ppf/RSnostim_750_emptyER_sm/sm3/s_00001/dat/ca.dat').data;
% d5=importdata('/data/kabir/output/ppf/RSnostim_750_emptyER_sm/sm5/s_00001/dat/ca.dat').data;
% d10=importdata('/data/kabir/output/ppf/RSnostim_750_emptyER_sm/sm10/s_00001/dat/ca.dat').data;

for k=sm
    dat=importdata(strcat('/data/kabir/output/ppf/RSnostim_750_emptyER_sm/sm',string(k),'/s_00001/dat/ca.dat'));
    d=dat.data;
    plot(d(:,1),d(:,4))
    hold on
end
legend(string(sm));

curve

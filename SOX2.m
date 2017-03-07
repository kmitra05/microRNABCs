%% Upload bcmir as 3 column vector of all deltaCT values across conditions. Upload SOX2gap as 3 column vector of
%% of deltaCT values
dbcmir=2.^(-bcmir);
dsox2gap=2.^(-SOX2gap);
for i=1:3
    norbcmir(:,i)=dbcmir(:,i)./dbcmir(6,i);
    norsox2gap(:,i)=dsox2gap(:,i)./dsox2gap(6,i);
end
meanbcmir=mean(norbcmir,2);
meansox2gap=mean(norsox2gap,2);
stdmir=std(norbcmir,0,2);
stdsox=std(norsox2gap,0,2);
h=errorbar(meanbcmir,meansox2gap,stdsox,stdsox,stdmir,stdmir,'*r')
hold on
h.Parent.XScale='log';
h.Parent.YScale='log';
xlabel('normalized microRNA signal (a.u)');
ylabel('normalized SOX2 mRNA signal (a.u)');
title('SOX2mRNA vs miRNA signal')
fit=polyfit(meanbcmir,meansox2gap,1);
x1=linspace(1,h.Parent.XLim(2));
f1 = polyval(fit,x1);
plot(x1,f1,'b--')
Rsq1 = 1 - sum((meansox2gap - (meanbcmir*fit(1))).^2)/sum((meansox2gap - mean(meansox2gap)).^2);
s1=sprintf('R^2 value = %2f',Rsq1);
text(10,700,s1)
hold off
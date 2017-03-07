%% Upload lucCT as 3 column vector of all deltaCT values across conditions. Upload luc as 2 column vector of
%% 2 runs of luciferase 
dlucCT=2.^-(lucCT);
for i=1:3
    normlucCT(:,i)=dlucCT(:,i)./dlucCT(9,i);
end
for i=1:2
    normluc(:,i)=luc(:,i)./luc(9,i);
end
meanluc=mean(normluc,2);
meanlucCT=mean(normlucCT,2);
stdluc=std(normluc,0,2);
stdlucCT=std(normlucCT,0,2);
figure
h=errorbar(meanlucCT,meanluc,stdluc,stdluc,stdlucCT,stdlucCT,'*r')
hold on
h.Parent.XScale='log';
h.Parent.YScale='log';
xlabel('normalized microRNA signal (a.u)');
ylabel('normalized luciferase signal (a.u)');
title('Luciferase vs miRNA signal')
fit=polyfit(meanlucCT,meanluc,1);
x1=linspace(2,h.Parent.XLim(2));
f1 = polyval(fit,x1);
plot(x1,f1,'b--')
Rsq1 = 1 - sum((meanluc - (meanlucCT*fit(1))).^2)/sum((meanluc - mean(meanluc)).^2);
s1=sprintf('R^2 value = %2f',Rsq1);
text(80,5000,s1)
hold off
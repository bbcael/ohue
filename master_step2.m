clear all; close all; clc; load T_H_Uh.mat;

%% show kappa has curvature
%{
for i = 1:200;
    Ti = cumsum(T(:,i));
    fit2 = fit(Ti,H,'poly2','Weights',Uh.^-2);
    P1(i) = fit2.p1;
end
sum(P1>0)./length(P1); % p = 0.985

[y x] = ecdf(P1);
plot(x,y,'linewidth',2,'color',[128 40 97]./256)
box on
set(gca,'fontsize',16,'ticklabelinterpreter','latex')
ylabel('CDF','interpreter','latex')
xlabel('Quadratic term of quadratic fit','interpreter','latex')
lgnd = legend('$p(x>0) = 0.985$'); set(lgnd,'interpreter','latex','location','northwest')
%}
%% find delta

d = -0.002:.0001:.011;

for j = 1:size(T,2);
    j
for i = 1:length(d);
    delta = d(i);
    tvar = 1+delta.*(0:51);
    Ti = cumsum(tvar'.*T(:,j));
    fit3 = fit(Ti,H,'poly3','Weights',Uh.^-2);
    ci3 = confint(fit3);
    z3(i) = abs(fit3.p1)./abs(ci3(2,1)-fit3.p1).*1.96; % z cubic
    fit2 = fit(Ti,H,'poly2','Weights',Uh.^-2);
    ci2 = confint(fit2);
    z2(i) = abs(fit2.p1)./abs(ci2(2,1)-fit2.p1).*1.96; % z cubic
    [fit1,gof] = fit(Ti,H,'poly1');
    k(i) = fit1.p1;
    rss(i) = gof.sse;
end
rss(z2>1.645) = max(rss);
rss(z3>1.645) = max(rss);
[~,indx] = min(rss);
D(j) = d(indx);
K(j) = k(indx);
Z2(j) = z2(indx);
Z3(j) = z3(indx);
end
clear all; close all; clc; % generate ensembles
load domlev_ishii_cheng.mat; % mat file with data from Cheng website
H = H-H(1,:); H = H-H(end,:)+mean(H(end,:)); % make in agreement at end as integrated backwards
C = cov(H');
M = mean(H');
h = mvnrnd(M,C,10000);
load T.mat; % HadCRUT5 temperature anomalies as .mat file
C = cov(T(1:48,:)');
M = mean(T(1:48,:)');
t = mvnrnd(M,C,10000);
h = h-h(:,1); % make in agreement in beginning after the fact as defined as anomaly w/r/t 1970
clearvars -EXCEPT h t; save ht.mat;

%%

clear all; close all; clc; load ht.mat; % curvature in residuals in null case

for i = 1:10000;
    Ti = cumsum(t(i,:));
    fit2 = fit(Ti',h(i,:)','poly2');
    P1(i) = fit2.p1;
end

%%

clear all; close all; clc; load ht.mat; % main calculation

d = -0.01:.0001:.025;

for j = 1:size(t,1);
    j
for i = 1:length(d);
    delta = d(i);
    tvar = 1+delta.*(0:47);
    Ti(i,:) = cumsum(tvar.*t(j,:));
end
H = repmat(h(j,:),size(Ti,1),1);
[r,m,b] = regression(Ti,H);
[~,indx] = max(r);
D(j) = d(indx);
K(j) = m(indx);
end

save multi_DK.mat;

%%

clear all; close all; clc; load ht.mat; % curvature in residuals in delta=/=0 case

d = -0.01:.0001:.025;

for j = 1:size(t,1);
    j
for i = 1:length(d);
    delta = d(i);
    tvar = 1+delta.*(0:47);
    Ti = cumsum(tvar.*t(j,:));
    fit2 = fit(Ti',h(j,:)','poly2');
    ci2 = confint(fit2);
    z2(i) = abs(fit2.p1)./abs(ci2(2,1)-fit2.p1).*1.96; % z cubic
    [fit1,gof] = fit(Ti',h(j,:)','poly1');
    rss(i) = gof.sse;
end
[~,indx] = min(rss);
Z2(j) = z2(indx);
end

%%

clear all; close all; clc; load ht.mat; % t-dep. vs. T-dep. calculation

d = -0.01:.001:.025;

for j = 1:size(t,1);
    j
for i = 1:length(d);
    delta = d(i);
    tvar = 1+delta.*(0:47);
    Ti(i,:) = cumsum(tvar.*t(j,:));
    Tvar = 1+(1/40).*delta.*t(j,:); % can adjust conversion factor for range of delta-T vs. delta-t sampled but this does not make a difference 
    TiT(i,:) = cumsum(Tvar.*t(j,:));
end
H = repmat(h(j,:),size(Ti,1),1);
[r,~,~] = regression(Ti,H);
[r2,~,~] = regression(TiT,H);
Rt(j) = max(r);
RT(j) = max(r2);
end

sum(RT<Rt)./length(RT)

%%

clear all; close all; clc; % figure 3
load multi_DK.mat;
t = 0:48;
kappa = K'./16.0886.*(1+D'.*t);
for i = 1:48;
    k(:,i) = prctile(kappa(:,i),[16 50 84]);
end

A1 = k(1,:)';
A2 = (k(2,:)-k(1,:))';
A3 = (k(3,:)-k(2,:))';
A1 = linspace(A1(1),A1(end),length(A1));
A2 = linspace(A2(1),A2(end),length(A2));
A3 = linspace(A3(1),A3(end),length(A3));

figure;
a = area(1970:2017,[A1; A2; A3]');
a(1).FaceColor = [1 1 1];
a(1).EdgeColor = [1 1 1];
a(2).EdgeColor = [1 1 1];
a(3).EdgeColor = [1 1 1];
a(2).FaceColor = [128 40 97]./256;
a(3).FaceColor = [128 40 97]./256;
a(2).FaceAlpha = 0.6;
a(3).FaceAlpha = 0.6;
hold on;
plot(1970:2017,A1+A2,'linewidth',2,'color',[128 40 97]./256)
set(gca,'ticklabelinterpreter','latex','fontsize',16)
ylabel('Ocean heat uptake efficiency $\kappa$ [W/m$^{2}$K]','interpreter','latex')
axis([1969.9 2017.1 .5 .76])
xlabel('Year','interpreter','latex')
lgnd = legend([a(2)],'$\pm1$ s.d.')
set(lgnd,'interpreter','latex','fontsize',16,'location','northwest')
%%

clear all; close all; clc; load multi_DK.mat; load ht.mat; % figure 1

T = cumsum((1+(0:47)'.*D)'.*t,2);
Tm = median(T);
Hm = median(h);
plot(0:34,9.334.*(0:34)+18.24,'k','linewidth',3)
hold on;
scatter(Tm,Hm,100,1970:2017,'filled');
scatter(Tm,Hm,105,1970:2017,'k');
clrbr = colorbar;
set(clrbr,'fontsize',16,'ticklabelinterpreter','latex','ytick',1970:10:2020)
ylabel(clrbr,'Year','interpreter','latex','fontsize',16)
turbomap;
colormap turbo;
xlabel('$\mathcal{T}_\delta$ [K y]','interpreter','latex','fontsize',16)
ylabel('$\mathcal{H}$ [ZJ]','interpreter','latex','fontsize',16)
axis([-1 36 -20 350])
set(gca,'fontsize',16,'ticklabelinterpreter','latex')
text(0,325,'$\kappa_{1970} = 0.58 \pm 0.08$ W/m$^{2}$K','interpreter','latex','fontsize',16)
text(3,300,'$\delta = 0.46 \pm 0.43$ \% y$^{-1}$','interpreter','latex','fontsize',16)

%%

clear all; close all; clc; load multi_DK.mat; %figure 2

scatterhist(K(D<max(D))./16.0886,47.*K(D<max(D))./16.0886.*D(D<max(D)),'Kernel','on','Color',[160 79 56]./256,'linewidth',2,'Marker','+','MarkerSize',[3],'direction','out','location','northeast')
set(gca,'ticklabelinterpreter','latex','fontsize',16,'xtick',.4:.2:.8,'ytick',-.2:.2:.4)
ylabel('$\kappa_{2017}-\kappa_{1970}$ [W/m$^{2}$K]','fontsize',16','interpreter','latex')
xlabel('$\kappa_{1970}$ [W/m$^{2}$K]','fontsize',16','interpreter','latex')
axis([.3 .95 -.3 .5])

%%

clear all; close all; clc; load multi_DK.mat; % 1.5/2 degrees calculation
f2x = 4+.3.*randn(1,10000);
l = 1.3+.44.*randn(1,10000);
tcr0 = f2x./(l+K./16.0886);
tcr1 = f2x./(l+K.*(1+51.*D)./16.0886);
y_1p5_0 = 70.*1.5./tcr0;
y_2_0 = 70.*2./tcr0;
y_1p5_1 = 70.*1.5./tcr1;
y_2_1 = 70.*2./tcr1;

ecdf(y_2_1-y_2_0)
hold on;
ecdf(y_1p5_1-y_1p5_0)

[median(y_2_1-y_2_0) std(y_2_1-y_2_0) median(y_1p5_1-y_1p5_0) std(y_1p5_1-y_1p5_0)]

%%

clear all; close all; clc; load ht.mat; % figure 4
T = cumsum(t,2);
[r,m,b] = regression(T,h);
R = h-(m.*T+b);
Rm = mean(R);
Ru = .5.*diff(prctile(R,[84 16]));
errorbar(1970:2017,Rm,Ru,'linewidth',2,'color',[85 102 80]./256)
set(gca,'ticklabelinterpreter','latex','fontsize',16)
ylabel('Residuals [ZJ]','interpreter','latex')
axis([1969 2018 -40 40])
hold on;
plot(1970:2017,(1:48).^2*.009377-.4785.*(1:48)+4.296,'k','linewidth',2)
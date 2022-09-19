%% generate ensembles
clear all; close all; clc; load HT.mat;
n = 200;
H = [Hc Hi Hn];
C = cov(H');
M = mean(H');
h = mvnrnd(M,C,n);
C = cov(T');
M = mean(T');
t = mvnrnd(M,C,n);
t = repmat(T,1,n./200)';
clearvars -EXCEPT h t;

%% calculate curvature when delta = 0

for i = 1:size(t,1);
    Ti = cumsum(t(i,:));
    fit2 = fit(Ti',h(i,:)','poly2');
    P1(i) = fit2.p1;
    i
end

%% find delta and kappa

d = -0.005:.0001:.02;
for j = 1:size(t,1);
    j
for i = 1:length(d);
    delta = d(i);
    tvar = 1+delta.*(0:size(t,2)-1);
    Ti(i,:) = cumsum(tvar.*t(j,:));
end
H = repmat(h(j,:),size(Ti,1),1);
[r,m,b] = regression(Ti,H);
[~,indx] = max(r);
ind0 = find(d==0);
D(j) = d(indx);
K(j) = m(indx)./16.09;
end

%% figure 1

figure;
T = cumsum((1+(0:49)'.*D)'.*t,2);
Tm = median(T);
Hm = median(h);
[~,m,b] = regression(Tm,Hm-Hm(1));
plot(0:43,m.*(0:43)+b,'k','linewidth',3)
hold on;
scatter(Tm,Hm-Hm(1),100,1970:2019,'filled');
scatter(Tm,Hm-Hm(1),105,1970:2019,'k');
clrbr = colorbar;
set(clrbr,'fontsize',16,'ticklabelinterpreter','latex','ytick',1970:10:2020)
ylabel(clrbr,'Year','interpreter','latex','fontsize',16)
turbomap;
colormap turbo;
xlabel('$\mathcal{T}_\delta$ [K y]','interpreter','latex','fontsize',16)
ylabel('$\mathcal{H}$ [ZJ]','interpreter','latex','fontsize',16)
axis([-1 45 -20 370])
set(gca,'fontsize',16,'ticklabelinterpreter','latex')
text(0,350,'$\kappa_{1970} = 0.48 \pm 0.08$ W/m$^{2}$K','interpreter','latex','fontsize',16)
text(3,325,'$\delta = 0.83 \pm 0.27$ \% y$^{-1}$','interpreter','latex','fontsize',16)


%% T-dep. vs. t-dep. calculation

d = -0.005:.0001:.02;

for j = 1:size(t,1);
    j
for i = 1:length(d);
    delta = d(i);
    tvar = 1+delta.*(0:size(t,2)-1);
    Ti(i,:) = cumsum(tvar.*t(j,:));
    Tvar = 1+50*delta.*t(j,:); % can adjust conversion factor for range of delta-T vs. delta-t sampled but this does not make a difference 
    TiT(i,:) = cumsum(Tvar.*t(j,:));
end
H = repmat(h(j,:),size(Ti,1),1);
[r,~,~] = regression(Ti,H);
[r2,~,~] = regression(TiT,H);
[Rt(j),it] = max(r);
[RT(j),iT] = max(r2);
dt(j) = d(it);
dT(j) = d(iT);
end

sum(RT<Rt)./length(RT)

%% 1.5/2C calculation

f2x = 4+.3.*randn(1,length(D));
l = 1.3+.44.*randn(1,length(D));
tcr0 = f2x./(l+K);
tcr1 = f2x./(l+K.*(1+49.*D));
y_1p5_0 = 70.*1.5./tcr0;
y_2_0 = 70.*2./tcr0;
y_1p5_1 = 70.*1.5./tcr1;
y_2_1 = 70.*2./tcr1;

ecdf(y_2_1-y_2_0)
hold on;
ecdf(y_1p5_1-y_1p5_0)

[median(y_2_1-y_2_0) std(y_2_1-y_2_0) median(y_1p5_1-y_1p5_0) std(y_1p5_1-y_1p5_0)]

%% figure 2

figure;
scatterhist(K.*(1+24.5.*D),49.*K.*D,'Kernel','on','Color',[160 79 56]./256,'linewidth',2,'Marker','+','MarkerSize',[3],'direction','out','location','northeast')
set(gca,'ticklabelinterpreter','latex','fontsize',16,'xtick',.4:.2:.8,'ytick',0:.1:.3)
ylabel('Change in OHUE, $\Delta\kappa$ [W/m$^{2}$K]','fontsize',16','interpreter','latex')
xlabel('Mean OHUE, $\overline\kappa$ [W/m$^{2}$K]','fontsize',16','interpreter','latex')
axis([.3 .95 -.05 .35])

%% figure 3

t = 0:50;
kappa = K'.*(1+D'.*t);
for i = 1:50;
    k(:,i) = prctile(kappa(:,i),[16 50 84]);
end

A1 = k(1,:)';
A2 = (k(2,:)-k(1,:))';
A3 = (k(3,:)-k(2,:))';
A1 = linspace(A1(1),A1(end),length(A1));
A2 = linspace(A2(1),A2(end),length(A2));
A3 = linspace(A3(1),A3(end),length(A3));

figure;
a = area(1970:2019,[A1; A2; A3]');
a(1).FaceColor = [1 1 1];
a(1).EdgeColor = [1 1 1];
a(2).EdgeColor = [1 1 1];
a(3).EdgeColor = [1 1 1];
a(2).FaceColor = [128 40 97]./256;
a(3).FaceColor = [128 40 97]./256;
a(2).FaceAlpha = 0.6;
a(3).FaceAlpha = 0.6;
hold on;
plot(1970:2019,A1+A2,'linewidth',2,'color',[128 40 97]./256)
set(gca,'ticklabelinterpreter','latex','fontsize',16)
ylabel('Ocean heat uptake efficiency $\kappa$ [W/m$^{2}$K]','interpreter','latex')
axis([1969.9 2019.1 .3 .8])
xlabel('Year','interpreter','latex')
lgnd = legend([a(2)],'$\pm1$ s.d.')
set(lgnd,'interpreter','latex','fontsize',16,'location','northwest')

%% figure 4

figure;
[y x] = ecdf(P1);
plot(x,y,'linewidth',2,'color',[160 79 56]./256)
box on
hold on;
set(gca,'fontsize',16,'ticklabelinterpreter','latex')
ylabel('Cumulative Distribution Function','interpreter','latex')
xlabel('Quadratic term of fit to $\delta = 0$ case','interpreter','latex')
axis([-.01 .06 0 1]);
hold on;
plot(-1:1,26/10000+0*(-1:1),'k','linewidth',.75)
plot(0*(0:1),0:1,'k','linewidth',.75)


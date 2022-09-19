%% generate ensembles
clear all; close all; clc; load HT.mat;
n = 10000;
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

%% find delta and kappa with/without delta

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
K0(j) = m(ind0)./16.09;
end

%% T-dep. vs. t-dep. calculation

d = -0.05:.0001:.02;

for j = 1:size(t,1);
    j
for i = 1:length(d);
    delta = d(i);
    tvar = 1+delta.*(0:size(t,2)-1);
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

%% 1.5/2C calculation

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

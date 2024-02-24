clc; clear;

load('Kistler.mat');
% %% remove biotites
% HD = HD(2:end,:);
% CP = CP(2:end,:);
% GA = GA(2:end,:);
%% remove Kspar
HD = [HD(1:2,:);HD(4:5,:)];
CP = [CP(1:3,:);CP(5:7,:)];
GA = [GA(1:2,:);GA(4:6,:)];


%% all data 
alldata = [HD;CP;GA];
x = alldata(:,1);
x_uc = alldata(:,2);
y = alldata(:,3);
y_uc = alldata(:,4);
[a_ALL,a_uc_ALL,b_ALL,b_uc_ALL,MSWD_ALL] = YorkFit(x,y,x_uc,y_uc);
% critical value
CriVal = chi2inv(0.95,length(x)-2)./(length(x)-2);

% calculate the age and uncertainty
lambda = 1.42e-11;
lambda_uc = lambda.*2e-2;
[t_ALL, t_uc_ALL] = IsochronAge(b_ALL, b_uc_ALL, lambda, lambda_uc);

% plot
f1 = figure(1);
f1.Position = [100 100 500 400];
x_lin = linspace(max(alldata(:,1)).*1.1,min(alldata(:,1)).*0.9,100);
y_lin = b_ALL.*x_lin + a_ALL;
e = errorbar(alldata(:,1),alldata(:,3),alldata(:,4),alldata(:,4),alldata(:,2),alldata(:,2),".","MarkerSize",20,"LineStyle","none"); hold on;
plot(x_lin,y_lin,'k--','LineWidth',1.5); hold on;
xlabel('(^{87}Rb/ ^{86}Sr)');
ylabel('(^{87}Sr/ ^{86}Sr)');
txt1 = ['MSWD = ' num2str(MSWD_ALL)];
txt2 = ['y = ' num2str(b_ALL) '±' num2str(b_uc_ALL) 'x + ' num2str(a_ALL) '±' num2str(a_uc_ALL)];
txt3 = ['Age = ' num2str(t_ALL./1e6) '±' num2str(t_uc_ALL./1e6) 'Ma'];
txt4 = ['Critical Value = ' num2str(CriVal) ' (p = 0.95)'];
text(15,0.88,txt2,'FontSize',14);
text(15,0.87,txt1,'FontSize',14);
text(15,0.86,txt4,'FontSize',14);
text(100,0.71,txt3,'FontSize',14);
title('All data isochrom (no K-feldspar)');
set(gca,'FontSize',14,'linewidth',1);


%% HD 
data = HD;
x = data(:,1);
x_uc = data(:,2);
y = data(:,3);
y_uc = data(:,4);
[a_HD,a_uc_HD,b,b_uc,MSWD] = YorkFit(x,y,x_uc,y_uc);
% critical value
CriVal = chi2inv(0.05,length(x)-2)./(length(x)-2);

% calculate the age and uncertainty
lambda = 1.42e-11;
lambda_uc = lambda.*2e-2;
[t, t_uc] = IsochronAge(b, b_uc, lambda, lambda_uc);

% % plot
% f2 = figure(2);
% f2.Position = [100 100 500 400];
% x_lin = linspace(max(data(:,1)).*1.1,min(data(:,1)).*0.9,100);
% y_lin = b.*x_lin + a;
% e = errorbar(data(:,1),data(:,3),data(:,4),data(:,4),data(:,2),data(:,2),".","MarkerSize",20,"LineStyle","none"); hold on;
% plot(x_lin,y_lin,'k--','LineWidth',1.5); hold on;
% xlabel('(^{87}Rb/ ^{86}Sr)');
% ylabel('(^{87}Sr/ ^{86}Sr)');
% txt1 = ['MSWD = ' num2str(MSWD)];
% txt2 = ['y = ' num2str(b) '±' num2str(b_uc) 'x + ' num2str(a) '±' num2str(a_uc)];
% txt3 = ['Age = ' num2str(t./1e6) '±' num2str(t_uc./1e6) 'Ma'];
% txt4 = ['Critical Value = ' num2str(CriVal) ' (p = 0.05)'];
% text(15,0.87,txt2,'FontSize',14);
% text(15,0.86,txt1,'FontSize',14);
% text(15,0.85,txt4,'FontSize',14);
% text(50,0.71,txt3,'FontSize',14);
% title('HD series (no K-feldspar)');
% set(gca,'FontSize',14,'linewidth',1);

%% CP
data = CP;
x = data(:,1);
x_uc = data(:,2);
y = data(:,3);
y_uc = data(:,4);
[a_CP,a_uc_CP,b,b_uc,MSWD] = YorkFit(x,y,x_uc,y_uc);
% critical value
CriVal = chi2inv(0.05,length(x)-2)./(length(x)-2);

% calculate the age and uncertainty
lambda = 1.42e-11;
lambda_uc = lambda.*2e-2;
[t, t_uc] = IsochronAge(b, b_uc, lambda, lambda_uc);

% % plot
% f3 = figure(3);
% f3.Position = [100 100 500 400];
% x_lin = linspace(max(data(:,1)).*1.1,min(data(:,1)).*0.9,100);
% y_lin = b.*x_lin + a;
% e = errorbar(data(:,1),data(:,3),data(:,4),data(:,4),data(:,2),data(:,2),".","MarkerSize",20,"LineStyle","none"); hold on;
% plot(x_lin,y_lin,'k--','LineWidth',1.5); hold on;
% xlabel('(^{87}Rb/ ^{86}Sr)');
% ylabel('(^{87}Sr/ ^{86}Sr)');
% txt1 = ['MSWD = ' num2str(MSWD)];
% txt2 = ['y = ' num2str(b) '±' num2str(b_uc) 'x + ' num2str(a) '±' num2str(a_uc)];
% txt3 = ['Age = ' num2str(t./1e6) '±' num2str(t_uc./1e6) 'Ma'];
% txt4 = ['Critical Value = ' num2str(CriVal) ' (p = 0.05)'];
% text(15,0.87,txt2,'FontSize',14);
% text(15,0.86,txt1,'FontSize',14);
% text(15,0.85,txt4,'FontSize',14);
% text(50,0.71,txt3,'FontSize',14);
% title('CP series (no K-feldspar)');
% set(gca,'FontSize',14,'linewidth',1);

%% GA
data = GA;
x = data(:,1);
x_uc = data(:,2);
y = data(:,3);
y_uc = data(:,4);
[a_GA,a_uc_GA,b,b_uc,MSWD] = YorkFit(x,y,x_uc,y_uc);
% critical value
CriVal = chi2inv(0.05,length(x)-2)./(length(x)-2);

% calculate the age and uncertainty
lambda = 1.42e-11;
lambda_uc = lambda.*2e-2;
[t, t_uc] = IsochronAge(b, b_uc, lambda, lambda_uc);

% % plot
% f4 = figure(4);
% f4.Position = [100 100 500 400];
% x_lin = linspace(max(data(:,1)).*1.1,min(data(:,1)).*0.9,100);
% y_lin = b.*x_lin + a;
% e = errorbar(data(:,1),data(:,3),data(:,4),data(:,4),data(:,2),data(:,2),".","MarkerSize",20,"LineStyle","none"); hold on;
% plot(x_lin,y_lin,'k--','LineWidth',1.5); hold on;
% xlabel('(^{87}Rb/ ^{86}Sr)');
% ylabel('(^{87}Sr/ ^{86}Sr)');
% txt1 = ['MSWD = ' num2str(MSWD)];
% txt2 = ['y = ' num2str(b) '±' num2str(b_uc) 'x + ' num2str(a) '±' num2str(a_uc)];
% txt3 = ['Age = ' num2str(t./1e6) '±' num2str(t_uc./1e6) 'Ma'];
% txt4 = ['Critical Value = ' num2str(CriVal) ' (p = 0.05)'];
% text(15,0.87,txt2,'FontSize',14);
% text(15,0.86,txt1,'FontSize',14);
% text(15,0.85,txt4,'FontSize',14);
% text(50,0.71,txt3,'FontSize',14);
% title('GA series (no K-feldspar)');
% set(gca,'FontSize',14,'linewidth',1);

%% Initial 87Sr/86Sr compare
a = [a_CP a_HD a_GA];
a_uc = [a_uc_CP a_uc_HD a_uc_GA];
w = 1./a_uc.^2;
Wmean = sum(a.*w)./sum(w);
Wmean_uc = sqrt(1./sum(w));
Weam_MSWD = sum( (a - Wmean).^2.*w )./(length(a) - 1);

f = figure(5);
f.Position = [100 100 500 400];
xaix = 1:1:3;
e = errorbar(xaix,a,a_uc,".","MarkerSize",20,"LineStyle","none"); hold on;
e.Cap.LineWidth = 1.5;
e.Bar.LineWidth = 1.5;
ylabel('(^{87}Sr/ ^{86}Sr)_{i}');
xticks(xaix);
samplename = ["CP" "HD" "GA"];
xticklabels(samplename);

xlim([0 4]);
xlimit = get(gca,'Xlim');
plot(xlimit,[Wmean Wmean],'b-','LineWidth',1.5); hold on;
plot(xlimit,[Wmean+Wmean_uc Wmean+Wmean_uc],'b--','LineWidth',1.5); hold on;
plot(xlimit,[Wmean-Wmean_uc Wmean-Wmean_uc],'b--','LineWidth',1.5); hold on;
legend('Samples','Weighted Mean');
ylim([0.7 0.73]);
set(gca,'FontSize',14,'linewidth',1);

crival = chi2inv(0.05,length(a)-1)./(length(a)-1)
%% Rb-Sr
a = 0.89123;
sig_a = 0.00012;
b = 132.34;
sig_b = b*0.01;
c = 0.70721;
lambda87 = 1.42e-11;
sig_lambda87 = lambda87.*2e-2;

pdf_a = makedist("Normal",'mu',a,'sigma',sig_a);
pdf_b = makedist("Normal",'mu',b,'sigma',sig_b);
pdf_lambda87 = makedist("Normal",'mu',lambda87,'sigma',sig_lambda87);

%% U-Pb
d = 0.015495;
sig_d = d.*0.02e-2;
lambda238 = 1.55125e-10;
sig_lambda238 = lambda238.*0.107e-2;

pdf_d = makedist("Normal",'mu',d,'sigma',sig_d);
pdf_lambda238 = makedist("Normal",'mu',lambda238,'sigma',sig_lambda238);

%% Monte Carlo 
% MC number
N = 5:1:1e4;
t_RbSr_stds = nan(size(N));
t_RbSr_aves = nan(size(N));
tt_RbSr_stds = nan(size(N));
tt_RbSr_aves = nan(size(N));
t_UPb_stds = nan(size(N));
t_UPb_aves = nan(size(N));
tt_UPb_stds = nan(size(N));
tt_UPb_aves = nan(size(N));
for I = 1:size(N,2)
   % Rb-Sr
   random_a = random(pdf_a,N(I),1);
   random_b = random(pdf_b,N(I),1);
   random_lambda87 = random(pdf_lambda87,N(I),1);
   t_RbSr = log( (random_a - c)./random_b + 1 )./lambda87;
   tt_RbSr = log( (random_a - c)./random_b + 1 )./random_lambda87;
   t_RbSr_stds(I) = std(t_RbSr);
   t_RbSr_aves(I) = mean(t_RbSr);
   tt_RbSr_stds(I) = std(tt_RbSr);
   tt_RbSr_aves(I) = mean(tt_RbSr);

   % U-Pb
   random_d = random(pdf_d,N(I),1);
   random_lambda238 = random(pdf_lambda238,N(I),1);
   t_UPb = log(random_d + 1)./lambda238;
   tt_UPb = log(random_d + 1)./random_lambda238;
   t_UPb_stds(I) = std(t_UPb);
   t_UPb_aves(I) = mean(t_UPb);
   tt_UPb_stds(I) = std(tt_UPb);
   tt_UPb_aves(I) = mean(tt_UPb);
   
   % threshold
      if I > 1000
        if abs(mean(tt_RbSr_aves(I-1000:I)) - tt_RbSr_aves(I))./ ...
           mean(tt_RbSr_aves(I-1000:I)) < 1e-6 && ...
           abs(mean(tt_UPb_aves(I-1000:I)) - tt_UPb_aves(I))./ ...
           mean(tt_UPb_aves(I-1000:I)) < 1e-6
           break
        end
      end
end
f1 = figure(1);
f1.Position = [100 100 1000 1000];
subplot(4,1,1)
plot(N,t_RbSr_aves./1e6,'LineWidth',1.5); hold on;
ylabel('age (Ma)');
title('Rb-Sr');
xlimit = get(gca,'Xlim');
plot(xlimit,[97.8552 97.8552],'r--','LineWidth',1.5);
legend('Monte Carlo Simulation','Analytical Solution');
set(gca,'FontSize',14,'linewidth',1);

subplot(4,1,2)
plot(N,t_RbSr_stds./1e6,'LineWidth',1.5); hold on;
ylabel('uncertainty (Ma)');
xlimit = get(gca,'Xlim');
plot(xlimit,[0.9799  0.9799 ],'r--','LineWidth',1.5);
set(gca,'FontSize',14,'linewidth',1);

subplot(4,1,3)
plot(N,t_UPb_aves./1e6,'LineWidth',1.5); hold on;
ylabel('age (Ma)');
title('U-Pb');
xlimit = get(gca,'Xlim');
plot(xlimit,[99.1212  99.1212],'r--','LineWidth',1.5);
set(gca,'FontSize',14,'linewidth',1);

subplot(4,1,4)
plot(N,t_UPb_stds./1e6,'LineWidth',1.5); hold on;
xlabel('Monte Carlo Number');
ylabel('uncertainty (Ma)');
xlimit = get(gca,'Xlim');
plot(xlimit,[0.0197  0.0197],'r--','LineWidth',1.5);
set(gca,'FontSize',14,'linewidth',1);

f2 = figure(2);
f2.Position = [100 100 1000 1000];
subplot(4,1,1)
plot(N,tt_RbSr_aves./1e6,'LineWidth',1.5); hold on;
ylabel('age (Ma)');
title('Rb-Sr');
xlimit = get(gca,'Xlim');
plot(xlimit,[97.8552 97.8552],'r--','LineWidth',1.5);
legend('Monte Carlo Simulation','Analytical Solution');
set(gca,'FontSize',14,'linewidth',1);

subplot(4,1,2)
plot(N,tt_RbSr_stds./1e6,'LineWidth',1.5); hold on;
ylabel('uncertainty (Ma)');
xlimit = get(gca,'Xlim');
plot(xlimit,[2.1887  2.1887],'r--','LineWidth',1.5);
set(gca,'FontSize',14,'linewidth',1);

subplot(4,1,3)
plot(N,tt_UPb_aves./1e6,'LineWidth',1.5); hold on;
ylabel('age (Ma)');
title('U-Pb');
xlimit = get(gca,'Xlim');
plot(xlimit,[99.1212  99.1212],'r--','LineWidth',1.5);
set(gca,'FontSize',14,'linewidth',1);

subplot(4,1,4)
plot(N,tt_UPb_stds./1e6,'LineWidth',1.5); hold on;
xlabel('Monte Carlo Number');
ylabel('uncertainty (Ma)');
xlimit = get(gca,'Xlim');
plot(xlimit,[0.1079 0.1079],'r--','LineWidth',1.5);
set(gca,'FontSize',14,'linewidth',1);

%% display results
disp(mean(t_RbSr_aves(I-1000:I)./1e6));
disp(mean(t_RbSr_stds(I-1000:I)./1e6));
disp(mean(t_UPb_aves(I-1000:I)./1e6));
disp(mean(t_UPb_stds(I-1000:I)./1e6));

disp(mean(tt_RbSr_aves(I-1000:I)./1e6));
disp(mean(tt_RbSr_stds(I-1000:I)./1e6));
disp(mean(tt_UPb_aves(I-1000:I)./1e6));
disp(mean(tt_UPb_stds(I-1000:I)./1e6));





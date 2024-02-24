N = 5:1:10000;
pdf1 = makedist("Normal",'mu',64.4,'sigma',1.4);
pdf2 = makedist("Normal",'mu',95.3,'sigma',3.5);
stds = nan(size(N));
aves = nan(size(N));
for I = 1:size(N,2)
    MC1 = random(pdf1,N(I),1);
    MC2 = random(pdf2,N(I),1);
    MC = MC1 + MC2;
    stds(I) = std(MC);
    aves(I) = mean(MC);
    if I > 100
    if (stds(I-1) - stds(I))./stds(I-1) < 1e-6
        break;
    end
    end
end
figure(1);
plot(N,stds,'LineWidth',1.5); hold on;
xlabel('Monte Carlo Number');
ylabel('uncertainty');
xlimit = get(gca,'Xlim');
sig = sqrt(1.4^2 + 3.5^2);
plot(xlimit,[sig sig],'r--','LineWidth',1.5);
legend('Monte Carlo Simulation','Analytical Solution');
set(gca,'FontSize',14,'linewidth',1);

figure(2);
plot(N,aves,'LineWidth',1.5); hold on;
xlabel('Monte Carlo Number');
ylabel('estimated sum value');

xlimit = get(gca,'Xlim');
plot(xlimit,[64.4+95.3 64.4+95.3],'r--','LineWidth',1.5);
legend('Monte Carlo Simulation','Analytical Solution');
set(gca,'FontSize',14,'linewidth',1);

figure(3);
diff_ave = abs(aves - (64.4+95.3));
plot(N,diff_ave,'LineWidth',1.5); hold on;
title('sum');
xlabel('Monte Carlo Number');
ylabel('Absolute difference between MC and analytical');
set(gca,'FontSize',14,'linewidth',1);


figure(4);
diff_std = abs(stds - sig);
plot(N,diff_std,'LineWidth',1.5); hold on;
title('uncertainty');
xlabel('Monte Carlo Number');
ylabel('Absolute difference between MC and analytical');
set(gca,'FontSize',14,'linewidth',1);

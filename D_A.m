load('NdSr.mat');
%% calculate CHUR
t = 82.7795e6;

Nd143_CHUR = 0.51264;
Sm147_CHUR = 0.1936;
lambda147 = 6.54e-12;
Nd143_i_CHUR = Nd143_CHUR - Sm147_CHUR.*(exp(lambda147.*t) - 1);

Rb87_CHUR = 0.089;
Sr87_CHUR = 0.7045;
lambda87 = 1.42e-11;
Sr87_i_CHUR = Sr87_CHUR - Rb87_CHUR.*(exp(lambda87.*t) - 1);

%% samples epsilon
% Sm-Nd
Sm147_sp = SmNd(1:end-1,1);
Sm147_sp_uc = SmNd(1:end-1,2);
Nd143_sp = SmNd(1:end-1,3);
Nd143_sp_uc = SmNd(1:end-1,4);

Nd143_i_sp = Nd143_sp - Sm147_sp.*(exp(lambda147.*t) - 1);
Nd143_i_uc_sp = sqrt( Nd143_sp_uc.^2 + (exp(lambda147.*t) - 1).^2.*Sm147_sp_uc.^2 );

eNd = ( Nd143_i_sp./Nd143_i_CHUR - 1).*1e4;
eNd_uc = Nd143_i_uc_sp.*1e4./Nd143_i_CHUR;
%% Rb-Sr
Rb87_sp = RbSr(1:end-1,1);
Rb87_sp_uc = RbSr(1:end-1,2);
Sr87_sp = RbSr(1:end-1,3);
Sr87_sp_uc = RbSr(1:end-1,4);

Sr87_i_sp = Sr87_sp - Rb87_sp.*(exp(lambda87.*t) - 1);
Sr87_i_uc_sp = sqrt( Sr87_sp_uc.^2 + (exp(lambda87.*t) - 1).^2.*Rb87_sp_uc.^2 );

eSr = ( Sr87_i_sp./Sr87_i_CHUR - 1).*1e4;
eSr_uc = Sr87_i_uc_sp.*1e4./Sr87_i_CHUR;

%% Fit all
[a,~,b,~,MSWD] = YorkFit(eSr,eNd,eSr_uc,eNd_uc);
CritVal = chi2inv(0.95,length(eSr)-2)./(length(eSr)-2);
f1 = figure(1);
f1.Position = [100 100 500 400];
e = errorbar(eSr,eNd,eNd_uc,eNd_uc,eSr_uc,eSr_uc,".","MarkerSize",20,"LineStyle","none"); hold on;
xlabel('\epsilon_{Sr}');
ylabel('\epsilon_{Nd}');
x = linspace(max(eSr).*1.1,min(eSr).*0.9,100);
y = b.*x + a;
plot(x,y,'k--','LineWidth',1.5); hold on;
txt = ['MSWD = ' num2str(MSWD) ];
text(25,-2.5,txt,'FontSize',14);
txt2 = ['Critical Value = ' num2str(CritVal) ' (95%)'];
text(25,-3,txt2,'FontSize',14);
%% Mantle 
disp(-a./b)
%% Fit Youlumne
[a,~,b,~,MSWD_Y] = YorkFit(eSr(1:5),eNd(1:5),eSr_uc(1:5),eNd_uc(1:5));
CritVal_Y = chi2inv(0.95,length(eSr(1:5))-2)./(length(eSr(1:5))-2);
e2 = errorbar(eSr(1:5),eNd(1:5),eNd_uc(1:5),eNd_uc(1:5),eSr_uc(1:5),eSr_uc(1:5),".","MarkerSize",20,"LineStyle","none"); hold on;
e2.Color = [0.8500 0.3250 0.0980];
xlabel('\epsilon_{Sr}');
ylabel('\epsilon_{Nd}');
x = linspace(max(eSr).*1.1,min(eSr).*0.9,100);
y = b.*x + a;
p = plot(x,y,'--','LineWidth',1.5); hold on;
p.Color = [0.8500 0.3250 0.0980];
txt = ['MSWD = ' num2str(MSWD_Y) ];
text(17,-8,txt,'FontSize',14);
txt2 = ['Critical Value = ' num2str(CritVal_Y) ' (95%)'];
text(17,-8.5,txt2,'FontSize',14);
set(gca,'FontSize',14,'linewidth',1);

%% Mantle 
disp(-a./b)

%% This is a function of GoF test to see if error follows 
% a certain distribution, pd 
function [h,p] = GoFtest(err,pd)
f1 = figure(1);
%% get chi of GoF
% number of bins, applying bin span = 0.4*std
m = round(( max(err) - min(err) )./( 0.4*pd.sigma )); 
h = histogram(err,m); hold on;
f = h.Values; % values in each bin
bin_edge = h.BinEdges; % bins edges, the left edge
bin_width = bin_edge(m)-bin_edge(m-1);

% calculate the different between histogram and assumed distribution pdf
chi = 0;
n = size(err,1);
F = zeros(size(f)); % frequency area calculated by cdf
for I = 1:m 
    F(I) = ( cdf(pd,bin_edge(I)+bin_width) - cdf(pd,bin_edge(I)) ).*n;
    chi = chi + ((f(I) - F(I)).^2)./F(I);
end

%% p value test
dof = m-1-pd.NumParameters; % degree of freedom
% the expectation (some called mean) of chi2 = dof
if chi > dof % chi2 value falls on the right
    p = 1 - cdf('chi2',chi,dof); % p value is to the right tail
else
    p = cdf('chi2',chi,dof); % p value is to the left tail
end
if p >= 0.05 % p-value is larger than the confidence 
    h = 0; % null hyphothesis cannot be rejected, 
           % i.e., error follows the presumed distribution
else
    h = 1;
end

%% Plot
figure(1);
length = bin_edge(m) - bin_edge(1);
% span = linspace(bin_edge(1)-length./5,2.*bin_edge(m)-bin_edge(m-1)+length./5,100);
span = linspace(-1.7*length./2,1.7*length./2,100);
plot(span,n*pdf(pd,span),'LineWidth',2); hold on;
ylabel('frequency');
set(gca,'FontSize',14,'linewidth',1);
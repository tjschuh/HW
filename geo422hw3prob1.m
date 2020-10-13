function X2 = geo422hw3prob1(N)
% X2 = GEO422HW3PROB1(N)
%
% INPUT:
%
% N    user-defined set of numbers of size N
%
% OUTPUT:
%
% X2   a statistic based on observed frequencies, f, in a
%      histogram of N observations, and predicted frequencies
%      F
%
% TESTED ON: 9.8.0.1451342 (R2020a) Update 5
%
% Written by tschuh@princeton.edu, 10/13/2020

% Generate an Nx1 array with random numbers drawn from a standard normal distribution
R = randn(N,1);

% Calculate 0:10:100 -th percentiles of R
percents = prctile(R,0:10:100);

% Calculate the frequencies
f = histc(R,[percents(1:end-1) percents(end) + 1]);

%f(end) = 1;
  
f = f(1:end-1);
%bar(percents,f)
% Some plotting
%bar(percents,[f ; 0])

% So now percents contains our integration limits
F = [[normcdf(percents(2:end))-normcdf(percents(1:end-1))]*N]';

for i = 1:size(f)
    X2(i) = ((f(i) - F(i))^2)/F(i);
end
X2 = sum(X2);

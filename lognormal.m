clc
clear 
close all

sigma = 9; % standard deviation
mu = 30; % mean
n = 1000;

%Dataset ----------------------------------------
%%-- in dB
x_dB = sigma*randn(n, 1) + mu;
y_dB = sort(x_dB);
%%-- in watts
x_watts = 10.^(x_dB/10);
y_watts = sort(x_watts);
%------------------------------------------------

%PDFs -------------------------------------------
%%-- Paramters
pdB.m = mean(y_dB);
pdB.s = std(y_dB);

pwatts.m0 = mean( y_watts );
pwatts.s0 = std( y_watts );

pwatts.m1 = mean( log(y_watts) );
pwatts.s1 = std( log(y_watts) );

%%-- Gaussian
gaussian_dB = ( 1/(pdB.s*sqrt(2*pi)) )*exp(-0.5.*((y_dB - pdB.m)/pdB.s).^2);
gaussian_watts = ( 1/(pwatts.s0 * sqrt(2*pi)) )*exp(-0.5.*((y_watts-pwatts.m0)/pwatts.s0).^2);

%   Lognormal
% lognormal_watts = lognpdf(y_watts, pwatts.m1, pwatts.s1);
lognormal_watts = ( 1./(pwatts.s1 * sqrt(2*pi).*y_watts )).*exp(-0.5.*((log(y_watts)-pwatts.m1)/pwatts.s1).^2);
%------------------------------------------------

%Mean and variance of a lognormal distribuition--
m = exp( pwatts.m1 + (0.5*pwatts.s1^2) );
v = exp( 2*pwatts.m1 + pwatts.s1^2)*(exp(pwatts.s1^2)- 1 );
%------------------------------------------------

%Mean and variance of a normal distribuition-----
muW = exp(2*log(m) - 0.5*(log(v + m^2)) );
mudB = 10*log10(exp(2*log(m) - 0.5*(log(v + m^2)) ))
sigmadB = 10*log10(exp(sqrt((-2*log(m)+log(v + m^2)))))
%------------------------------------------------

%Plot--------------------------------------------
subplot(2,1,1);
histogram(x_dB, 'Normalization', 'pdf', 'BinWidth',1);
hold on
plot(y_dB, gaussian_dB, '-r', 'LineWidth',3)
hold on
plot(pdB.m, 0:max(gaussian_dB)/256:max(gaussian_dB), '*r')
hold on
plot(10*log10(pwatts.m0), 0:max(gaussian_dB)/256:max(gaussian_dB), '.k')
hold on
plot(mudB, 0:max(gaussian_dB)/256:max(gaussian_dB), '.y')

subplot(2,1,2);
histogram(x_watts, 'Normalization', 'pdf', 'BinWidth',1);
hold on ;
plot(y_watts, gaussian_watts, '-r', 'LineWidth',3)
hold on ;
plot(y_watts, lognormal_watts, '-y', 'LineWidth',3)
hold on ;
plot(pwatts.m0, 0:max(lognormal_watts)/256:max(lognormal_watts), '.k')
hold on ;
plot(muW, 0:max(lognormal_watts)/256:max(lognormal_watts), '.y')


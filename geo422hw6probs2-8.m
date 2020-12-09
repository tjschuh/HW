%Problem 2
%saved Frederik's 2 coulumn YEARLY.plt file as a .txt file
%1st column is calendar year, 2nd column is measurement
yearly = importdata('YEARLY.txt');

subplot(221)
plot(yearly(:,1),yearly(:,2))
title('YEARLY data')
xlabel('year')
ylabel('mysterious measurement')
xlim([1690 2017])

subplot(222)
periodogram(yearly(:,2))
title('Periodogram PSD Estimate')
xlabel({'Normalized Frequency','(\times\pi rad/sample)'})
ylabel({'Power/frequency','(dB/(rad/sample))'})

%Problem 3
%saved Frederik's 3 coulumn DECADAL.plt file as a .txt file
%1st column is # of years before calendar year 1950,
%2nd column is measurement, 3rd column is measurement uncertainty
decadal = importdata('DECADAL.txt');

subplot(223)
plot(decadal(:,1),decadal(:,2))
title('DECADAL data')
xlabel({'# of years before','calendar year 1950'})
ylabel('mysterious measurement')
xlim([-300 11700])

subplot(224)
periodogram(decadal(:,2))
title('Periodogram PSD Estimate')
xlabel({'Normalized Frequency','(\times\pi rad/sample)'})
ylabel({'Power/frequency','(dB/(rad/sample))'})
%% 
%Problem 4
%popular data windows are:
%Bartlett, Blackman, Chebwin, Hamming, Hann, Kaiser

%exploring windowing on YEARLY data
figure
periodogram(yearly(:,2),bartlett(length(yearly(:,2))))
title({'Periodogram PSD Estimate','YEARLY - BARTLETT'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-30 60])

figure
periodogram(yearly(:,2),blackman(length(yearly(:,2))))
title({'Periodogram PSD Estimate','YEARLY - BLACKMAN'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-30 60])

figure
periodogram(yearly(:,2),chebwin(length(yearly(:,2))))
title({'Periodogram PSD Estimate','YEARLY - CHEBWIN'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-30 60])

figure
periodogram(yearly(:,2),hamming(length(yearly(:,2))))
title({'Periodogram PSD Estimate','YEARLY - HAMMING'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-30 60])

figure
periodogram(yearly(:,2),hann(length(yearly(:,2))))
title({'Periodogram PSD Estimate','YEARLY - HANN'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-30 60])

figure
periodogram(yearly(:,2),kaiser(length(yearly(:,2))))
title({'Periodogram PSD Estimate','YEARLY - KAISER'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-30 60])

%%exploring windowing on DECADAL data
figure
periodogram(decadal(:,2),bartlett(length(decadal(:,2))))
title({'Periodogram PSD Estimate','DECADAL - BARTLETT'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-60 60])

figure
periodogram(decadal(:,2),blackman(length(decadal(:,2))))
title({'Periodogram PSD Estimate','DECADAL - BLACKMAN'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-60 60])

figure
periodogram(decadal(:,2),chebwin(length(decadal(:,2))))
title({'Periodogram PSD Estimate','DECADAL - CHEBWIN'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-60 60])

figure
periodogram(decadal(:,2),hamming(length(decadal(:,2))))
title({'Periodogram PSD Estimate','DECADAL - HAMMING'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-60 60])

figure
periodogram(decadal(:,2),hann(length(decadal(:,2))))
title({'Periodogram PSD Estimate','DECADAL - HANN'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-60 60])

figure
periodogram(decadal(:,2),kaiser(length(decadal(:,2))))
title({'Periodogram PSD Estimate','DECADAL - KAISER'})
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/(rad/sample))')
ylim([-60 60])
%% 
%Problem 6
%making a PSD of YEARLY data without using periodogram

%period
T=10;

%data
data=yearly(:,2);

H=hamming(length(data));
H=H/sqrt([H'*H]);
PSD=abs(fft((data(:)-mean(data)).*H*sqrt(length(data)),length(data))).^2;
selekt=1:floor(length(data)/2)+1;
PSD=PSD(selekt);
xsint=1;
X=(selekt - 1)'/xsint/length(data);

%plotting the PSD I generated and the PSD created
%from periodogram earlier in the assignment
subplot(211)
plot(X,10*log10(PSD))
set(gca,'XTick',[0 1/T 1/2/xsint]);
grid on
set(gca,'TickDir','out','TickLength',[0.02 0.025])
title('Power Spectral Density Estimate NOT using Periodogram')
xlabel('Normalized Frequency (\times 2\pi rad/year)')
ylabel('Power/frequency')

subplot(212)
periodogram(yearly(:,2),hamming(length(yearly(:,2))))
xlabel('Normalized Frequency (\times \pi rad/year)')
%% 
%Problems 7 and 8

%specific set of frequencies between 0.05 and 0.15
%chose this set because our frequency fo interest is ~0.1
freq = [0.05:0.001:0.15];

%for loop goes through every frequency from our array and
%computes the sine and cosine coefficients among other things
for i = 1:size(freq,2)
    %G is our Fourier Series coefficient matrix
    %G = [a0 amcos(2pift) bmsin(2pift)]
    G = [ones(size(yearly,1),1) ...
        cos(2*pi*freq(i)*yearly(:,1)) ...
        sin(2*pi*freq(i)*yearly(:,1))];

    %mhat uses G to compute our coefficients m1, m2, and m3
    %mhat = (G^T*G)^-1 * G^T * d
    mhat(:,i) = ((G'*G)\G')*yearly(:,2);
    
    %G*mhat = d, where d is our estimated data
    %as opposed to the true yearly data
    estD(:,i) = G*mhat(:,i);
    
    %compute the residuals which are found
    %using the estimated data - the true data
    residual(:,i) = estD(:,i) - yearly(:,2);
    
    %compute the sums of squares of the residuals
    sumsq(i,:) = sum((residual(:,i)).^2);
    
    %compute the sums of squares of the expansion coefficients
    %for the particular frequencies under consideration
    sqexp(i,:) = sum(mhat(2,i)^2 + mhat(3,i)^2);
    
    %compute variance of the residuals
    resvar(i,1) = var(residual(:,i));
end

%compute variance of real data
truevar = var(yearly(:,2));

%compute the variance ratio that we want
varratio = resvar./truevar;

%plot the expansion coefficients as a function of frequency
figure
plot(freq,mhat(1,:),'k',freq,mhat(2,:),'r',freq,mhat(3,:),'b')
title('Expansion Coefficients as a function of Frequency')
xlabel('Normalized Frequency (\times 2\pi rad/year)')
ylabel('Expansion Coefficient Value')
legend('m1','m2','m2')

%plot the sums of squares of the residuals on top and
%plot the sums of squares of the expansion coefficients on bottom
figure
subplot(211)
plot(freq,sumsq)
title('Sums of Squares of the Residuals as a function of Normalized Frequency')
xlabel('Normalized Frequency (\times 2\pi rad/year)')
ylabel({'Sums of Squares of Residuals','(Estimated Data - True Data)'})
subplot(212)
plot(freq,sqexp)
title('Sums of Squares of the Expansion Coefficients as a function of Normalized Frequency')
xlabel('Normalized Frequency (\times 2\pi rad/year)')
ylabel('Sums of Squares of Expansion Coefficients')

%plot the variance ratio vs frequency
figure
plot(freq,varratio)
title({'Ratio of the Variance of the Residuals to the',...
    'Variance of the Orginal Signal vs Frequency'})
xlabel('Normalized Frequency (\times 2\pi rad/year)')
ylabel({'Ratio of the Variance of the Residuals',...
    'to the Variance of the Original Signal'})
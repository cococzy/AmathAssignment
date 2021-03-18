%% Part 1 GNR Guitar Frequencies
clear; close all; clc;
[y1, Fs1] = audioread('GNR.m4a');
tr_gnr = length(y1)/Fs1; % record time in seconds

Ltime = tr_gnr;
le = length((1:length(y1))/Fs1);
ti = linspace(0,Ltime,le+1); 
t = ti(1:le);
k = (1/Ltime)*[0:le/2-1 -le/2:-1];
kshift = fftshift(k);
y = y1';

val = 500; 

tau = 0:0.1:Ltime;
filt_Yft_spec = zeros(length(y), length(tau));

for j = 1:length(tau)
    filter = exp(-val*(t - tau(j)).^2);
    Yf = filter .* y;
    Yft = fft(Yf);

    [Max, Ind] = max(abs(Yft));
    [Max_Ind] = ind2sub(size(Yft), Ind);
    Max_Val = abs(k(Max_Ind));
    
    fft_tau = 0.0001;
    fft_filt = exp(-fft_tau*(k - Max_Val).^2);
    filt_Yft = fft_filt .* Yft;
   
    filt_Yft_spec(:,j) = (fftshift(abs(filt_Yft)));
end

figure(1);
pcolor(tau, kshift, abs(filt_Yft_spec));
colormap(hot)
shading interp
ylim([200 900])
colorbar
xlabel('Time (t)'), ylabel('Music notes')
yticks([277.2, 311.1, 370, 415.3, 554.4, 698.5, 740])
yticklabels({'C4#=277.2', 'D4#=311.1', 'F4#=370', 'G4#=415.3', 'C5#=554.4', 'F5=698.5', 'F5#=740'})

%% Part 2 - Floyd Bass Fequencies
clear; close all; clc;
[y2, Fs2] = audioread('Floyd.m4a');
tr_floyd = length(y2)/Fs2;

Ltime = tr_floyd;
le = length((1:length(y2)-1)/Fs2);
ti = linspace(0,Ltime,le+1); t = ti(1:le);
k = (1/Ltime)*[0:le/2-1 -le/2:-1];
kshift = fftshift(k);

y = y2(1:length(y2)-1)';

val = 100;
tau = 0:0.5:10;
filt_Yft_spec = zeros(length(y), length(tau));

for j = 1:length(tau)
    gabor = exp(-val*(t - tau(j)).^2);
    Yf = gabor .* y;
    Yft = fft(Yf);

    [Max, Ind] = max(abs(Yft));
    [Max_Ind] = ind2sub(size(Yft), Ind);
    Max_Val = abs(k(Max_Ind));
    
    fft_tau = 0.1;
    fft_filt = exp(-fft_tau*(k - Max_Val).^2);
    filt_Yft = fft_filt .* Yft;

    filt_Yft_spec(:,j) = (fftshift(abs(filt_Yft)));

end

figure(2);
pcolor(tau,kshift,log(abs(filt_Yft_spec)+1));
colormap(hot)
shading interp
ylim([0 150])
colorbar
xlabel('Time (t)'), ylabel('Music notes')
yticks([82.4, 92.5, 98, 110, 123.5])
yticklabels({'E2=82.4', 'G2b=92.5', 'G2=98', 'A2=110', 'B2=123.5'})

%% Part 3 - Floyd Guitar Frequencies
% clear; close all; clc;
[y2, Fs2] = audioread('Floyd.m4a');
tr_floyd = length(y2)/Fs2;

Ltime = tr_floyd;
le = length((1:length(y2)-1)/Fs2);
ti = linspace(0,Ltime,le+1); t = ti(1:le);
k = (1/Ltime)*[0:le/2-1 -le/2:-1];
kshift = fftshift(k);

y = y2(1:length(y2)-1)';

prefilt_Yft = fftshift(fft(y));
prefilt_Yft(abs(prefilt_Yft) >= 300) = 0;
prefilter_Y = ifft(ifftshift(prefilt_Yft));

val = 100;
tau = 30:0.5:Ltime;
filt_Yft_spec = zeros(length(y), length(tau));

for j = 1:length(tau)
    gabor = exp(-val*(t - tau(j)).^2);
    Yf = gabor .* prefilter_Y;
    Yft = fft(Yf);
    
    [Max, Ind] = max(abs(Yft));
    [Max_Ind] = ind2sub(size(Yft), Ind);
    Max_Val = abs(k(Max_Ind));
    
    fft_tau = 0.1;
    fft_filt = exp(-fft_tau*(k - Max_Val).^2);
    filt_Yft = fft_filt .* Yft;

    filt_Yft_spec(:,j) = (fftshift(abs(filt_Yft)));

end

figure(3);
pcolor(tau,kshift,filt_Yft_spec);
colormap(hot)
shading interp
ylim([200 600])
colorbar
xlabel('Time (t)'), ylabel('Music notes')
yticks([220, 247, 293.7, 370, 440.3, 493.9, 587.3])
yticklabels({'A3=220', 'B3=247', 'D4=293.7', 'F4#=370', 'A4=440.3','B4=493.9','D5=587.3'})

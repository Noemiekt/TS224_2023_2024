%% MOREL TOM KEGL NOEMIE groupe 4
clear; close all; clc;

%% Initialisation

signal_data = load("fcno03fz.mat");

[N, x] = size(signal_data.fcno03fz);

b = randn(N,1);
RSB_db = 15;
SNR = 10^(RSB_db / 10);

window = 256; 
Nfft = window / 2; 
fft = 1024; 

%% Code
% Ajout du bruit au signal de parole
puissance_signal = mean(signal_data.fcno03fz.^2);
puissance_bruit = puissance_signal/SNR;

bruit_ajuste = sqrt(puissance_bruit).*b;

signal_avec_bruit = signal_data.fcno03fz + bruit_ajuste;

% Calcul des spectro
[spectro_sans_bruit, f_sans_bruit, t_sans_bruit] = spectrogram(signal_data.fcno03fz, window, Nfft, fft);
[spectro_avec_bruit, f_avec_bruit, t_avec_bruit] = spectrogram(signal_avec_bruit, window, Nfft, fft);


%% Affichage

% Affichage du signal bruité
figure;
subplot(2,1,1);
plot(signal_data.fcno03fz);
title('Signal sans bruit');
xlabel('Échantillons');
ylabel('Amplitude');

subplot(2,1,2);
imagesc(t_sans_bruit, f_sans_bruit, 10*log10(abs(spectro_sans_bruit)));
axis xy;
title('Spectrogramme du signal sans bruit');
xlabel('Temps (s)');
ylabel('Fréquence (Hz)');

% soundsc(signal_avec_bruit);

figure;
subplot(2,1,1);
plot(signal_avec_bruit);
title('Signal avec bruité');
xlabel('Échantillons');
ylabel('Amplitude');

subplot(2, 1, 2);
imagesc(t_avec_bruit, f_avec_bruit, 10*log10(abs(spectro_avec_bruit)));
axis xy;
title('Spectrogramme du signal avec bruit');
xlabel('Temps (s)');
ylabel('Fréquence (Hz)');
colorbar;



%% MOREL TOM KEGL NOEMIE groupe 4
clear; close all; clc;

%% Initialisation

signal_data = load("fcno03fz.mat");

signal = signal_data.fcno03fz;
[N, x] = size(signal);

b = randn(N,1);
RSB_db = 15;
SNR = 10^(RSB_db / 10);

% Initialisation des recouvrements 
window = 256;
pourcentage_recouvrement = 1/2;
recouvrement =  window *pourcentage_recouvrement; % 50%

Nfft = window / 2; 
fft = 1024; 


%% Code
% Ajout du bruit au signal de parole
puissance_signal = mean(signal_data.fcno03fz.^2);
puissance_bruit = puissance_signal/SNR;

bruit_ajuste = sqrt(puissance_bruit).*b;

signal_avec_bruit = signal + bruit_ajuste;


% Calcul des trames avec recouvrement
seg = [];
for i = 1:recouvrement:(N - window + 1)
    segment = signal(i:i+window-1).*hamming(window);
    seg = [seg; segment];
end
% 
% % Reconstruction du signal
% signal_reconstruit = [];
% for l = 1:recouvrement:N + recouvrement
% 
%     m = l:(l+recouvrement-1);
%     k = (l+recouvrement*pourcentage_recouvrement):(l+(3*pourcentage_recouvrement)*recouvrement -1);
% 
%     seg_sum = (seg(m, :) + seg(k, :)).';
%     wk_sum = (m + k);
% 
%     reconstructed_segment = seg_sum ./ wk_sum;
%     signal_reconstruit = [signal_reconstruit, reconstructed_segment];
% end

signal_reconstruit = zeros(size(signal)); 
wk = zeros(size(signal));

for i = 1:recouvrement:(N - window + 1)
    segment = signal(i:i+window-1);
    signal_reconstruit(i:i+window-1) = signal_reconstruit(i:i+window-1) + segment;
    wk(i:i+window-1) = wk(i:i+window-1) + length(segment);
end

% Normaliser le signal reconstruit en divisant par le nombre de chevauchements
signal_reconstruit = signal_reconstruit ./ wk;


%% Affichage

% figure;
figure;
plot(signal);
title('signal avec recouvrement à 50%');
xlabel('Temps (s)');
ylabel('valeur des trames');

figure;
plot(signal_reconstruit');
title('signal avec recouvrement à 50%');
xlabel('Temps (s)');
ylabel('valeur des trames');

soundsc(signal_reconstruit);

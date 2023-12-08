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

% Initialisation pour H
L = floor(recouvrement/2);
M = recouvrement+1-L;

H = zeros(L, M);

seuil =  3.5e3;


%% Code
% Ajout du bruit au signal de parole
puissance_signal = mean(signal_data.fcno03fz.^2);
puissance_bruit = puissance_signal/SNR;
bruit_ajuste = sqrt(puissance_bruit).*b;
signal_avec_bruit = signal + bruit_ajuste;

% Calcul des trames avec recouvrement
seg = [];
signal_rec = zeros(size(signal_avec_bruit))';



for i = 1:recouvrement:(N - window + 1)
    current_position = i;
    segment = signal_avec_bruit(i:i+window-1).*hamming(window);
    seg = [seg; segment];
    
    for l=1:L
        H(l, :) = segment(l:M+l-1);
    end

    [U, Sigma, V] = svd(H);
    
    K = 0;
   
    % Parcourir tous les éléments de la matrice
    K = sum(Sigma(:) > seuil);
    Sigmas = Sigma(1:K,1:K);
    Us = U(1:K,1:K);
    Vs = V(1:K,1:K);

    Hs = Us*Sigmas*Vs';

    s_antidiag = 0;
    if K ~= 0
        
        [m, n] = size(Hs);
        
        % moyennes des antidiagonales
        moyennes = zeros(1, m + n - 1);
      
        for k = 1:(m + n - 1)
            antidiagonale = diag(flipud(Hs), k - m);
            moyennes(k) = mean(antidiagonale);
        end
        result = zeros(m, n);
        
        for l = 1:m
            for j = 1:n
                result(l, j) = moyennes(l + j - 1);
            end
        end
        
        premiere_colonne = result(:, 1);
        derniere_ligne = result(end, :);
        
        signal_resul = [premiere_colonne', derniere_ligne];

        signal_rec(current_position:current_position+length(signal_resul)-1) = signal_rec(current_position:current_position+length(signal_resul)-1) + signal_resul;
        current_position = current_position + recouvrement; 

        
    end
    
    
end


signal_reconstruit = zeros(size(signal)); 
wk = zeros(size(signal));

N = length(signal_rec);

for i = 1:recouvrement:(N - window + 1)
    segment = signal_rec(i:i+window-1)';
    signal_reconstruit(i:i+window-1) = signal_reconstruit(i:i+window-1) + segment;
    wk(i:i+window-1) = wk(i:i+window-1) + length(segment);
end

% Normaliser le signal reconstruit en divisant par le nombre de chevauchements
signal_reconstruit = signal_reconstruit ./ wk;



%% Affichage

% figure;
figure;
plot(signal_reconstruit);
title('signal avec recouvrement à 50%');
xlabel('Temps (s)');
ylabel('valeur des trames');

soundsc(signal_reconstruit);

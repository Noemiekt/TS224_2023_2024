%% MOREL TOM KEGL NOEMIE groupe 4
clear; close all; clc;

%% Initialisation

sigma2 = 4;
N = 1000;
b = randn(N,1)*sigma2;

windowSize = 50;
halfWindow = floor(windowSize / 2);
segmentLength = 256;
overlap = 128;
dirac = zeros(1,N);
dirac(500) = 1;

correlogramme = zeros(1,N);

%% Code

% Estimateurs et spectres de puissance
% Fonction d'autocorrélation théorique
autocorrelation_theorique = sigma2 * dirac;

% Estimateur 1 : Fonction d'autocorrélation empirique
autocorrelation_estimee1 = xcorr(b, 'biased');

% Estimateur 2 : Fonction d'autocorrélation biaisée
autocorrelation_estimee2 = xcorr(b, 'unbiased');

% Estimateur 3 : Fonction d'autocorrélation sans biais (biais retiré)
autocorrelation_estimee3 = xcorr(b, 'coeff');

% Spectre de puissance du bruit blanc
spectre_puissance = fftshift(abs(fft(b)).^2) / N;

% Densité spectrale de puissance du bruit
densite_spectrale_puissance = sigma2 * ones(1, N);

%Periodigrammes

% Lissage pour le Périodogramme de Daniell
daniell = zeros(size(spectre_puissance));
for k = 1:N
    startIdx = max(k - halfWindow, 1);
    endIdx = min(k + halfWindow, N);
    daniell(k) = mean(spectre_puissance(startIdx:endIdx));
end

% Périodogramme de Bartlett
numSegments = floor(length(b) / segmentLength);
Pxx_Bartlett = zeros(N, 1);
for k = 1:numSegments
    segment = b((k-1)*segmentLength + 1:k*segmentLength);
    windowed = hamming(segmentLength).*segment;
    fftSegment = fft(windowed, N);
    Pxx_Bartlett = Pxx_Bartlett + abs(fftSegment).^2;
end
Pxx_Bartlett = Pxx_Bartlett / (segmentLength * numSegments);

% Périodogramme de Welch
step = segmentLength - overlap;
numSegments = 1 + floor((length(b) - segmentLength) / step);
Pxx_Welch = zeros(N, 1);
for k = 0:numSegments-1
    segment = b(k*step + 1:k*step + segmentLength);
    windowed = hamming(segmentLength).*segment;
    fftSegment = fft(windowed, N);
    Pxx_Welch = Pxx_Welch + abs(fftSegment).^2;
end
Pxx_Welch = Pxx_Welch / (segmentLength * numSegments);

% Convert to 'twosided' spectrum if needed
Pxx_Bartlett = [Pxx_Bartlett(N/2+1:end); Pxx_Bartlett(1:N/2)];
Pxx_Welch = [Pxx_Welch(N/2+1:end); Pxx_Welch(1:N/2)];


% Calcul de l'auto-corrélation pour différents décalages
for tau = 0:N-1
    sum = 0;
    for t = 1:N-tau
        sum = sum + b(t) * b(t + tau);
    end
    correlogramme(tau + 1) = sum / (N - tau); 
end


%% Affichage

figure;

% Fonction d'autocorrélation théorique
subplot(3,1,1);
plot(autocorrelation_theorique);
title('Fonction d''autocorrélation théorique');
xlabel('Décalage \tau'); 
ylabel('Autocorrélation');

% Fonction d'autocorrélation estimée (biaisée)
subplot(3,1,2);
plot(autocorrelation_estimee1);
title('Fonction d''autocorrélation estimée (biaisée)');
xlabel('Décalage \tau'); 
ylabel('Autocorrélation'); 

% Fonction d'autocorrélation estimée (sans biais)
subplot(3,1,3);
plot(autocorrelation_estimee2);
title('Fonction d''autocorrélation estimée (non biaisée)');
xlabel('Décalage \tau');
ylabel('Autocorrélation');

figure;

subplot(2,1,1);
plot(spectre_puissance);
title('Spectre de puissance du bruit blanc');
ylim([-50 150]);

subplot(2,1,2);
plot(densite_spectrale_puissance);
title('Densité spectrale de puissance du bruit');

% Comparaison
figure;
subplot(3,1,1);
plot(daniell);
title('Périodogramme de Daniell');

subplot(3,1,2);
plot(Pxx_Bartlett);
title('Périodogramme de Bartlett');

subplot(3,1,3);
plot(Pxx_Welch);
title('Périodogramme de Welch');

% Tracer le corrélogramme
figure;
lags = 0:N-1; % Vecteur des décalages
plot(lags, correlogramme);
xlabel('Décalages');
ylabel('Auto-corrélation');
title('Corrélogramme');

close all;
clear all;

% Variables d'initialisation

N = 16;     % Nombre de porteuses   
Nbits = 1000;   % Nombre de bits total
Nactives = 16;  % Nombre de porteuses actives
lignes_garde = 2;   % Nombre de lignes de garde
EBN0dB = 0:0.1:6;  % Valeurs de Eb/N0 en dB
EBN0 = (EBN0dB./10).^10;    % Valeurs de Eb/N0
M = 4;           % Modulation QPSK

% Mapping QPSK

Mapping = zeros(N,Nbits);
Ligne = zeros(1,Nbits*N*2);

for i=1:N
    A = randi([0,1],1, Nbits*2);    % On crée un vecteur de 0 et de 1
    Reel = -2*A(1:2:end) + 1;        % On crée la partie réelle constituée de -1 et +1
    Imag = -2*A(2:2:end) + 1;        % On crée la partie imaginaire constituée de -1i et +1i
    Mapping(i,:) = Reel + 1i*Imag;
end

MappingLigne = reshape(Mapping,1,[]);

reelmap = sign(real(MappingLigne));
imaginairemap = sign(imag(MappingLigne));

MappingBits = zeros(1,N*Nbits*2);
MappingBits(1:2:end) = (reelmap+1)/2;
MappingBits(2:2:end) = (imaginairemap+1)/2;


% Canal AWGN

Xe = ifft(Mapping,N);

garde = Xe(N-lignes_garde+1:end,:); 
Xe_garde = [garde;Xe];  % On introduit ici un préfixe cyclique "au dessus" de Xe


ofdm_lineaire = reshape(Xe_garde, 1, []);

Px = mean(abs(ofdm_lineaire).^2);

for i=1:61      % 61 : Nombre de points de EBN0dB

    % Création du signal reçu

    sigmacarre = ((Px)/(2*log2(M)*EBN0(i)));
    random = randn(1,length(ofdm_lineaire));
    noise = sqrt(sigmacarre).*random;
    signalrecu = ofdm_lineaire + noise;

    % Démodulation

    Y_reshape = reshape(signalrecu,size(Xe_garde));      % On reshape la matrice en ligne
    Xe_sansprefixe = Y_reshape(lignes_garde+1:N+lignes_garde,:);      % On supprime le préfixe cyclique
    Y_recu = fft(Xe_sansprefixe,N);
    Yligne = reshape(Y_recu,1,[]);

    % Démapping

    reel = sign(real(Yligne));
    imaginaire = sign(imag(Yligne));

    YBits = zeros(1,N*Nbits*2);
    YBits(1:2:end) = (reel+1)/2;
    YBits(2:2:end) = (imaginaire+1)/2;

    TEB(i) = mean(MappingBits ~= YBits);
end

% Tracé du TEB

TEB_theorique = 2*qfunc(sqrt(2*log2(M)*10.^(EBN0dB/10))*sin(pi/M))/log2(M); % TEB QPSK Théorique
semilogy(EBN0dB,TEB,'b-')
hold on
semilogy(EBN0dB,TEB_theorique,'r-')


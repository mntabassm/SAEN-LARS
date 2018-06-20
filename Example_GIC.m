%% Example - explaining the usage of c-LARS-GIC method.
%% Finding the sparisty level (or model order) with a Uniform Linear Array (ULA)
%% in Compressed Beamforming Application.
%
% This example is for second simulation setup (i.e., Fig. 2) in the paper,
% when the number of sensors in the ULA is n = 40.
% Code by Muhammad Naveed Tabassum and Esa Ollila, Aalto University. <2018>

% Firstly, generate the design matrix (i.e., dictionary) X for 'p' possible
% look directions (angles) using a ULA composed 'n' sensor elements.
% The ULA has an inter-element spacing of half a wavelength.

clearvars;  close all hidden;   clc;

p = 90;    % possible look directions
n = 40;    % #of Sensors in ULA
TH_deg = 180*(0:p-1)/p - 90;
X = (1/sqrt(n)) * exp( 1i*pi*(0:n-1)'*sin(deg2rad(TH_deg)) );
sdX = sqrt( sum(X.*conj(X)) );
X = bsxfun(@rdivide, X, sdX);
load('seed_data_gic.mat', 's','K','srcLoc','srcPow','SNRdB');

%%
sigVar = sum(srcPow.^2) / K; % signal variance = average source power = (|s1|^2+|s2|^2+...+|sK|^2)/K
noiseVar = 10^(-0.1*SNRdB)*sigVar; % noise variance based on SNR-level in dB.
[~, srcIndx] = ismember(srcLoc, TH_deg);
Xs = X(:,srcIndx);
M = 1000;
rng(s);
nD = 0;	% no. of times exact sparsity order detected
for mc = 1:M
    srcTrueK = srcPow.*exp(1i*unifrnd(0,2*pi, [K 1])); 	% signal = |s|exp(j*U(0,2*pi))
    cscg_noise = sqrt(noiseVar/2)*(randn(n, 1) + 1i*randn(n, 1));	% circularly symmetric complex gaussian (CSCG) noise
    y = Xs*srcTrueK + cscg_noise;  % Noisy observations
    
    % Finding sparsity order K using c-LARS-GIC
    [lam,A,B] = clarswlasso(y,X, 0,0);
    [estK,~,~] = clarsgic(y,X, A,B,2); % e.g., for gamma = 2
    
    if estK==K
        nD = nD + 1;
        fprintf('(%d / %d) correct detections of sparsity level\n', nD,mc);
    end
end
PD = nD/M;
fprintf('Probability of Detection, PD = %.3f\n', PD);

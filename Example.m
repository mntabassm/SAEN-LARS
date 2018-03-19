%% Example
%% Direction of Arrival Estimation with a Uniform Linear Array in Compressed Beamforming Application.
%
% This example shows 3 cases, where SAEN successfully recovers the exact
% support in all three cases. However, Lasso fails to recover exact support
% in all three cases. Nevertheless, due to other \alpha than unity, elastic
% net was able to recover exact support once.
% NOTE: These 3 cases are for set-up 4 taken from 1000 cases of Table IV in the paper.
% Code by Muhammad Naveed Tabassum and Esa Ollila, Aalto University. <2018>

% Firstly, generate the design matrix (i.e., dictionary) X for 'p' possible
% look directions (angles) using a uniform linear array (ULA) composed 'n'
% sensor elements. The ULA has an inter-element spacing of half a wavelength.

clearvars;  close all hidden;   clc;

p = 180;    % possible look directions
n = 40;     % #of Sensors in ULA
TH_deg = 180*(0:p-1)/p - 90;
X = (1/sqrt(n)) * exp( 1i*pi*(0:n-1)'*sin(deg2rad(TH_deg)) );
sdX = sqrt( sum(X.*conj(X)) );
X = bsxfun(@rdivide, X, sdX);
load('seed_data.mat');

%%
sigVar = sum(srcPow.^2) / K; % signal variance = average source power = (|s1|^2+|s2|^2+...+|sK|^2)/K
noiseVar = 10^(-0.1*SNRdB)*sigVar; % noise variance based on SNR-level in dB.
[~, srcIndx] = ismember(srcLoc, TH_deg);
Xs = X(:,srcIndx);
m = 51;
AL = linspace(1,5e-1, m);  % Grid of elastic net tunning parameter with 'm' values.X = (1/sqrt(n))*exp(1i*pi*(0:(n-1))'*sin(TH));
ALGOz = {'SAEN', 'EN', 'Lasso'};
rng(s);
nERS = zeros(length(ALGOz), 1);	% no. of exact recover of support (ERS)
for mc = 1:3
    srcTrueK = srcPow.*exp(1i*unifrnd(0,2*pi, [K 1])); 	% signal = |s|exp(j*U(0,2*pi))
    cscg_noise = sqrt(noiseVar/2)*(randn(n, 1) + 1i*randn(n, 1));	% circularly symmetric complex gaussian (CSCG) noise
    y = Xs*srcTrueK + cscg_noise;  % Noisy observations
    
    % SAEN
    [betaK_saen] = saen(y,X,K,AL);
    nzi = find(abs(betaK_saen));
    if sum(ismember(srcIndx, nzi))==K && length(nzi)<=K
        nERS(1) = nERS(1) + 1;
    end
    
    % EN
    [betaK_en,~,~] = cpwwen(y,X,K,AL);
    nzi = find(abs(betaK_en));
    if sum(ismember(srcIndx, nzi))==K && length(nzi)<=K
        nERS(2) = nERS(2) + 1;
    end
    
    % Lasso
    [betaK_lasso,~,~] = cpwwen(y,X,K);
    nzi = find(abs(betaK_lasso));
    if sum(ismember(srcIndx, nzi))==K && length(nzi)<=K
        nERS(3) = nERS(3) + 1;
    end
end
table(nERS, 'rowNames', ALGOz, 'VariableNames', {'Exact_support_recovery'})

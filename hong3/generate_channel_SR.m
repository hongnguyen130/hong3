function ch = generate_channel_SR(parGenCH)
% ============================================================
% Generate normalized small-scale fading channels only
% (NO path loss included)
% ============================================================

% if isfield(parGenCH, 'seed')
%     rng(parGenCH.seed);
% end

M = parGenCH.M;

K_BR = parGenCH.K_BR;
K_RU = parGenCH.K_RU;
K_RE = parGenCH.K_RE;

% ---------- BS -> RIS : Rician ----------
hBI_los  = ones(M,1);
hBI_nlos = (randn(M,1) + 1j*randn(M,1)) / sqrt(2);
ch.hBI = sqrt(K_BR/(K_BR+1))*hBI_los + sqrt(1/(K_BR+1))*hBI_nlos;

% ---------- RIS -> U1 : Rician ----------
g1_los  = ones(M,1);
g1_nlos = (randn(M,1) + 1j*randn(M,1)) / sqrt(2);
ch.g1 = sqrt(K_RU/(K_RU+1))*g1_los + sqrt(1/(K_RU+1))*g1_nlos;

% ---------- RIS -> U2 : Rician ----------
g2_los  = ones(M,1);
g2_nlos = (randn(M,1) + 1j*randn(M,1)) / sqrt(2);
ch.g2 = sqrt(K_RU/(K_RU+1))*g2_los + sqrt(1/(K_RU+1))*g2_nlos;

% ---------- RIS -> E : Rician ----------
gE_los  = ones(M,1);
gE_nlos = (randn(M,1) + 1j*randn(M,1)) / sqrt(2);
ch.gE = sqrt(K_RE/(K_RE+1))*gE_los + sqrt(1/(K_RE+1))*gE_nlos;

% ---------- BS -> U1/U2/E : Rayleigh ----------
ch.hB1 = (randn + 1j*randn) / sqrt(2);
ch.hB2 = (randn + 1j*randn) / sqrt(2);
ch.hBE = (randn + 1j*randn) / sqrt(2);

end
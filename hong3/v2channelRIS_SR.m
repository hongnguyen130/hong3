function coef = v2channelRIS_SR(parCh)
% ============================================================
% channelRIS_SR
% Compute effective RIS-related coefficients for both
% Block A (SCA/CVX) and Block B (PGA/backtracking).
%
% ASSUMPTION
%    - hBI, g1, g2, gE, hB1, hB2, hBE are normalized small-scale fading only
%    - Path loss is included HERE, not in generate_channel_SR.m
%
% INPUT: parCh with fields
%   phi     : Mx1 RIS phase vector, |phi_m| = 1
%
%   hBI     : Mx1 BS-IRS (normalized small-scale fading)
%   g1      : Mx1 IRS-U1 small-scale channel
%   g2      : Mx1 IRS-U2 small-scale channel
%   gE      : Mx1 IRS-E  small-scale channel
%
%   hB1     : scalar BS-U1 direct channel
%   hB2     : scalar BS-U2 direct channel
%   hBE     : scalar BS-E  direct channel
%
%   P0
%   d0, d1, d2, dE, dIRS1, dIRS2, dRSE
%   m1, m2
%   NBSa, N1a, N2a, NEa   % optional, only NBSa used here
%
% OUTPUT: coef struct
%   lambda1, lambda2, lambdaE
%   A1, B1, A2, B2, C2, S1, S1bar
%   K1, K2, KE, KC2
% ============================================================

%-----------unpack------%
phi = parCh.phi(:);
M=parCh.M();
hBI = parCh.hBI(:)/sqrt(M);
g1  = parCh.g1(:);
g2  = parCh.g2(:);
gE  = parCh.gE(:);
%fprintf('g1=%d,g2=%d,gE=%d\n',g1,g2,gE);
hB1 = parCh.hB1;
hB2 = parCh.hB2;
hBE = parCh.hBE;
disp(hBI);
P0 = parCh.P0;

d0    = parCh.d0;
d1    = parCh.d1;
d2    = parCh.d2;
dE    = parCh.dE;
dIRS1 = parCh.dIRS1;
dIRS2 = parCh.dIRS2;
dRSE  = parCh.dRSE;

m1 = parCh.m1;
m2 = parCh.m2;

NBSa = parCh.NBSa;
N2a  = parCh.N2a;
M=parCh.M;

% ============================================================
% 1) Cascade channels lambda (NO path loss inside lambda)
% ============================================================
coef.lambda1 = (g1 .* hBI).' * phi/sqrt(M);
coef.lambda2 = (g2 .* hBI).' * phi/sqrt(M);
coef.lambdaE = (gE .* hBI).' * phi/sqrt(M);
disp(phi);
disp(coef);
% ============================================================
% 2) Precompute path-loss scaling constants
% ============================================================
coef.K1  = P0 * (dIRS1)^(-m2) * (d0)^(-m1);
coef.K2  = P0 * (dIRS2)^(-m2) * (d0)^(-m1);
coef.KE  = P0 * (dRSE )^(-m2) * (d0)^(-m1);
%for C2 term in SIC denominator at U2

coef.KC2 = (dIRS2)^(-m2) * NBSa;

% ============================================================
% 3) Effective coefficients
% ============================================================
% -------U1--------
coef.A1 = coef.K1 * abs(coef.lambda1)^2;
coef.B1 = P0 * (d1)^(-m1) * abs(hB1)^2;

% U2
coef.A2 = coef.K2 * abs(coef.lambda2)^2;
coef.B2 = P0 * (d2)^(-m1) * abs(hB2)^2;

% SIC noise-related term at U2
coef.C2 = coef.KC2 * abs(coef.lambda2)^2 + N2a;

% Eve
coef.S1    = coef.KE * abs(coef.lambdaE)^2 + P0 * (dE)^(-m2) * abs(hBE)^2;
coef.S1bar = P0 * (dE)^(-m2) * abs(hBE)^2;

end
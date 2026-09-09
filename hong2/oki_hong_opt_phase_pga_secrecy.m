function out = oki_hong_opt_phase_pga_secrecy(parB)
% ============================================================
% Block B: optimize RIS phase phi by PGA + backtracking
%
% Objective:
%   maximize F(phi) = SSR(phi) - etaLeak * max(REworst(phi)-epsLeak,0)^2
%
% Variables:
%   phi \in C^{M x 1}, |phi_n| = 1
%
% Inputs in parB:
%   phi0          : initial RIS phase vector (Mx1 complex)
%   anpha1        : power split for U1
%   anpha2        : power split for U2
%   anphaAN       : power split for AN
%
%   rho           : time-switching coefficient
%   P0            : BS transmit power
%
%   d0, d1, d2, dE, dIRS1, dIRS2, dRSE
%   m1, m2
%
%   N1a, N2a, NEa, NBSa
%
%   hBI           : BS-IRS small-scale channel (Mx1)
%   g1            : IRS-U1 small-scale channel (Mx1)
%   g2            : IRS-U2 small-scale channel (Mx1)
%   gE            : IRS-E small-scale channel (Mx1)
%
%   hB1           : direct BS-U1 channel
%   hB2           : direct BS-U2 channel
%   hBE           : direct BS-E channel
%
%   epsLeak       : epsilon leakage threshold
%   etaLeak       : penalty coefficient
%
%   mu0, beta, lsMax, maxIter, tol
%
% Outputs:
%   out.phi
%   out.F_hist
%   out.iters
%   out.met       : true metrics at final phi
% ============================================================

% ---------- initialization ----------
phi = parB.phi0(:);
phi = exp(1j*angle(phi));   % project to unit modulus

mu0   = parB.mu0;
beta  = parB.beta;
lsMax = parB.lsMax;
maxIt = parB.maxIter;
tol   = parB.tol;

F_hist = zeros(maxIt,1);

for it = 1:maxIt

    [F0, gphi, met0] = local_eval_F_grad(parB, phi);
    F_hist(it) = F0;

    % stopping by objective improvement
    if it > 1
        if abs(F_hist(it) - F_hist(it-1)) <= tol * max(1, abs(F_hist(it-1)))
            break;
        end
    end

    % ascent direction
    d = gphi;

    % if gradient too small -> stop
    if norm(d, 2) <= tol
        break;
    end

    % ---------- backtracking ----------
    mu = mu0;
    accept = false;

    for ls = 1:lsMax
        phi_try = phi + mu*d;
        phi_try = exp(1j*angle(phi_try));   % projection onto |phi_n|=1

        [F_try, ~, ~] = local_eval_F_grad(parB, phi_try);

        if F_try >= F0
            accept = true;
            break;
        end

        mu = beta * mu;
    end

    if ~accept
        break;
    end

    % variable-change stopping
    if norm(phi_try - phi, 2) <= tol
        phi = phi_try;
        break;
    end

    phi = phi_try;
end

% ---------- final evaluation ----------
[~, ~, met] = local_eval_F_grad(parB, phi);

out.phi   = phi;
out.F_hist = F_hist(1:it);
out.iters = it;
out.met   = met;

end


% ============================================================
% Local function: evaluate objective + gradient + metrics
% ============================================================
function [F, gphi, met] = local_eval_F_grad(parB, phi)

% ---------- unpack ----------
anpha1  = parB.anpha1;
anpha2  = parB.anpha2;
anphaAN = parB.anphaAN;

rho = parB.rho;
P0  = parB.P0;

d0    = parB.d0;
d1    = parB.d1;
d2    = parB.d2;
dE    = parB.dE;
dIRS1 = parB.dIRS1;
dIRS2 = parB.dIRS2;
dRSE  = parB.dRSE;

m1 = parB.m1;
m2 = parB.m2;

N1a  = parB.N1a;
N2a  = parB.N2a;
NEa  = parB.NEa;
NBSa = parB.NBSa;

hBI = parB.hBI(:);
g1  = parB.g1(:);
g2  = parB.g2(:);
gE  = parB.gE(:);

hB1 = parB.hB1;
hB2 = parB.hB2;
hBE = parB.hBE;

epsLeak = parB.epsLeak;
etaLeak = parB.etaLeak;

ln2 = log(2);

% ============================================================
% Cascade channels lambda1, lambda2, lambdaE
% IMPORTANT:
% lambda are "pure cascade channels", NO path loss included
% ============================================================
lambda1 = (g1 .* hBI).' * phi;
lambda2 = (g2 .* hBI).' * phi;
lambdaE = (gE .* hBI).' * phi;

X1 = abs(lambda1)^2;
X2 = abs(lambda2)^2;
XE = abs(lambdaE)^2;

% ============================================================
% Effective coefficients A1, A2, B1, B2, C2, S1, S1bar
% ============================================================
A1 = P0 * (dIRS1)^(-m1) * (d0)^(-m1) * X1;
B1 = P0 * (d1)^(-m2) * abs(hB1)^2;

A2 = P0 * (dIRS2)^(-m1) * (d0)^(-m1) * X2;
B2 = P0 * (d2)^(-m2) * abs(hB2)^2;

C2 = (dIRS2)^(-m1) * NBSa * X2 + N2a;

S1    = P0 * (dRSE)^(-m1) * (d0)^(-m1) * XE + P0 * (dE)^(-m2) * abs(hBE)^2;
S1bar = P0 * (dE)^(-m2) * abs(hBE)^2;

% ============================================================
% SINRs (according to your latest formulas)
% ============================================================

% ---- U1 ----
SINR_U1_x1_IRS    = (anpha1 * (A1 + B1)) / (anpha2 * (A1 + B1) + N1a);
SINR_U1_x1_nonIRS = (anpha1 * B1)        / (anpha2 * B1 + N1a);

% ---- U2 SIC for x1 ----
SINR_U2_x1_IRS = (anpha1 * (A2 + B2)) / (anpha2 * (A2 + B2) + C2);

% ---- U2 decode x2 ----
SINR_U2_x2_IRS = (anpha2 * B2) / N2a;

% ---- Eve decode x1 ----
SINR_E_x1_IRS    = (anpha1 * S1)    / ((anpha2 + anphaAN) * S1    + NEa);
SINR_E_x1_nonIRS = (anpha1 * S1bar) / ((anpha2 + anphaAN) * S1bar + NEa);

% ---- Eve decode x2 ----
SINR_E_x2_IRS    = (anpha2 * S1)    / ((anpha1 + anphaAN) * S1    + NEa);
SINR_E_x2_nonIRS = (anpha2 * S1bar) / ((anpha1 + anphaAN) * S1bar + NEa);

% ============================================================
% Rates
% ============================================================
RU1 = rho * log2(1 + SINR_U1_x1_IRS) + (1-rho) * log2(1 + SINR_U1_x1_nonIRS);
RU2 = rho * log2(1 + SINR_U2_x2_IRS);

RE1 = rho * log2( (S1 + NEa)    / ((anpha2 + anphaAN)*S1    + NEa) ) ...
    + (1-rho) * log2( (S1bar + NEa) / ((anpha2 + anphaAN)*S1bar + NEa) );

RE2 = rho * log2( (S1 + NEa)    / ((anpha1 + anphaAN)*S1    + NEa) ) ...
    + (1-rho) * log2( (S1bar + NEa) / ((anpha1 + anphaAN)*S1bar + NEa) );

Rsc1 = max(RU1 - RE1, 0);
Rsc2 = max(RU2 - RE2, 0);
SSR  = Rsc1 + Rsc2;

REworst = max(RE1, RE2);
leak_violation = max(REworst - epsLeak, 0);

F = SSR - etaLeak * leak_violation^2;

% ============================================================
% Gradient wrt phi*
% Core identities:
%   d|lambda|^2/dphi* = conj(v)*lambda
% where lambda = v^T phi, v = g .* hBI
% ============================================================
v1 = g1 .* hBI;
v2 = g2 .* hBI;
vE = gE .* hBI;

dX1_dphi = conj(v1) * lambda1;
dX2_dphi = conj(v2) * lambda2;
dXE_dphi = conj(vE) * lambdaE;

% ---- derivative of effective coefficients wrt phi* ----
K1  = P0 * (dIRS1)^(-m1) * (d0)^(-m1);
K2  = P0 * (dIRS2)^(-m1) * (d0)^(-m1);
KC2 = (dIRS2)^(-m1) * NBSa;
KE  = P0 * (dRSE)^(-m1)  * (d0)^(-m1);

dA1_dphi = K1  * dX1_dphi;
dA2_dphi = K2  * dX2_dphi;
dC2_dphi = KC2 * dX2_dphi;
dS1_dphi = KE  * dXE_dphi;

% ============================================================
% Chain rule
% ============================================================
dF_dphi = zeros(size(phi));

% ----------------- RU1 contribution -----------------
if (RU1 - RE1) > 0
    % IRS part depends on A1
    Num = anpha1 * (A1 + B1);
    Den = anpha2 * (A1 + B1) + N1a;
    SINR = Num / Den;

    dSINR_dA1 = (anpha1 * Den - Num * anpha2) / (Den^2);
    dRU1IRS_dA1 = rho * (1/ln2) * (1/(1 + SINR)) * dSINR_dA1;

    dF_dphi = dF_dphi + dRU1IRS_dA1 * dA1_dphi;
end

% ----------------- RU2 contribution -----------------
if (RU2 - RE2) > 0
    % RU2 = rho*log2(1 + anpha2*B2/N2a)
    % does NOT depend on phi in your current model
    % so no contribution
end

% ----------------- RE1 contribution -----------------
if (RU1 - RE1) > 0
    % RE1 = rho log2((S1+NEa)/((a2+aAN)S1+NEa))
    %      +(1-rho)log2((S1bar+NEa)/((a2+aAN)S1bar+NEa))
    % only IRS part depends on S1

    dRE1_dS1 = rho * (1/ln2) * ( 1/(S1 + NEa) ...
               - (anpha2 + anphaAN)/((anpha2 + anphaAN)*S1 + NEa) );

    dF_dphi = dF_dphi - dRE1_dS1 * dS1_dphi;
end

% ----------------- RE2 contribution -----------------
if (RU2 - RE2) > 0
    dRE2_dS1 = rho * (1/ln2) * ( 1/(S1 + NEa) ...
               - (anpha1 + anphaAN)/((anpha1 + anphaAN)*S1 + NEa) );

    dF_dphi = dF_dphi - dRE2_dS1 * dS1_dphi;
end

% ============================================================
% Penalty gradient from REworst
% ============================================================
if leak_violation > 0
    if RE1 >= RE2
        dREw_dS1 = rho * (1/ln2) * ( 1/(S1 + NEa) ...
                  - (anpha2 + anphaAN)/((anpha2 + anphaAN)*S1 + NEa) );
    else
        dREw_dS1 = rho * (1/ln2) * ( 1/(S1 + NEa) ...
                  - (anpha1 + anphaAN)/((anpha1 + anphaAN)*S1 + NEa) );
    end

    dPenalty_dphi = - 2 * etaLeak * leak_violation * dREw_dS1 * dS1_dphi;
    dF_dphi = dF_dphi + dPenalty_dphi;
end

gphi = dF_dphi;

% ============================================================
% Return metrics
% ============================================================
met.RU1 = RU1;
met.RU2 = RU2;
met.RE1 = RE1;
met.RE2 = RE2;
met.REworst = REworst;
met.Rsc1 = Rsc1;
met.Rsc2 = Rsc2;
met.SSR  = SSR;
met.A1 = A1; met.A2 = A2; met.C2 = C2; met.S1 = S1;
met.lambda1 = lambda1;
met.lambda2 = lambda2;
met.lambdaE = lambdaE;

end
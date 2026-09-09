function [F, gphi] = hong_eval_F_and_grad_phi(parB, phi)
% ============================================================
% Return:
%   F    : objective value (scalar)
%   gphi : gradient wrt phi (Mx1 complex)
% di dung huong, giu dung nguyen tac va ko nhay qua da
% NOTE:
% - Here we show a clean template:
%   1) compute lambdas (cascaded channels)
%   2) compute rates RU1, RU2, RE1, RE2
%   3) compute SSR = [RU1-RE1]^+ + [RU2-RE2]^+
%   4) F = SSR (or SSR - t*xi*Ptot)
%   5) gradient via chain rule (IRS parts only)
% ============================================================

rho = parB.rho;

% fixed anpha (from Block A)
a1 = parB.anpha1;
a2 = parB.anpha2;
a3 = parB.anphaAN;

% noises
N1 = parB.N1a;
N2 = parB.N2a;
NE = parB.NEa;

% cascaded vectors (Mx1)
% lambda = v^T phi, where v = g .* hBI (elementwise)
v1 = parB.v1;   % for U1: v1 = g1 .* hBI
v2 = parB.v2;   % for U2: v2 = g2 .* hBI
vE = parB.vE;   % for E : vE = gE .* hBI

% direct links gains already packed (scalars >=0)
B1 = parB.B1;   % P0*(d1)^(-m2)*|h1_0|^2   (or your equivalent)
B2 = parB.B2;   % P0*(d2)^(-m2)*|h2_0|^2
BE = parB.BE;   % P0*(dE)^(-m2)*|hE0|^2

% IRS coefficients (scalars >=0)
K1 = parB.K1;   % P0*(dIRS1)^(-m1)*(d0)^(-m1)
K2 = parB.K2;   % P0*(dIRS2)^(-m1)*(d0)^(-m1)
KE = parB.KE;   % P0*(dRSE)^(-m1)*(d0)^(-m1)

% -------- lambdas --------
lam1 = v1.' * phi;     % scalar
lam2 = v2.' * phi;
lamE = vE.' * phi;

X1 = abs(lam1)^2;
X2 = abs(lam2)^2;
XE = abs(lamE)^2;

% -------- build A terms (IRS received power part) --------
A1 = K1 * X1;
A2 = K2 * X2;
AE = KE * XE;

% ============================================================
% Example SINRs (you MUST match exactly your SINR model)
% Here I follow your earlier algebra style:
% U1 decoding x1:
SINR_U1_x1_IRS    = (a1*(A1+B1)) / (a2*(A1+B1) + N1);
SINR_U1_x1_nonIRS = (a1*(B1))    / (a2*(B1)    + N1);

% U2 decoding x2:
SINR_U2_x2_IRS    = (a2*(A2+B2)) / (N2);     % adjust if interference exists
SINR_U2_x2_nonIRS = (a2*(B2))    / (N2);

% E decoding x1 and x2 (template like you used)
SINR_E_x1_IRS     = (a1*(AE+BE)) / ((a2+a3)*(AE+BE) + NE);
SINR_E_x1_nonIRS  = (a1*(BE))    / ((a2+a3)*(BE)    + NE);

SINR_E_x2_IRS     = (a2*(AE+BE)) / ((a1+a3)*(AE+BE) + NE);
SINR_E_x2_nonIRS  = (a2*(BE))    / ((a1+a3)*(BE)    + NE);

% -------- rates --------
RU1 = rho*log2(1+SINR_U1_x1_IRS)    + (1-rho)*log2(1+SINR_U1_x1_nonIRS);
RU2 = rho*log2(1+SINR_U2_x2_IRS)    + (1-rho)*log2(1+SINR_U2_x2_nonIRS);

RE1 = rho*log2(1+SINR_E_x1_IRS)     + (1-rho)*log2(1+SINR_E_x1_nonIRS);
RE2 = rho*log2(1+SINR_E_x2_IRS)     + (1-rho)*log2(1+SINR_E_x2_nonIRS);

Rsc1 = max(RU1 - RE1, 0);
Rsc2 = max(RU2 - RE2, 0);

SSR = Rsc1 + Rsc2;

% Dinkelbach inner objective (Ptot fixed wrt phi, so optional)
if isfield(parB,'t') && isfield(parB,'xi') && isfield(parB,'Ptot')
    F = SSR - parB.t*parB.xi*parB.Ptot;
else
    F = SSR;
end

% ============================================================
% Gradient part (only IRS terms depend on phi)
% Key identity:
%   lam = v^T phi
%   d|lam|^2 / dphi* = conj(v) * lam
% We'll do chain rule: dF/dX1, dF/dX2, dF/dXE then map to phi
% ============================================================

% If secrecy is clipped to 0, gradient should be 0 for that part
dSSR_dX1 = 0;
dSSR_dX2 = 0;
dSSR_dXE = 0;

ln2 = log(2);

% ----- helper: d/dx log2(1+SINR(x)) = (1/ln2) * (1/(1+SINR)) * dSINR/dx -----

% 1) RU1 IRS depends on X1 via A1 = K1*X1
if (RU1 - RE1) > 0
    % SINR_U1_x1_IRS = a1*(A1+B1)/(a2*(A1+B1)+N1)
    Num = a1*(A1+B1);
    Den = a2*(A1+B1) + N1;
    sinr = Num/Den;

    % derivative wrt A1
    dsinr_dA1 = (a1*Den - Num*a2) / (Den^2);

    dRU1_dA1 = rho*(1/ln2)*(1/(1+sinr))*dsinr_dA1;
    dSSR_dX1 = dSSR_dX1 + dRU1_dA1 * K1;
end

% 2) RU2 IRS depends on X2 via A2 = K2*X2
if (RU2 - RE2) > 0
    % SINR_U2_x2_IRS example used above: a2*(A2+B2)/N2
    sinr = (a2*(A2+B2))/N2;
    dsinr_dA2 = a2/N2;

    dRU2_dA2 = rho*(1/ln2)*(1/(1+sinr))*dsinr_dA2;
    dSSR_dX2 = dSSR_dX2 + dRU2_dA2 * K2;
end

% 3) RE1 IRS depends on XE via AE = KE*XE (if Rsc1>0 contributes negative)
if (RU1 - RE1) > 0
    Num = a1*(AE+BE);
    Den = (a2+a3)*(AE+BE) + NE;
    sinr = Num/Den;

    dsinr_dAE = (a1*Den - Num*(a2+a3)) / (Den^2);
    dRE1_dAE  = rho*(1/ln2)*(1/(1+sinr))*dsinr_dAE;

    dSSR_dXE = dSSR_dXE - dRE1_dAE*KE;  % minus because SSR has -RE1
end

% 4) RE2 IRS depends on XE similarly (if Rsc2>0 contributes negative)
if (RU2 - RE2) > 0
    Num = a2*(AE+BE);
    Den = (a1+a3)*(AE+BE) + NE;
    sinr = Num/Den;

    dsinr_dAE = (a2*Den - Num*(a1+a3)) / (Den^2);
    dRE2_dAE  = rho*(1/ln2)*(1/(1+sinr))*dsinr_dAE;

    dSSR_dXE = dSSR_dXE - dRE2_dAE*KE;
end

% Map dSSR/dX to gradient wrt phi*:
% dX/dphi* = conj(v) * lam
gphi = dSSR_dX1 * (conj(v1) * lam1) ...
     + dSSR_dX2 * (conj(v2) * lam2) ...
     + dSSR_dXE * (conj(vE) * lamE);

end
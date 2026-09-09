function met = v2hong_eval_true_objective(parEval)
% ============================================================
% hong_eval_true_objective
% Evaluate TRUE metrics for the current solution
% using the common channel-coefficient function channelRIS_SR.m
%
% INPUT: parEval
%   Required fields:
%       phi
%       anpha1, anpha2, anphaAN
%       rho
%       P0
%
%       hBI, g1, g2, gE
%       hB1, hB2, hBE
%
%       d0, d1, d2, dE, dIRS1, dIRS2, dRSE
%       m1, m2
%
%       N1a, N2a, NEa, NBSa
%
%   Optional:
%       PcBS, PcIRS, M, xi
%
% OUTPUT: met struct
%   met.RU1, met.RU2
%   met.RE1, met.RE2, met.REworst
%   met.Rsc1, met.Rsc2
%   met.Rs1        : secrecy sum-rate
%   met.Ptot
%   met.SEE
%   met.A1, met.A2, met.B1, met.B2, met.C2, met.S1, met.S1bar
%   met.lambda1, met.lambda2, met.lambdaE
% ============================================================

% ---------- unpack ----------
anpha1  = parEval.anpha1;
anpha2  = parEval.anpha2;
anphaAN = parEval.anphaAN;

rho = parEval.rho;
P0  = parEval.P0;

N1a = parEval.N1a;
N2a = parEval.N2a;
NEa = parEval.NEa;

% ---------- compute channel coefficients ----------

% parCh.phi   = parEval.phi;
% 
% parCh.hBI   = parEval.hBI;
% parCh.g1    = parEval.g1;
% parCh.g2    = parEval.g2;
% parCh.gE    = parEval.gE;
% 
% parCh.hB1   = parEval.hB1;
% parCh.hB2   = parEval.hB2;
% parCh.hBE   = parEval.hBE;
% 
% parCh.P0    = parEval.P0;
% 
% parCh.d0    = parEval.d0;
% parCh.d1    = parEval.d1;
% parCh.d2    = parEval.d2;
% parCh.dE    = parEval.dE;
% parCh.dIRS1 = parEval.dIRS1;
% parCh.dIRS2 = parEval.dIRS2;
% parCh.dRSE  = parEval.dRSE;
% 
% parCh.m1    = parEval.m1;
% parCh.m2    = parEval.m2;
% 
% parCh.NBSa  = parEval.NBSa;
% parCh.N2a   = parEval.N2a;

%coef = channelRIS_SR(parCh);

A1    = parEval.A1;
B1    = parEval.B1;
A2    = parEval.A2;
B2    = parEval.B2;
C2    = parEval.C2;
S1    = parEval.S1;
S1bar = parEval.S1bar;
Pc=parEval.Pc;
% ---------- TRUE SINRs ----------
% U1 decoding x1
SINR_U1_x1_IRS    = (anpha1 * (A1 + B1)) / (anpha2 * (A1 + B1) + N1a);
SINR_U1_x1_nonIRS = (anpha1 * B1)        / (anpha2 * B1 + N1a);

% U2 SIC decoding x1
SINR_U2_x1_IRS = (anpha1 * (A2 + B2)) / (anpha2 * (A2 + B2) + C2);

% U2 decoding x2
SINR_U2_x2_IRS = anpha2*( A2+ B2) /N2a ;
SINR_U2_x2_nonIRS = (anpha2 * B2) / N2a;

% Eve decoding x1
SINR_E_x1_IRS    = (anpha1 * S1)    / ((anpha2 + anphaAN) * S1    + NEa);
SINR_E_x1_nonIRS = (anpha1 * S1bar) / ((anpha2 + anphaAN) * S1bar + NEa);

% Eve decoding x2
SINR_E_x2_IRS    = (anpha2 * S1)    / ((anpha1 + anphaAN) * S1    + NEa);
SINR_E_x2_nonIRS = (anpha2 * S1bar) / ((anpha1 + anphaAN) * S1bar + NEa);

% ---------- TRUE rates ----------
RU1 = rho * log2(1 + SINR_U1_x1_IRS) + (1-rho) * log2(1 + SINR_U1_x1_nonIRS);

% Theo công th?c b?n ?ang dùng hi?n t?i
RU2 = rho * log2(1 + SINR_U2_x2_IRS)+(1-rho) * log2(1 + SINR_U2_x2_nonIRS);

% Leakage rates at Eve
RE1 = rho * log2( (S1 + NEa)    / ((anpha2 + anphaAN)*S1    + NEa) ) ...
    + (1-rho) * log2( (S1bar + NEa) / ((anpha2 + anphaAN)*S1bar + NEa) );

RE2 = rho * log2( (S1 + NEa)    / ((anpha1 + anphaAN)*S1    + NEa) ) ...
    + (1-rho) * log2( (S1bar + NEa) / ((anpha1 + anphaAN)*S1bar + NEa) );

REworst = max(RE1, RE2);

% ---------- secrecy rates ----------
Rsc1 = max(RU1 - RE1, 0);
Rsc2 = max(RU2 - RE2, 0);

SSR = Rsc1 + Rsc2;   % secrecy sum-rate

% ---------- total power ----------
% default optional parameters


xi = 1;

% if isfield(parEval, 'PcBS')
%     PcBS = parEval.PcBS;
% end
% if isfield(parEval, 'PcIRS')
%     PcIRS = parEval.PcIRS;
% end
% if isfield(parEval, 'M')
%     M = parEval.M;
% end
% if isfield(parEval, 'xi')
%     xi = parEval.xi;
% end

% total consumed power
Ptx  = P0 * (anpha1 + anpha2 + anphaAN);
%Ptot = xi * Ptx + PcBS + rho * M * PcIRS;
Ptot=xi*Ptx+Pc;

% ---------- SEE ----------
SEE = SSR / max(Ptot, 1e-12);

% ---------- output ----------
met = struct();

met.RU1 = RU1;
met.RU2 = RU2;

met.RE1 = RE1;
met.RE2 = RE2;
met.REworst = REworst;

met.Rsc1 = Rsc1;
met.Rsc2 = Rsc2;
met.SSR  = SSR;

met.Ptx  = Ptx;
met.Ptot = Ptot;
met.SEE  = SEE;

% store SINRs for debugging
met.SINR_U1_x1_IRS    = SINR_U1_x1_IRS;
met.SINR_U1_x1_nonIRS = SINR_U1_x1_nonIRS;
met.SINR_U2_x1_IRS    = SINR_U2_x1_IRS;
met.SINR_U2_x2_IRS    = SINR_U2_x2_IRS;
met.SINR_E_x1_IRS     = SINR_E_x1_IRS;
met.SINR_E_x1_nonIRS  = SINR_E_x1_nonIRS;
met.SINR_E_x2_IRS     = SINR_E_x2_IRS;
met.SINR_E_x2_nonIRS  = SINR_E_x2_nonIRS;

% store coefficients for debugging
% met.A1 = A1; met.B1 = B1;
% met.A2 = A2; met.B2 = B2;
% met.C2 = C2;
% met.S1 = S1; met.S1bar = S1bar;
% 
% met.lambda1 = coef.lambda1;
% met.lambda2 = coef.lambda2;
% met.lambdaE = coef.lambdaE;

end
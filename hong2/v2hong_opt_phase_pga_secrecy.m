function out = v2hong_opt_phase_pga_secrecy(parB)
% ============================================================
% hong_opt_phase_pga_secrecy
% Block B: optimize RIS phase phi using PGA + backtracking
%
% Objective:
%   maximize F(phi) = SSR(phi) - etaLeak * max(REworst(phi)-epsLeak,0)^2
%
%Robust debug version: numerical 
% INPUT parB:
%   phi0
%   anpha1, anpha2, anphaAN
%   rho, P0
%   hBI, g1, g2, gE
%   hB1, hB2, hBE
%   d0, d1, d2, dE, dIRS1, dIRS2, dRSE
%   m1, m2
%   N1a, N2a, NEa, NBSa
%   PcBS, PcIRS, M, xi,Pc   (optional, for full metric compatibility)
%
%   mu0, beta, lsMax, maxIter, tol
%   epsLeak, etaLeak
%
% OUTPUT:
%   out.phi
%   out.F_hist
%   out.iters
%   out.met
% ============================================================

% ---------- initialization ----------
phi = parB.phi0;
phi = exp(1j*angle(phi));   % projection to feasible set |phi_n|=1
theta=angle(phi);
Pc=parB.Pc;
%fprintf('Pc=%f',Pc);

mu0    = parB.mu0;
beta   = parB.beta;
lsMax  = parB.lsMax;
maxIter = parB.maxIter;
tol    = parB.tol;
m1=parB.m1;
m2=parB.m2;
fd_eps=1e-5;%finite-difference step
gradTol=1e-8; 

F_hist = zeros(maxIter,1);
theta_hist=zeros(length(theta),maxIter);
for it = 1:maxIter

    %[F0, gphi, met0] = local_eval_F_grad(parB, phi);
    % ----- current objective -----
    phi = exp(1j*theta);
    [F0, met0] = local_eval_F(parB, phi);
    F_hist(it) = F0;
    theta_hist(:,it) = theta;
    
    % ----- numerical gradient wrt theta -----
    gth = zeros(length(theta),1);

    for m = 1:length(theta)
        theta_p = theta;
        theta_m = theta;

        theta_p(m) = theta_p(m) + fd_eps;
        theta_m(m) = theta_m(m) - fd_eps;

        phi_p = exp(1j*theta_p);
        phi_mi = exp(1j*theta_m);

        Fp = local_eval_F_only(parB, phi_p);
        Fm = local_eval_F_only(parB, phi_mi);

        gth(m) = (Fp - Fm)/(2*fd_eps);
    end
    
    gnorm = norm(gth);

    %fprintf('[PGA] it=%d, F=%.6f, ||grad||=%.3e\n', it, F0, gnorm);

     % ----- stop if gradient too small -----
    if gnorm <= gradTol
        %fprintf('[PGA] Stop: gradient too small.\n');
        break;
    end
    
%     % stopping by objective improvement
%     if it > 1
%         if abs(F_hist(it) - F_hist(it-1)) <= tol * max(1, abs(F_hist(it-1)))
%             break;
%         end
%     end
% 
%     % stopping by gradient norm
%     if norm(gphi, 2) <= tol
%         break;
%     end

%     % ascent direction
%     d = gphi;

    % ---------- backtracking (line search) on theta ----------
    mu = mu0;
    accept = false;

    for ls = 1:lsMax
%         phi_try = phi + mu*d;
%         phi_try = exp(1j*angle(phi_try));   % project back to unit-modulus
        theta_try = theta + mu*gth;          % ascent
        theta_try = wrapToPi_local(theta_try);

        phi_try = exp(1j*theta_try);
        F_try = local_eval_F_only(parB, phi_try);

        if F_try >= F0
            accept = true;
            break;
        end

        mu = beta * mu;
    end

    if ~accept
        %fprintf('[PGA] Line search failed. Stop.\n');
        break;
    end
%---update-----
%     if norm(phi_try - phi, 2) <= tol
%         phi = phi_try;
%         break;
%     end

   % phi = phi_try;
   if norm(theta_try - theta) <= tol
        theta = theta_try;
        %fprintf('[PGA] Stop: theta update below tol.\n');
        break;
    end

    theta = theta_try;
end

% ---------- final evaluation ----------
phi = exp(1j*theta);
[F_final, met] = local_eval_F(parB, phi);
%[~, ~, met] = local_eval_F_grad(parB, phi);

out.phi    = phi;
out.theta=theta;
out.theta_hist=theta_hist(:,1:it);
out.F_hist = F_hist(1:it);
out.iters  = it;
out.met    = met;
out.F_final=F_final;

end


% ============================================================
% Evaluate objective only
% ============================================================
function F = local_eval_F_only(parB, phi)
    [F, ~] = local_eval_F(parB, phi);
end


% ============================================================
% Evaluate objective and gradient wrt phi*
% ============================================================
function [F, met] = local_eval_F(parB, phi)

% ---------- leakage penalty params ----------
epsLeak = parB.epsLeak;
etaLeak = parB.etaLeak;

rho = parB.rho;
P0  = parB.P0;

anpha1  = parB.anpha1;
anpha2  = parB.anpha2;
anphaAN = parB.anphaAN;

N1a = parB.N1a;
N2a = parB.N2a;
NEa = parB.NEa;
NBSa = parB.NBSa;

ln2 = log(2);

% ============================================================
% TRUE evaluation using the common evaluator
% ============================================================
parEval = struct();

parEval.phi      = phi;
parEval.anpha1   = anpha1;
parEval.anpha2   = anpha2;
parEval.anphaAN  = anphaAN;

parEval.rho      = parB.rho;
parEval.P0       = parB.P0;

parEval.hBI      = parB.hBI;
parEval.g1       = parB.g1;
parEval.g2       = parB.g2;
parEval.gE       = parB.gE;

parEval.hB1      = parB.hB1;
parEval.hB2      = parB.hB2;
parEval.hBE      = parB.hBE;

parEval.d0       = parB.d0;
parEval.d1       = parB.d1;
parEval.d2       = parB.d2;
parEval.dE       = parB.dE;
parEval.dIRS1    = parB.dIRS1;
parEval.dIRS2    = parB.dIRS2;
parEval.dRSE     = parB.dRSE;
parEval.m1=parB.m1;
parEval.m2=parB.m2;
parEval.NBSa = parB.NBSa;
parEval.N2a  = parB.N2a;
% parEval.A1 = parB.A1;
% parEval.A2 = parB.A2;
% parEval.B1 = parB.B1;
% parEval.B2 = parB.B2;
% parEval.C2 = parB.C2;
% parEval.S1 = parB.S1;
% parEval.S1bar = parB.S1bar;
coef=v2channelRIS_SR(parEval);
parEval.A1 = coef.A1;
parEval.A2 = coef.A2;
parEval.B1 = coef.B1;
parEval.B2 = coef.B2;
parEval.C2 = coef.C2;
parEval.S1 = coef.S1;
parEval.S1bar = parB.S1bar;

parEval.m1       = parB.m1;
parEval.m2       = parB.m2;

parEval.N1a      = parB.N1a;
parEval.N2a      = parB.N2a;
parEval.NEa      = parB.NEa;
parEval.NBSa     = parB.NBSa;
parEval.Pc=parB.Pc;

if isfield(parB,'PcBS'),  parEval.PcBS  = parB.PcBS;  end
if isfield(parB,'PcIRS'), parEval.PcIRS = parB.PcIRS; end
if isfield(parB,'M'),     parEval.M     = parB.M;     end
if isfield(parB,'xi'),    parEval.xi    = parB.xi;    end

met = v2hong_eval_true_objective(parEval);

% ============================================================
% Objective for Block B
% ============================================================
leak_violation = max(met.REworst - epsLeak, 0);
F = met.SSR - etaLeak * leak_violation^2;
end

% ============================================================
% wrap angle to (-pi, pi]
% ============================================================
function x = wrapToPi_local(x)
    x = mod(x + pi, 2*pi) - pi;
end

% % ============================================================
% % Gradient wrt phi*
% % ============================================================
% % Use the same lambda definitions as channelRIS_SR
% v1 = parB.g1(:) .* parB.hBI(:);
% v2 = parB.g2(:) .* parB.hBI(:);
% vE = parB.gE(:) .* parB.hBI(:);
% 
% lambda1 = v1.' * phi;
% lambda2 = v2.' * phi;
% lambdaE = vE.' * phi;
% 
% dX1_dphi = conj(v1) * lambda1;   % d|lambda1|^2 / dphi*
% dX2_dphi = conj(v2) * lambda2;   % d|lambda2|^2 / dphi*
% dXE_dphi = conj(vE) * lambdaE;   % d|lambdaE|^2 / dphi*
% 
% % Effective path-loss scaling
% K1  = P0 * (parB.dIRS1)^(-parB.m1) * (parB.d0)^(-parB.m1);
% K2  = P0 * (parB.dIRS2)^(-parB.m1) * (parB.d0)^(-parB.m1);
% KE  = P0 * (parB.dRSE )^(-parB.m1) * (parB.d0)^(-parB.m1);
% KC2 = (parB.dIRS2)^(-parB.m1) * NBSa;
% 
% dA1_dphi = K1  * dX1_dphi;
% dA2_dphi = K2  * dX2_dphi;
% dS1_dphi = KE  * dXE_dphi;
% dC2_dphi = KC2 * dX2_dphi;
% 
% % Read current coefficients from met
% 
% 
% % ============================================================
% % Derivatives of Rs wrt A1, A2, C2, S1
% % ============================================================
% dRs_dA1 = 0;
% dRs_dA2 = 0;
% dRs_dC2 = 0;
% dRs_dS1 = 0;
% 
% % ---------- Rsc1 = max(RU1-RE1,0) ----------
% if (met.RU1 - met.RE1) > 0
% 
%     % RU1 = rho*log2(1 + a1(A1+B1)/(a2(A1+B1)+N1))
%     Num = anpha1*(A1 + B1);
%     Den = anpha2*(A1 + B1) + N1a;
%     SINR = Num / Den;
% 
%     dSINR_dA1 = (anpha1*Den - Num*anpha2) / (Den^2);
%     dRU1_dA1 = rho * (1/ln2) * (1/(1 + SINR)) * dSINR_dA1;
% 
%     % RE1 = rho*log2((S1+NE)/((a2+aAN)S1+NE)) + ...
%     dRE1_dS1 = rho * (1/ln2) * ...
%         ( 1/(S1 + NEa) - (anpha2 + anphaAN)/((anpha2 + anphaAN)*S1 + NEa) );
% 
%     dRs_dA1 = dRs_dA1 + dRU1_dA1;
%     dRs_dS1 = dRs_dS1 - dRE1_dS1;
% end
% 
% % ---------- Rsc2 = max(RU2-RE2,0) ----------
% if (met.RU2 - met.RE2) > 0
% 
%     % RU2 = rho*log2(1 + anpha2*B2/N2a)
%     % => with current formula, RU2 does NOT depend on phi
%     % so dRU2_dA2 = 0, dRU2_dC2 = 0
% 
%     % RE2 = rho*log2((S1+NE)/((a1+aAN)S1+NE)) + ...
%     dRE2_dS1 = rho * (1/ln2) * ...
%         ( 1/(S1 + NEa) - (anpha1 + anphaAN)/((anpha1 + anphaAN)*S1 + NEa) );
% 
%     dRs_dS1 = dRs_dS1 - dRE2_dS1;
% end
% 
% % map to phi
% gphi = dRs_dA1 * dA1_dphi ...
%      + dRs_dA2 * dA2_dphi ...
%      + dRs_dC2 * dC2_dphi ...
%      + dRs_dS1 * dS1_dphi;
% 
% % ============================================================
% % Add leakage penalty gradient
% % ============================================================
% if leak_violation > 0
% 
%     if met.RE1 >= met.RE2
%         dREworst_dS1 = rho * (1/ln2) * ...
%             ( 1/(S1 + NEa) - (anpha2 + anphaAN)/((anpha2 + anphaAN)*S1 + NEa) );
%     else
%         dREworst_dS1 = rho * (1/ln2) * ...
%             ( 1/(S1 + NEa) - (anpha1 + anphaAN)/((anpha1 + anphaAN)*S1 + NEa) );
%     end
% 
%     gphi = gphi - 2*etaLeak*leak_violation * dREworst_dS1 * dS1_dphi;
% end
% 
% end
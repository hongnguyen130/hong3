function met = hong_eval_true_objective(par)
% Evaluate true secrecy + SEE for given (phi, P1,P2,PAN)

% unpack
phi  = par.phi;
P1   = par.P1;  P2 = par.P2;  PAN = par.PAN;

% effective channels: h_i = h_Bi + g_i.' * phi
h1 = par.hB1 + par.g11.'*phi;
h2 = par.hB2 + par.g22.'*phi;
hE = par.hBE + par.gEE.'*phi;

a1 = abs(h1)^2; a2 = abs(h2)^2; aE = abs(hE)^2;

sig1 = par.sig1; sig2 = par.sig2; sigE = par.sigE;
alpha = par.alpha;

% ---------- Legit rates (NOMA) ----------
R1 = log2( (a1*(P1+P2)+sig1) / (a1*P2 + sig1) );     % U1
R2 = log2( 1 + (a2*P2)/sig2 );                       % U2 after SIC

% SIC check (optional)
gamma1   = (a1*P1)/(a1*P2 + sig1);
gamma2to1= (a2*P1)/(a2*P2 + sig2);
met.SIC_ok = (gamma2to1 >= gamma1);

% ---------- Eaves rates (for U1 only in Step 1) ----------
cE  = sigE + alpha*PAN;
RE1 = log2( (aE*(P1+P2)+cE) / (aE*P2 + cE) );

% secrecy (Step 1: U1 only)
Rs1 = max(R1 - RE1, 0);

% total power
Ptx = P1 + P2 + PAN;
Ptot = par.xi*Ptx + par.Pc + par.Mris*par.Pris;

% metrics
met.R1 = R1; met.RE1 = RE1; met.Rs1 = Rs1;
met.Ptx = Ptx; met.Ptot = Ptot;

% SEE (Step 1: based on Rs1)
met.SEE = Rs1 / max(Ptot, 1e-12);

end

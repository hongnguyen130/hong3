function out = hong_opt_phase_pga_secrecy_u1(par)
%OPT_PHASE_PGA_SECRECY_U1  Block-B: optimize RIS phases (unit-modulus) via projected gradient ascent
%
%  This block assumes powers are FIXED: P1, P2, PAN are given.
%  It updates phi to maximize secrecy of U1:
%       Rs1(phi) = [R1(phi) - RE1(phi)]^+
%  where
%       R1  = log2( (|h1|^2 (P1+P2)+sig1) / (|h1|^2 P2 + sig1) )
%       RE1 = log2( (|hE|^2 (P1+P2)+cE) / (|hE|^2 P2 + cE) )
%       cE = sigE + alpha*PAN
%
%  Effective channels (RIS-reflection, SISO):
%       h1(phi) = hB1 + g1' * phi
%       hE(phi) = hBE + gE' * phi
%  IMPORTANT: here we use g' * phi (no conjugate) consistently.
%  If your model uses g^H phi, then set g = conj(g) before passing in.

% -------------------- Unpack --------------------
phi  = par.phi0;          % (M x 1) complex, |phi_n|=1 initial
hB1  = par.hB1;           % scalar complex (direct BS->U1)
hBE  = par.hBE;           % scalar complex (direct BS->E)
g1   = par.g1;            % (M x 1) complex (cascaded RIS path to U1) such that g1' * phi
gE   = par.gE;            % (M x 1) complex (cascaded RIS path to E)

xi = par.xi;
anpha1   = par.anpha1;
anpha2   = par.P2;
PAN  = par.PAN;
Ptx = anpha1 + anpha2 + anphaAN;
Ptot = par.xi*Ptx + par.Pc + par.Mris*par.Pris;

sig1 = par.sig1;          % effective noise at U1 (include jammer if you model it)
sigE = par.sigE;          % noise at E
alpha = par.alpha;        % AN-to-E coefficient

maxIter = par.maxIter;
tol     = par.tol;

% step size + backtracking
mu0      = par.mu0;       % initial step size
beta     = par.beta;      % step shrink factor, e.g., 0.5
lsMax    = par.lsMax;     % max line-search steps

% safety
tiny = 1e-12;

M = length(phi);
obj_hist = zeros(maxIter,1);

% -------------------- Helper: objective Rs1 --------------------
    function [SEE1,Rs] = SEE_u1(phi_in)
        h1 = hB1 + (g1.' * phi_in);          % scalar
        hE = hBE + (gE.' * phi_in);          % scalar
        x  = abs(h1)^2;
        y  = abs(hE)^2;

        R1  = log2( (x*(anpha1+anpha2) + sig1) / (x*anpha2 + sig1) );
        cE  = sigE + alpha*PAN;
        RE1 = log2( (y*(anpha1+anpha2) + cE) / (y*anpha2 + cE) );

        Rs = max(R1 - RE1, 0);
        SEE1=Rs/max(Ptot,1e-12);
    end

% -------------------- Helper: Wirtinger gradient wrt phi* --------------------
    function grad = grad_secrecy_u1(phi_in)
        % Compute h and |h|^2
        h1 = hB1 + (g1.' * phi_in);
        hE = hBE + (gE.' * phi_in);
        x  = abs(h1)^2;
        y  = abs(hE)^2;

        % dR1/dx
        A1 = x*(anpha1+anpha2) + sig1;
        B1 = x*anpha2      + sig1;
        dR1_dx = (1/log(2)) * ( (anpha1+anpha2)/max(A1,tiny) - (anpha2)/max(B1,tiny) );

        % dRE1/dy
        cE = sigE + alpha*PAN;
        AE = y*(anpha1+anpha2) + cE;
        BE = y*anpha2      + cE;
        dRE1_dy = (1/log(2)) * ( (anpha1+anpha2)/max(AE,tiny) - (anpha2)/max(BE,tiny) );

        % dx/dphi* and dy/dphi* (Wirtinger):
        % if h = hB + g'phi, then ∂|h|^2/∂phi* = g * h
        dx_dphic = g1 * h1;    % (N x 1)
        dy_dphic = gE * hE;    % (N x 1)

        % secrecy gradient wrt phi* (ignore [·]^+ clipping for gradient; stable in practice)
        grad = dR1_dx * dx_dphic - dRE1_dy * dy_dphic;
    end

% -------------------- Main loop: projected gradient ascent --------------------
[SEE_old,Rs_old] = SEE_u1(phi);

for k = 1:maxIter
    obj_hist(k) = SEE_old;

    grad = grad_secrecy_u1(phi);    % gradient wrt phi*

    % ascent direction
    mu = mu0;
    improved = false;

    for ls = 1:lsMax
        phi_try = phi + mu * grad;

        % project to unit-modulus
        phi_try = exp(1j*angle(phi_try));

       [SEE_try,Rs_try] = SEE_u1(phi_try);

        if SEE_try >= SEE_old + 1e-8   % sufficient increase (very mild)
            improved = true;
            break;
        end
        mu = mu * beta;
    end

    if ~improved
        % no improvement found -> stop
        obj_hist = obj_hist(1:k);
        break;
    end

    % update
    phi_new = phi_try;
    SEE_new  = SEE_try;
    Rs_new=Rs_old;

    % convergence check
    if norm(phi_new - phi, 2) <= tol
        phi = phi_new;
        SEE_old = SEE_new;
        Rs_old = Rs_new;
        obj_hist = obj_hist(1:k);
        break;
    end

    phi = phi_new;
    SEE_old = SEE_new;
    Rs_old = Rs_new;
end

% outputs
out.phi = phi;
[SEE1,Rs1]=SEE_u1(phi);
out.SEE1 = SEE1;
out.Rs1=Rs1;
out.obj_hist = obj_hist;

end

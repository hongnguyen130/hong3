function out = hong_opt_phase_pga_secrecy(parB)
% ============================================================
% Block B: Projected Gradient Ascent (phase-only RIS)
% maximize F(phi) with |phi_n| = 1
% using backtracking line-search
% ============================================================

phi   = parB.phi0;              % Mx1 complex, |phi|=1
mu0   = parB.mu0;               % initial step
beta  = parB.beta;              % shrink factor (0.5 typical)
lsMax = parB.lsMax;             % max line-search steps
tol   = parB.tol;               % tolPGA
maxIt = parB.maxIter;

% Armijo parameters (you can tune)
c1 = 1e-4;

% Ensure feasible init
phi = exp(1j*angle(phi));

F_hist = zeros(maxIt,1);

for it = 1:maxIt

    % ---- evaluate objective + gradient at current phi ----
    [F0, g] = hong_eval_F_and_grad_phi(parB, phi);  % g: gradient wrt phi (complex)

    F_hist(it) = F0;

    % ---- stopping criterion ----
    if it > 1
        if abs(F_hist(it) - F_hist(it-1)) <= tol*max(1,abs(F_hist(it-1)))
            break;
        end
    end

    % ---- ascent direction ----
    d = g;  % gradient ascent

    % ---- backtracking line search ----
    mu = mu0;
    accept = false;

    for ls = 1:lsMax
        phi_try = phi + mu*d;
        phi_try = exp(1j*angle(phi_try)); % projection onto |phi_n|=1

        F_try = hong_eval_F_and_grad_phi(parB, phi_try); % allow 1-output version
        % Armijo (ascent): F_try >= F0 + c1*mu*Re(g^H d)
        rhs = F0 + c1*mu*real(g' * d);

        if F_try >= rhs
            accept = true;
            break;
        end
        mu = beta*mu;
    end

    if ~accept
        % If cannot find improvement, stop early
        break;
    end

    % ---- update ----
    phi = phi_try;

end

out.phi = phi;
out.F_hist = F_hist(1:it);
out.iters = it;

end
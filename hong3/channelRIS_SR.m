function ch = channelRIS_SR(par)
%GEN_CHANNELS_SISO_RIS Generate SISO channels with RIS reflection + path loss
%
% Outputs:
%   ch.h1_eff, ch.h2_eff, ch.hE_eff : effective complex channels (scalar)
%   ch.h_BU1, ch.h_BU2, ch.h_BE     : direct links (scalar)
%   ch.h_BI                         : BS->RIS link (N x 1)
%   ch.h_IU1, ch.h_IU2, ch.h_IE     : RIS->node links (N x 1)
%   ch.a1, ch.a2, ch.aE             : effective power gains |h_eff|^2
%
% Model:
%   h_i_eff = h_Bi + (h_Ii.' .* h_BI.').' * phi
%           = h_Bi + ( (h_Ii .* h_BI).' * phi )
% where phi is N x 1, |phi_n|=1 (reflection coefficients).
%
% Path loss:
%   PL(d) = C0 * d^{-alpha}   (linear scale)
% Small-scale:
%   Rayleigh: CN(0,1)
%   Rician  : sqrt(K/(K+1))*LoS + sqrt(1/(K+1))*NLoS

% ------------------ Unpack ------------------
M      = par.M;              % number of RIS elements
phi    = par.phi;            % N x 1, |phi|=1

% distances (meters)
d_BI   = par.d_BI;           % BS -> RIS
d_BU1  = par.d_BU1;          % BS -> U1
d_BU2  = par.d_BU2;          % BS -> U2
d_BE   = par.d_BE;           % BS -> E
d_IU1  = par.d_IU1;          % RIS -> U1
d_IU2  = par.d_IU2;          % RIS -> U2
d_IE   = par.d_IE;           % RIS -> E

% path loss exponents
alpha_BI  = par.alpha_BI;
alpha_BU  = par.alpha_BU;    % used for BS->U1/U2
alpha_BE  = par.alpha_BE;
alpha_IU  = par.alpha_IU;    % used for RIS->U1/U2/E

% reference gain (linear). Example: C0 = 10^(-PL0_dB/10).
C0_BI = par.C0_BI;
C0_BU = par.C0_BU;
C0_BE = par.C0_BE;
C0_IU = par.C0_IU;

% fading type: 'rayleigh' or 'rician'
%fading = par.fading;         % string
K_BI   = getfield_def(par,'K_BI',0);   %ok<GFLD> (only used if rician)
K_IU   = getfield_def(par,'K_IU',0);
K_BU   = getfield_def(par,'K_BU',0);
K_BE   = getfield_def(par,'K_BE',0);

% ------------------ Path loss (linear) ------------------
PL_BI  = C0_BI * (d_BI ^ (-alpha_BI));
PL_BU1 = C0_BU * (d_BU1^ (-alpha_BU));
PL_BU2 = C0_BU * (d_BU2^ (-alpha_BU));
PL_BE  = C0_BE * (d_BE ^ (-alpha_BE));
PL_IU1 = C0_IU * (d_IU1^ (-alpha_IU));
PL_IU2 = C0_IU * (d_IU2^ (-alpha_IU));
PL_IE  = C0_IU * (d_IE ^ (-alpha_IU));

% ------------------ Small-scale fading ------------------
% Helper to generate CN(0,1)
cn = @(m,n) (randn(m,n)+1j*randn(m,n))/sqrt(2);

%switch lower(fading)
 %   case 'rayleigh'
%   h_BI  = sqrt(PL_BI)  * cn(M,1);
%   h_IU1 = sqrt(PL_IU1) * cn(M,1);
%   h_IU2 = sqrt(PL_IU2) * cn(M,1);
%   h_IE  = sqrt(PL_IE)  * cn(M,1);
h_BU1 = sqrt(PL_BU1) * cn(1,1);
h_BU2 = sqrt(PL_BU2) * cn(1,1);
h_BE  = sqrt(PL_BE)  * cn(1,1);
%case 'rician'
% Simple LoS model: deterministic all-ones (you can replace with steering vectors)
LoS_vec = ones(M,1);

h_BI  = sqrt(PL_BI)  * ( sqrt(K_BI/(K_BI+1))*LoS_vec + sqrt(1/(K_BI+1))*cn(M,1) );
h_IU1 = sqrt(PL_IU1) * ( sqrt(K_IU/(K_IU+1))*LoS_vec + sqrt(1/(K_IU+1))*cn(M,1) );
h_IU2 = sqrt(PL_IU2) * ( sqrt(K_IU/(K_IU+1))*LoS_vec + sqrt(1/(K_IU+1))*cn(M,1) );
h_IE  = sqrt(PL_IE)  * ( sqrt(K_IU/(K_IU+1))*LoS_vec + sqrt(1/(K_IU+1))*cn(M,1) );

% h_BU1 = sqrt(PL_BU1) * ( sqrt(K_BU/(K_BU+1))*1 + sqrt(1/(K_BU+1))*cn(1,1) );
% h_BU2 = sqrt(PL_BU2) * ( sqrt(K_BU/(K_BU+1))*1 + sqrt(1/(K_BU+1))*cn(1,1) );
% h_BE  = sqrt(PL_BE)  * ( sqrt(K_BE/(K_BE+1))*1 + sqrt(1/(K_BE+1))*cn(1,1) );

%     otherwise
%         error('Unknown fading type. Use ''rayleigh'' or ''rician''.');
% end

% ------------------ Effective channels with RIS ------------------
% Cascaded term: sum_n (h_Ii(n) * h_BI(n) * phi(n))
% Note: (h_Ii .* h_BI).' * phi returns scalar
h1_eff = h_BU1 + ( (h_IU1 .* h_BI).' * phi );
h2_eff = h_BU2 + ( (h_IU2 .* h_BI).' * phi );
hE_eff = h_BE  + ( (h_IE  .* h_BI).' * phi );

% ------------------ Pack outputs ------------------
ch.h_BI  = h_BI;
ch.g1 = h_IU1;
ch.g2 = h_IU2;
ch.gE  = h_IE;

ch.hB1 = h_BU1;
ch.hB2 = h_BU2;
ch.hBE  = h_BE;

ch.h1_eff = h1_eff;
ch.h2_eff = h2_eff;
ch.hE_eff = hE_eff;

ch.a1 = abs(h1_eff)^2;
ch.a2 = abs(h2_eff)^2;
ch.aE = abs(hE_eff)^2;

end

% --------- small helper to safely get optional fields ----------
function v = getfield_def(s, name, def)
if isfield(s, name)
    v = s.(name);
else
    v = def;
end
end

function out = hong_opt_power_sca_secrecy(parA)

phi  = parA.phi;   
rho  = parA.rho;   
P0  = parA.P0;   
sig1 = parA.sig1;
sig2 = parA.sig2;
N1a=sig1 ; % effective noise at U1 (including jammer if you model it)
N2a=sig2; % effective noise at U2
sigE = parA.sigE; % noise at E
NEa=sigE;
alpha = parA.alpha; % AN-to-E coefficient

A1=parA.A1;
B1=parA.B1;
A2=parA.A2;
B2=parA.B2;
C2=parA.C2;
S1=parA.S1;
S1bar=parA.S1bar;
gamma0=parA.gamma0;
gamma1=parA.gamma1;
gamma2=parA.gamma2;



Pc=parA.Pc;
t  = parA.t;     % Dinkelbach parameter
xi = parA.xi;    % PA inefficiency coefficient
epsLeak=parA.epsLeak;
rho=parA.rho;
maxIter = 50;
tol     = parA.tol;
Rmin=1e-5;

% initialization
anpha1  = parA.init.anpha1;
anpha2  = parA.init.anpha2;
anphaAN = parA.init.anphaAN;

% Safety small number to avoid log(0)
tiny = 1e-12;

% ---initial point for SCA-----
    anpha1_k = anpha1;
    anpha2_k = anpha2;
    anphaAN_k = anphaAN;
for k = 1:maxIter 
   %2) Linearization pints
    %---RU1------%
    V_U1_IRS_k = max(anpha2_k*(A1+B1) + N1a,tiny);
    V_U1_non_k = max(anpha2_k*B1 + N1a,tiny);
    
%     V_U1_IRS_k=max(V_U1_IRS_k,tiny);
%     V_U1_non_k=max(V_U1_non_k, tiny);
    
    %---RU2------%
    %V_U2_IRS_k = anpha2_k*(A2+B2) + C2;
    %V_U2_IRS_k=max(V_U2_IRS_k,tiny);
    %----RE1----%
    
    V_EU1_IRS_k = max( (anpha2_k + anphaAN_k)*S1    + NEa, tiny );
    V_EU1_non_k = max( (anpha2_k + anphaAN_k)*S1bar + NEa, tiny );
    %---RE2----
    V_EU2_IRS_k = max( (anpha1_k + anphaAN_k)*S1    + NEa, tiny );
    V_EU2_non_k = max( (anpha1_k + anphaAN_k)*S1bar + NEa, tiny );
    
    cvx_clear
    % ------3)------------------ Solve convex surrogate --------------------
    cvx_begin quiet
        variables anpha1v anpha2v anphaANv s1 s2 tau        
        %-----RU1 surrogate---------%
        V_U1_IRS=anpha2v*(A1+B1)+N1a;
        V_U1_non=anpha2v*B1+N1a;
        rU1_IRS_surr = log((anpha1v+anpha2v)*(A1+B1) + N1a)/log(2) ...
             + ( -log(V_U1_IRS_k)/log(2) ...
                 - (1/log(2))*(1/V_U1_IRS_k)*(V_U1_IRS - V_U1_IRS_k) );

        rU1_non_surr = log((anpha1v+anpha2v)*B1 + N1a)/log(2) ...
             + ( -log(V_U1_non_k)/log(2) ...
                 + (1/log(2))*(1/V_U1_non_k)*(V_U1_non - V_U1_non_k) );
             
             
        RU1_surr = rho*rU1_IRS_surr + (1-rho)*rU1_non_surr;  
        %Tinh surrogate RE1-------%
        V_EU1_IRS = (anpha2v + anphaANv)*S1    + NEa;     % affine
        V_EU1_non = (anpha2v + anphaANv)*S1bar + NEa;     % affine
        
        rEU1_IRS_surr = log(S1 + NEa)/log(2) ...
              + ( -log(V_EU1_IRS_k)/log(2) ...
                  - (1/log(2))*(1/V_EU1_IRS_k)*(V_EU1_IRS - V_EU1_IRS_k) );

        rEU1_non_surr = log(S1bar + NEa)/log(2) ...
              + ( -log(V_EU1_non_k)/log(2) ...
                  - (1/log(2))*(1/V_EU1_non_k)*(V_EU1_non - V_EU1_non_k) );
        
        RE1_surr=rho*rEU1_IRS_surr + (1-rho)*rEU1_non_surr;
        %Tinh surrogate cua RU2------%
        %V2_IRS = anpha2v*(A2+B2) + C2;
        %V2_IRS=max(V2_IRS,tiny);
        
        rU2_IRS_surr = log(anpha2v*(A2+B2) + N2a)/log(2)- ( log(N2a)/log(2));
        rU2_non_surr = log(anpha2v*B2 + N2a)/log(2) - log(N2a)/log(2);
        
        RU2_surr=rho*rU2_IRS_surr + (1-rho)*rU2_non_surr;
        
        % Tinh surrogate cua RE2----%
        V_EU2_IRS = (anpha1v + anphaAN)*S1    + NEa;
        V_EU2_non = (anpha1v + anphaAN)*S1bar + NEa;
        
        rEU2_IRS_surr = log(S1 + NEa)/log(2) ...
              + ( -log(V_EU2_IRS_k)/log(2) ...
                  - (1/log(2))*(1/V_EU2_IRS_k)*(V_EU2_IRS - V_EU2_IRS_k) );

        rEU2_non_surr = log(S1bar + NEa)/log(2) ...
              + ( -log(V_EU2_non_k)/log(2) ...
                  - (1/log(2))*(1/V_EU2_non_k)*(V_EU2_non - V_EU2_non_k) );
        
        RE2_surr=rho*rEU2_IRS_surr + (1-rho)*rEU2_non_surr;
        
        %%%%%%%%%%%%%%%%
        Rsc1_surr=RU1_surr-RE1_surr;
        Rsc2_surr=RU2_surr-RE2_surr;
        SSR_surr=s1+s2;
        
        P1v=anpha1v*P0;
        P2v=anpha2v*P0;
        PANv=anphaANv*P0;
        
        %-------Dinkelbach inner objective-------
         
        Ptot=P1v + P2v+ PANv+Pc;
        
        maximize( SSR_surr - t*xi*Ptot);
        %neu muong dùng scalar hoa leakage -ko dung eps, có the
        %maximize( SSR_surr - t*xi*Ptot-wLeak*tau);
            
        subject to
        %----------------------power split constraints-------------
            anpha1v <= 0.95; 
            anpha1v>=0.4;
            anpha2v >= 0; 
            anpha2v<=0.1;
            anpha2v>=0.05;
            anphaANv >= 0;
            anpha1v + anpha2v + anphaANv <= 1;
            anpha1v >= anpha2v;  % power ordering for NOMA (weak user gets more power)
            %P1v>=0.5*epsP;
            anphaANv <= 0.5;
            %------------secrecy slack- ham max[ơ]+------
            s1>=0;
            s2>=0;
            s1<=Rsc1_surr;
            s2<=Rsc2_surr;
            %-----------Worst-case leakage epigraph-----%
            tau>=RE1_surr;
            tau>=RE2_surr;
            tau<=epsLeak;%bat neu chay epsilon-constraint.
            %-------SINR QoS/SIC constraints------%
            %(1)U1 decoding x1
            (A1+B1)*(anpha1v - gamma1*anpha2v) >= gamma1*N1a;
            %(2)U2 decoding x1 for SIC
            (A2+B2)*(anpha1v - gamma0*anpha2v) >= gamma0*C2;
            % C2*(anpha1 - gamma0*anpha2) >= gamma0*(N2a + anphaAN*D2);
            % U2 decoding x2
            (A2+B2)*anpha2v >= gamma2*N2a;
         

    cvx_end
%     fprintf('P1=%.3f P2=%.3f PAN=%.3f Rs_surr=%.4f\n',...
%         P1v, P2v, PANv, Rs1_surr);
% 
%     fprintf('k=%g,. P1v=%g,P2v=%g,PANv=%g\n',k,P1v,P2v,PANv);
    % Update
    cvx_status_last=cvx_status;
    
    anpha1_new  = max(anpha1v,0);
    anpha2_new  = max(anpha2v,0);
    anphaAN_new = max(anphaANv,0);

    % track objective value (cvx_optval is the surrogate objective value)
    obj_hist(k) = cvx_optval;
    anpha1_k = anpha1_new;
    anpha2_k = anpha2_new;
    anphaAN_k = anphaAN_new;
    obj_hist = obj_hist(1:k);
    % convergence check
    if norm([anpha1_new-anpha1_k; anpha2_new-anpha2_k; anphaAN_new-anphaAN_k], 2) <= tol
       break;
    end
end
%==============================
%outputs
%==============================
out.anpha1 = anpha1_k;
out.anpha2 = anpha2_k;
out.anphaAN = anphaAN_k;
out.tau=tau;
out.iters=k;
out.cvx_status=cvx_status_last;
out.obj_hist = obj_hist;

% Also output the "true" (non-surrogate) secrecy of U1 at the final powers (useful for Dinkelbach update)
% R1_true  = log2( (a1*(anpha1+anpha2)+sig1) / (a1*anpha2 + sig1) );
% RE1_true = log2( (aE*(anpha1+anpha2)+sigE+alpha*anphaAN) / (aE*anpha2 + sigE+alpha*anphaAN) );
% out.Rs1_true = max(R1_true - RE1_true, 0);

end

clear;clc;close all;
M=32;
% distances (m)
rho=0.7;


h_BS = 15;     %# chieu cao tram goc (m)
dx_IRS = 25;   %# toa Ä‘o X cua IRS (m)
h_IRS = 5;   %# chieu cao IRS (m)

dx1 = 30;     % # toa Ä‘o ngang nguoi dung 1 (m)
dx2 = 27;     % # toa Ä‘o ngang nguoi dung 2 (m)
dxE = 29;      %# toa Ä‘o ngang may nghe len (m)

%khoang cach tu BS den IRS
d_BI  = sqrt(((h_BS - h_IRS)^2 + dx_IRS^2));

%Khoang cach tu BS den cac nguoi dung 1,2 va eavesdropper
d_BU1 = sqrt(h_BS^2 + dx1^2);   % U1 far
d_BU2 = sqrt(h_BS^2 + dx2^2) ;    % U2 near
d_BUE  = sqrt(dxE^2 + h_BS^2); %Eavesdropper

%Khoang cach tu IRS den U1,U2,E phu thuoc vao vi tri tuong doi

if dx_IRS < dx2
    d_IU1 = sqrt((-dx_IRS + dx1)^2 + h_IRS^2);
    d_IU2 = sqrt((-dx_IRS + dx2)^2 + h_IRS^2);
    d_IE = sqrt((-dx_IRS + dxE)^2 + h_IRS^2);
elseif dx_IRS < dx1 && dx_IRS >= dx2
    d_IU1 = sqrt((-dx_IRS + dx1)^2 + h_IRS^2);
    d_IU2 = sqrt((dx_IRS - dx2)^2 + h_IRS^2);
    d_IE = sqrt((-dx_IRS + dxE)^2 + h_IRS^2);
elseif dx_IRS >= dxE
    d_IU1 = sqrt((dx_IRS - dx1)^2 + h_IRS^2);
    d_IU2 = sqrt((dx_IRS - dx2)^2 + h_IRS^2);
    d_IE = sqrt((dx_IRS - dxE)^2 + h_IRS^2);
else
    %# TrÆ°á»?ng há»£p cÃ²n láº¡i (IRS náº±m giá»¯a dx2 vÃ  dxE)
    d_IU1 = sqrt((dx_IRS - dx1)^2 + h_IRS^2);
    d_IU2 = sqrt((dx_IRS - dx2)^2 + h_IRS^2);
    d_IE = sqrt((dx_IRS - dxE)^2 + h_IRS^2); 
end


noise_awgn_dBm = -65;
sigma_awgn = 10^(noise_awgn_dBm / 10) * 1e-3;  %# cong suat nhieu (W)
sigma1=sigma_awgn;
sigma2=sigma_awgn;
sigmaE=sigma_awgn;

alpha=4; %jammer anh huong E nhu the nao

Pmax_set=[20 30 40];%dBm
SEE_prop=zeros(length(Pmax_set),1);%proposed model
SEE_noRIS=zeros(length(Pmax_set),1);%noRIS
SEE_rand=zeros(length(Pmax_set),1);%ris phase shift is random
SEE_ssr=zeros(length(Pmax_set),1);% maximizing SSR instead of SEE
%epsP=P0;
sig1 = sigma1 ;     
sig2 = sigma2 ;
sigE = sigmaE;   
NBSa=sig1;



epsLeak=1;%bao mat tot 1%bao mat vua, 5 -bao mat long
parCh.rho=rho;
parCh.M = M;%SO Phan tu IRS


parCh.d0    = d_BI;%bs-ris
parCh.d1    = d_BU1;%BS-U1
parCh.d2    = d_BU2;%BS-U2
parCh.dE    = d_BUE;%BS-UE
parCh.dIRS1 = d_IU1;%RIS-U1
parCh.dIRS2 = d_IU2;%RIS-U2
parCh.dRSE  = d_IE;%RIS-UE

% path loss exponents
parCh.m1 = 2.5;
parCh.m2 = 3.8;

% Rician coefficient
parGenCh.M=M;
parGenCh.K_BR = 5;
parGenCh.K_RU = 8;
parGenCh.K_RE = 8;
parGenCh.seed=1;
% parCh.K_BE = 2;

parCh.NBSa=NBSa;
parCh.N2a=sig2;
maxDink = 10;
maxAO=20;
tolA = 1e-4;
tolD=1e-4;
tolSCA= 1e-4;
tolPGA=1e-3;
maxSCA=5; %5-10-20
maxPGA=50;%M<50; M=64-100 -->80; m>=128 -->100

etaLeak=10;%penalty for leakage in Block B 10-20 khuyen cao
gamma0=0.5;
gamma1=0.5;
gamma2=0.5;


P_cU1=10 ^ (10 / 10) / 1000;
P_cU2=10 ^ (10 / 10) / 1000;
P_cBS=10 ^ (15 / 10) / 1000;
Pc_R=10^ (15 / 10) / 1000;
Pris=10 ^ (1 / 10) / 1000;%#1dBm cÃ´ng suat tieu hao mach Ä‘ien tren mot phan tu IRS
P_M_ris=M*Pris;
Pc=P_cU1+P_cU2+P_cBS+P_M_ris; 
Ns=10;
%tham so Dinklbach
beta=0.4;
lsMax=15;%max line search step
mu0=0.15;%initial step size

SEE_pmax=zeros(length(Pmax_set),1);
RE_worst_pmax=zeros(length(Pmax_set),1);
SEE_AO_samples=zeros(Ns,length(Pmax_set));
RE_AO_samples=zeros(Ns,length(Pmax_set));
for i=1:length(Pmax_set)
    SEE_hist={}; 
    RE_worst_hist={};
    epsilon=10^(Pmax_set(i)/10)*10^(-3);% doi sang watt
    fprintf('Pmax=%d\n',epsilon); 
    P0=epsilon;
    parCh.P0=P0;
    epsP=P0;
    % SEE_hist_prop = zeros(Ns,1);
    % SEE_hist_noRIS = zeros(Ns,1);
    % SEE_hist_rand = zeros(Ns,1);
    % SEE_hist_ssr = zeros(Ns,1);
    for itSample=1:Ns
        % ====== generate channels ======
        fprintf('itSample=%d\n',itSample);
        ch=generate_channel_SR(parGenCh);
    %     fprintf('ch=%d',ch);
        %%tinh het cac tham so tu ham nay
        parCh.g1=ch.g1;
        parCh.g2=ch.g2;
        parCh.gE=ch.gE;
        %fprintf('parCh.g1=%d, g2=%d,gE=%d',parCh.g1,parCh.g2,parCh.gE);
        parCh.hB1=ch.hB1;
        parCh.hB2=ch.hB2;
        parCh.hBE=ch.hBE;
        phi = exp(1j*2*pi*rand(M,1));   % random initial RIS phases--> dam bao tinh unit-modulus
        %fprintf('phi=%d',phi);

        phi0=phi;
        %parCh=ch;
        parCh.phi=phi;
        parCh.hBI=ch.hBI;


        coef=v2channelRIS_SR(parCh);
        parCommon = struct();
        parCommon.A1    = coef.A1;
        %fprintf('\n coef.A1=%d',coef.A1);
        parCommon.B1    = coef.B1;
        parCommon.A2    = coef.A2;
        parCommon.B2    = coef.B2;
        parCommon.C2    = coef.C2;
        parCommon.S1    = coef.S1;
        parCommon.S1bar = coef.S1bar;



     % ---- noise & jammer ----
        parCommon.sig1  = sig1;
        parCommon.sig2  = sig2;
        parCommon.sigE  = sigE;
        parCommon.alpha = alpha;   % jammer-to-E coefficient
        t=0; % Dinkelbach parameter
        xi=1;
    % ---- power consumption model ----
        parCommon.xi   = xi;        % PA inefficiency
        parCommon.Pc   = Pc;        % circuit power
        parCommon.Mris = M;         % number of RIS elements
        parCommon.Pris = Pris;      % power per RIS element
        anpha1  =(0.5+0.5*rand);
        anpha2  = (0.01+0.19*rand);
        anphaAN = 1-anpha1-anpha2;

        parCh.phi = exp(1j*2*pi*rand(parCh.M,1));   % random initial RIS phases
        solan=0;

        for itD = 1:maxDink
          improveTol  = 1e-8;     % objective acceptance tolerance
          F_prev = -inf;
            for itAO = 1:maxAO
                solan=solan+1;

                %fprintf(' solan=%g,.anpha1=%g,anpha2=%g,anphaAN=%g\n',solan,anpha1,anpha2,anphaAN);
                % ===================== Block A (power) =====================
                anpha1_old=anpha1; 
                anpha2_old=anpha2;
                anphaAN_old=anphaAN;    %             
                phi_old=phi; 


                parA=parCommon;

                parA.phi = phi; 
                parA.epsLeak=epsLeak;          % gioi han do di cho phep
                parA.t = t;
                parA.hBI=ch.hBI;
                parA.g1=ch.g1;
                parA.g2= ch.g2;
                parA.gE= ch.gE;
                parA.hB1=ch.hB1;
                parA.hB2=ch.hB2;
                parA.hBE=ch.hBE;

                parA.gamma0=gamma0;
                parA.gamma1=gamma1;
                parA.gamma2=gamma2;

                parA.init.anpha1=anpha1;
                parA.init.anpha2=anpha2;
                parA.init.anphaAN=anphaAN;



                parA.alpha=alpha;
                parA.epsP=epsP;
                parA.P0=P0;

                parA.sig1=sig1; 
                parA.sig2=sig2;
                parA.sigE=sigE;          


                parA.maxIter=maxSCA;
                parA.tol=tolSCA;

                parA.rho=rho;
                parA.Pc=Pc;

                outA = hong_opt_power_sca_secrecy(parA);
                anpha1 = outA.anpha1; 
                anpha2 = outA.anpha2;
                anphaAN = outA.anphaAN;
                %fprintf('sau BlockA: anpha1=%g, anpha 2=%g, anphaAN=%g\n',anpha1, anpha2, anphaAN);
                % ===================== Block B (phase) =====================
                %parB=struct();
                parB=parA;
                parB.phi0 = phi;
                parB.anpha1=anpha1;
                parB.anpha2=anpha2;
                parB.anphaAN=anphaAN;

                parB.hBI=ch.hBI;
                parB.g1=ch.g1;
                parB.g2= ch.g2;
                parB.gE= ch.gE;
                parB.hB1=ch.hB1;
                parB.hB2=ch.hB2;
                parB.hBE=ch.hBE;


                parB.d0    = d_BI;%bs-ris
                parB.d1    = d_BU1;%BS-U1
                parB.d2    = d_BU2;%BS-U2
                parB.dE    = d_BUE;%BS-UE
                parB.dIRS1 = d_IU1;%RIS-U1
                parB.dIRS2 = d_IU2;%RIS-U2
                parB.dRSE  = d_IE;%RIS-UE

                parB.A1    = coef.A1;
                parB.B1    = coef.B1;
                parB.A2    = coef.A2;
                parB.B2    = coef.B2;
                parB.C2    = coef.C2;
                parB.S1    = coef.S1;
                parB.S1bar = coef.S1bar;

                parB.m1 = 2.2;
                parB.m2 = 3.0;
                parB.beta=beta;
                parB.lsMax=lsMax;
   
                parB.alpha=alpha;
                parB.xi=xi;

                parB.Mris=M;
                parB.Pris=Pris;
                parB.mu0=mu0;
                parB.maxIter=maxPGA;
                parB.tol=tolPGA;
                parB.rho=rho;
                parB.Pc=Pc;
                parB.epsLeak=epsLeak;
                parB.etaLeak=etaLeak;
                parB.P0=P0;
                parB.N1a=sig1;
                parB.N2a=sig2;
                parB.NEa=sigE;
                parB.NBSa=sig1;


                outB = v2hong_opt_phase_pga_secrecy(parB);

                phi = outB.phi; 

    
                %----update parCh: phi ----
                parCh.phi=phi;

    %             fprintf('Sau Block A and B: anpha1=%g, anpha 2=%g, anphaAN=%g\n',anpha1, anpha2, anphaAN);
    %             fprintf('phi=%f,\n phi_old=%f',phi,phi_old);
                coefnew=v2channelRIS_SR(parCh);
                % ===================== Evaluate TRUE objective =====================

                parEval = parCommon;
                parEval.phi=phi; 
                parEval.rho=rho;
                parEval.P0=P0;
                parEval.N1a=sig1;
                parEval.N2a=sig2;
                parEval.NEa=sigE;
                parEval.NBSa=sig1;

                parEval.anpha1=anpha1;
                parEval.anpha2=anpha2;
                parEval.anphaAN=anphaAN;
                parEval.Pc=Pc;

                parEval.A1    = coefnew.A1;
                parEval.B1    = coefnew.B1;
                parEval.A2    = coefnew.A2;
                parEval.B2    = coefnew.B2;
                parEval.C2    = coefnew.C2;
                parEval.S1    = coefnew.S1;
                parEval.S1bar = coefnew.S1bar;
    %Chua tinh den su thay doi cua cua lamda tai day
                met = v2hong_eval_true_objective(parEval);

                % Dinkelbach inner objective value:
                F = met.SSR - t*met.Ptot;                
                % AO stopping: improvement small or variable change small
                dp = norm([anpha1-anpha1_old; anpha2-anpha2_old; anphaAN-anphaAN_old],2);
                dphi = norm(phi-phi_old,2);
                if itAO>1
                    if abs(F-F_prev) <= tolA* max (1, abs(F_prev))
                        %fprintf('AO converged at iter %d\n', itAO);
                        break;
                    end 
                end
                if dp <= 1e-4 && dphi <= 1e-3
                    %fprintf('check dp and dphi AO converged at iter %d\n', itAO);
                    break;
                end 

                F_prev = F;       

            end
            % ===================== Update Dinkelbach =====================           
            F_D=met.SSR-t*met.Ptot;
            %residual_hist(itSample, itD) = abs(F_D);
            t_new = met.SSR / met.Ptot;            
            %Dinkelbach stopping
            if abs(F_D) <= tolD                
                SEE_D=met.SEE;
                break;
            else                
                t=met.SSR/met.Ptot;
            end   
            t_old=t;
            if met.SSR<1e-8
               t=t_old;
            else
               t=t_new;
            end

        end
        outEval=v2hong_eval_true_objective(parEval);%tinh gia tri
        SEE_now=outEval.SEE;
        RE_worst_now=outEval.REworst;

        % RE_worst_hist.append(RE_worst_now);
        % SEE_pmax_hist.append(SEE_now);

        SEE_AO_samples(itSample, i) = SEE_now;
        RE_AO_samples(itSample, i)  = RE_worst_now;
    end
    SEE_AO_mean = mean(SEE_AO_samples, 1);
    SEE_AO_std  = std(SEE_AO_samples, 0, 1);
    
    RE_AO_mean = mean(RE_AO_samples, 1);
    RE_AO_std  = std(RE_AO_samples, 0, 1);
    
    save('AO_points_vs_Pmax.mat', ...
        'Pmax_set', ...
        'SEE_AO_samples', 'RE_AO_samples', ...
        'SEE_AO_mean', 'SEE_AO_std', ...
        'RE_AO_mean', 'RE_AO_std');

end

% figure;
% plot(Pmax_set, SEE_pmax, '-o', 'LineWidth', 1.5);
% grid on;
% 
% xlabel('P_{max}');
% ylabel('Average SEE');
% title('SEE versus transmit power');
% %%%%luu du lieu%%%%
% 
% save('SEE_Pmax_Ns_50result_dinkelbach.mat','SEE_prop','Ns','Pmax_set')
% disp('Simulation DONE & saved.');




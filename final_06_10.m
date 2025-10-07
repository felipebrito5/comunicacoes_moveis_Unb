% frequencia da portadora em GHz
fc_GHz = 3;
% componentes de multipercurso
N = 100;

% para UMi LoS
% tabela 1 (estatísticas de log10(sigma_tau [s]))
media = -0.24*log10(1+fc_GHz) - 7.14;       
desvio_padrao = 0.38;                          

% tabela 2
f_proporcionalidade = 3.0;                     

% gerando uma amostra aleatoria 
amostra_ale = normrnd(media, desvio_padrao, [N,1]);

% convertendo amostra aleatoria para segundos (sigma_tau)
sigma_tau_s = 10.^amostra_ale;                 

% para media da media 
media_tau = f_proporcionalidade*sigma_tau_s;

% gerando componentes exponenciais 
tau = exprnd(media_tau, [N,1]);                        

% desvio-padrão do sombreamento (dB) para UMi LoS
sigma_db = 4;
% termos de sombreamento por percurso (larga escala), em dB
xi_db = normrnd(0, sigma_db, [N,1]);

% normalizando e ordenando os atrasos
tau_orden = sort(tau - min(tau));
histogram(tau_orden);

% fator de Rice em dB (média 9 dB, desvio 5 dB para UMi LoS)
KR_dB = normrnd(9,5);
% converte K_R para escala linear
KR = 10^(KR_dB/10);

% potência difusa preliminar por percurso (n>=2): decaimento exponencial no atraso × sombreamento
alpha_hat2 = exp(-((f_proporcionalidade-1)/f_proporcionalidade).*(tau_orden ./ sigma_tau_s)) .* 10.^(sigma_db/10);
         

% Soma da parte difusa preliminar (exclui LoS n=1)
Omega_hat = sum(alpha_hat2(2:end));
if Omega_hat <= eps
    alpha_hat2(2:end) = 1/(N-1);                   
    Omega_hat = 1;
end

% aloca vetor final de potências normalizadas
alpha2 = zeros(N,1);
% potência do percurso LoS 
alpha2(1) = KR/(KR+1);     
% difusa normalizada p/ somar 1/(KR+1)
alpha2(2:end) = (alpha_hat2(2:end)./Omega_hat) / (KR+1); 

% ganho total do canal (≈1 após normalização)
Omega_c = sum(alpha2);

% checagem: razão entre LoS e soma da difusa (deve ≈ K_R)
KR_check = alpha2(1)/sum(alpha2(2:end));

tau_us = 1e6*tau_orden;

%item 3

% potências relativas (evita log(0))
ratio = max(eps, alpha2) ./ max(alpha2);

% ---------- Azimute (θ_n) ----------
% Estatísticas de σ_theta (Tabela 5, UMi LoS)
% média (log10 de graus)
mu_sth_log  = -0.08*log10(1+fc_GHz) + 1.73;  
% desvio-padrão (log10 de graus)
sig_sth_log =  0.014*log10(1+fc_GHz) + 0.28;  
% amostra de σθ;log
sth_log  = normrnd(mu_sth_log, sig_sth_log); 
% σθ em graus (escala linear)
sth_deg  = 10^(sth_log); 
% σθ em radianos
sth_rad  = deg2rad(sth_deg);                        

% Ângulos iniciais de azimute (em rad) equação 8
theta_p  = 1.42 * sth_rad .* sqrt( -log(ratio) );

% Sinais e flutuações
% U_n ∈ {-1,1}
Utheta   = randsample([-1,1], N, true);            
Ytheta   = normrnd(0, sth_rad/7, [N,1]);            

% Azimute final e ajuste LoS (θ1 = 0)
theta    = Utheta(:).*theta_p(:) + Ytheta;        
theta    = theta - theta(1);                       

% ---------- Elevação (φ_n) ----------
% Estatísticas de σ_phi (Tabela 6, UMi LoS)
% média (log10 de graus)
mu_sph_log  = -0.1*log10(1+fc_GHz) + 0.73;
% desvio-padrão (log10 de graus)
sig_sph_log = -0.04*log10(1+fc_GHz) + 0.34;       
sph_log  = normrnd(mu_sph_log, sig_sph_log); 
% σϕ em graus (escala linear)
sph_deg  = 10^(sph_log);    
% σϕ em radianos
sph_rad  = deg2rad(sph_deg);                        

% Ângulos iniciais de elevação (em rad)
phi_p   = -sph_rad .* log(ratio);                   

% Sinais e flutuações
Uphi    = randsample([-1,1], N, true);             
Yphi    = normrnd(0, sph_rad/7, [N,1]);            

% Elevação final e ajuste LoS (define φ̄ arbitrário em [0, π/2])
phi_bar = pi/6;                                    
phi     = Uphi(:).*phi_p(:) + Yphi; 
% LoS centrado em φ̄
phi     = phi - phi(1) + phi_bar;                   

% ---- visualizar espectros angulos de chegada azimute ----
theta_deg = rad2deg(theta);
y = alpha2;    
% evita zeros na escala log
y(y<=0) = min(y(y>0))*0.1;              

% topo baseado apenas na difusa (exclui n=1), com 20% de folga
y_top = 1.2 * max(y(2:end));
yl = [1e-8, y_top];

figure; hold on
% stems pretos só da difusa
s = stem(theta_deg(2:end), y(2:end), 'Color',[0 0 0], 'LineStyle','-');
s.Marker = '^'; s.MarkerSize = 6; s.MarkerFaceColor = [0 0 0]; s.MarkerEdgeColor = [0 0 0];

% LoS como seta azul sendo vizada
plot([0 0], yl, 'b', 'LineWidth', 2);
plot(0, yl(2), '^', 'MarkerSize', 9, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');

set(gca,'YScale','log','YLim',yl,'XLim',[-110 110],'YTick',10.^(-8:0));
grid on; grid minor
xlabel('Ângulos de chegada em azimute (°)')
ylabel('Potência')
title('PDP por azimute')

% (opcional) anotar o quanto o LoS é maior que a difusa
delta_dB = 10*log10(y(1)/max(y(2:end)));

hold off

% ---- visualizar espectros angulos de chegada elevação ----
% elevação em graus
phi_deg = rad2deg(phi);        
y = alpha2; 
y(y<=0) = min(y(y>0))*0.1;     

figure; hold on
% stems pretos só para a difusa (exclui o tap 1 / LoS)
s = stem(phi_deg(2:end), y(2:end), 'Color',[0 0 0], 'LineStyle','-');
s.Marker = '^'; s.MarkerSize = 6; s.MarkerFaceColor = [0 0 0]; s.MarkerEdgeColor = [0 0 0];

% seta azul no ângulo de elevação do LoS
yl  = [1e-8 1];                 
xL  = [0 90];                   
plot([phi_deg(1) phi_deg(1)], yl, 'b', 'LineWidth', 2);
plot(phi_deg(1), yl(2), '^', 'MarkerSize', 9, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');

set(gca,'YScale','log','YLim',yl,'XLim',xL,'YTick',10.^(-8:0));
grid on; grid minor
xlabel('Ângulos de chegada em elevação (°)')
ylabel('Potência')
title('PDP por elevação')
hold off



 
% converte elevação de rad 
phi_deg   = rad2deg(phi);     

% garante azimute em [0, 2π)
th = mod(theta, 2*pi);         
r  = phi_deg;                 

% --- Figura polar: LoS como seta azul, demais raios em vermelho ---
figure
pax = polaraxes; hold on

% limites e marcações do gráfico polar (0 a 90 graus de elevação)
rlim([0 90]);                  
rticks([0 20 40 60 80]);       
thetaticks(0:30:330);         

% seta azul do caminho LoS 
polarplot([th(1) th(1)], [0 r(1)], 'b', 'LineWidth', 1.8);
polarplot(th(1), r(1), 'bo', 'LineWidth', 1.8, 'MarkerFaceColor','none');

% espalha os demais multipercursos como pontos vermelhos
% (tamanho fixo; abaixo há uma versão ponderada pela potência)
polarscatter(th(2:end), r(2:end), 22, 'r', 'filled');

title('Direções de chegada (azimute \theta vs elevação \phi)')
grid on; hold off

% item 4
% --- Comprimento de onda na sua fc ---
c  = 3e8;                        
fc = fc_GHz*1e9;                   
lambda = c/fc;                     

% --- Vetores r_n a partir de theta, phi (já obtidos na Q3) ---
ux = cos(theta).*sin(phi);
uy = sin(theta).*sin(phi);
uz = cos(phi);
R  = [ux uy uz];                   % Nx3 com cada linha = r_n

% --- (a) Dois vetores velocidade: v=5 e 50 m/s ---
v_speeds = [5, 50];                


% e AZIMUTE ortogonal ao caminho LoS para que a componente LoS tenha ν≈0.
% 90° -> movimento no plano x-y
phi_v   = pi/2;      
% 90° em relação ao azimute do LoS
theta_v = theta(1) + pi/2;         
vhat = [cos(theta_v)*sin(phi_v); ...
        sin(theta_v)*sin(phi_v); ...
        cos(phi_v)];             

% --- (b) Cálculo dos desvios Doppler ν_n para cada velocidade ---
proj = R * vhat;                  

nu_5  = (v_speeds(1)/lambda) * proj;  
nu_50 = (v_speeds(2)/lambda) * proj; 

% --- (c) Espectros de potência no domínio Doppler (dois casos) ---
y = alpha2; 
y(y<=0) = min(y(y>0))*0.1;        

fDmax = v_speeds./lambda;          
yl = [1e-8 1];                      

figure
subplot(1,2,1)
s = stem(nu_5(2:end),  y(2:end),  'Color',[0 0 0], 'LineStyle','-'); hold on
s.Marker='^'; s.MarkerSize=6; s.MarkerFaceColor=[0 0 0]; s.MarkerEdgeColor=[0 0 0];
plot([0 0], yl, 'b','LineWidth',1.6)                                % “seta” no zero
set(gca,'YScale','log','YLim',yl,'XLim',[-fDmax(1) fDmax(1)]);
grid on; grid minor
xlabel('\nu (Hz)'); ylabel('Potência'); title('v_{rx} = 5 m/s')

subplot(1,2,2)
s = stem(nu_50(2:end), y(2:end), 'Color',[0 0 0], 'LineStyle','-'); hold on
s.Marker='^'; s.MarkerSize=6; s.MarkerFaceColor=[0 0 0]; s.MarkerEdgeColor=[0 0 0];
plot([0 0], yl, 'b','LineWidth',1.6)
set(gca,'YScale','log','YLim',yl,'XLim',[-fDmax(2) fDmax(2)]);
grid on; grid minor
xlabel('\nu (Hz)'); ylabel('Potência'); title('v_{rx} = 50 m/s')

%----------item 5 
% escolhe uma direção de velocidade para o receptor
    v_rx = 5;  
    % movimento no plano horizontal
    phi_v = pi/2;   
    % ortogonal ao LoS -> nu(LoS) ~ 0
    theta_v = theta(1)+pi/2;  
    vhat = [cos(theta_v)*sin(phi_v); sin(theta_v)*sin(phi_v); cos(phi_v)];
    nu = (v_rx/lambda) * (R*vhat);   % ν_n (Hz)  [N x 1]
    
    % -------------------- Q5(a): Pulsos retangulares --------------------
dt_list = [1e-7, 1e-5, 1e-3]; 
% nº de amostras no eixo do tempo para cada pulso
Nt = 1e5;                    

% Função pulso retangular unitário 
% s(t; dt) = 1 em [0, dt], e 0 fora
rect = @(t,dt) double((t>=0) & (t<=dt));

% -------------------- Q5(b): Sinal recebido para cada δt --------------------
results = cell(numel(dt_list),1);
for k = 1:numel(dt_list)
    dt = dt_list(k);
    t = linspace(0, 5*dt, Nt).';              
    r = complex(zeros(Nt,1));                  

    % termo fixo de fase por percurso 
    phase0 = exp(-1j*2*pi*(((fc_GHz*10^9)+nu).*tau_orden));  

    % soma multipercurso 
    for n = 1:numel(alpha2)
        % máscara temporal onde s(t-τ_n) = 1
        mask = (t >= tau_orden(n)) & (t <= (tau_orden(n)+dt));
        if any(mask)
            r(mask) = r(mask) + alpha2(n) * phase0(n) .* exp(1j*2*pi*nu(n).*t(mask));
        end
    end

    results{k} = struct('t',t,'r',r,'dt',dt);
end

% -------------------- Q5(c): Plota os três casos no mesmo canal --------------------
figure

if ~exist('sigma_tau','var')
    Omega_c = sum(alpha2);
    tau_bar = (alpha2.'*tau_orden)/Omega_c;
    sigma_tau = sqrt((alpha2.'*((tau_orden - tau_bar).^2))/Omega_c);  % [s]
end

for k = 1:numel(results)
    % eixo temporal [s] em [0, 5*dt]
    t  = results{k}.t;  
    % sinal recebido complexo 
    r  = results{k}.r; 
    % largura do pulso 
    dt = results{k}.dt;          

    % pulso retangular unitário: s(t)=1 em [0, dt], 0 fora
    s_tx = double((t>=0) & (t<=dt));

    % largura de banda aproximada do pulso (para o título)
    Bw = 1/dt;                  
    
    figure; hold on
    % pulso transmitido (step/retângulo azul)
    stairs(t, s_tx, 'b', 'LineWidth', 1.6);

    % sinal recebido (módulo em vermelho)
    plot(t, abs(r), 'r', 'LineWidth', 1.4);

    % estética parecida com a figura de referência
    grid on
    xlim([0, 5*dt]);
    ylim([0, 1.2*max([1; abs(r)])]);  
    xlabel('Tempo absoluto — t (s)')
    ylabel('|{\itr}\_{}(t)|','Interpreter','tex')  
    title(sprintf('Sinal Recebido, \\delta t = %.0e s,  B_w \\approx %.0e Hz,  \\sigma_\\tau \\approx %.0f ns', ...
                  dt, Bw, 1e9*sigma_tau), 'Interpreter','tex')
    legend({'Sinal Transmitido','Sinal Recebido'}, 'Location','northeast')
    hold off
end




%-----

figure
s = stem(tau_us, alpha2, 'Color',[0 0 0], 'LineStyle','-'); hold on
s.Marker = '^';
s.MarkerSize = 5;
s.MarkerFaceColor = [0 0 0];
s.MarkerEdgeColor = [0 0 0];
s.BaseLine.Visible = 'on';
s.BaseLine.BaseValue = 1e-8;

ax = gca;
ax.YScale = 'log';
ax.YLim   = [1e-8 1];
ax.XLim   = [0 0.8];
ax.YTick  = 10.^(-8:0);
grid on; grid minor

xlabel('Domínio de Atraso — \tau (\mus)')
ylabel('PDP')
title('Potência Multipercurso')

plot([0 0], ax.YLim, 'b', 'LineWidth', 1.2)  
hold off

% (g) atraso médio e espalhamento de atraso (RMS) com base em alpha2 e tau
tau_bar   = (alpha2.' * tau_orden) / Omega_c;                   % \bar{\tau}
sigma_tau = sqrt( (alpha2.' * ((tau_orden - tau_bar).^2)) / Omega_c );  % \sigma_\tau

% imprime métricas-chave
disp(struct( ...
    'KR_dB', KR_dB, ...
    'KR', KR, ...
    'KR_check', KR_check, ...
    'Omega_c', Omega_c, ...
    'tau_bar_us', 1e6*tau_bar, ...
    'sigma_tau_us', 1e6*sigma_tau ...
))
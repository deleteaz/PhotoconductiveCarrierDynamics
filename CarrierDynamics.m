clc;clear;close all;
% 双峰0.7eV and 0.9eV 推荐参数（截至2025年10月）
% Ep = 0.8;       % [eV]
% E0 = 1.5;       % [eV]
% Ea = 1.5;       % [eV]
% N0 = 1e17;      % 
% phi0 = 1.0;     % [eV]
% phi1 = 0.18;    % [eV]
% phi2 = 0.06;    % [eV]
% Tc = 30         % [摄氏度]
% ------------------------------------------------ %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 氧空位模型演化（常温探究）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 基础参数
mu = 1;              % [cm^2/(V*s)]
S  = 100e-7 * 10e-1; % [cm^2]
e  = 1.60217662e-19; % [C]
I  = 1e-9;           % [A]
epsilon = 200;       % [V/cm]
% 载流子浓度
n0 = I / (mu * e * epsilon * S) * 1; % [cm^-3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 可调参数
Epp = 0.10;     % 2价氧空位深度 [eV]
Ep  = 0.85;     % 1价氧空位深度 [eV]
E0  = 1.28;     % 0*价氧空位深度[eV]
Es  = 1.28;     % 施主深度      [eV]

phi_pp = 0.06;  % 2价到1价激活能 [eV] % 现为一种概率因子
phi_p  = 0.115;  % 1价到0*价激活能[eV] % 现为一种概率因子
phi_0  = 1.00;  % 0*价到0价激活能[eV] % 现为一种概率因子

sigma_pp = 0.2; % 2价半高宽
sigma_p  = 0.2; % 1价半高宽
sigma_0  = 0.2; % 0*价半高宽
sigma_s  = 0.2; % 施主半高宽

Ns = 1e17;      % 施主浓度
Tc = 30;        % [摄氏度]
Erange = 0.5;   % E+ΔE的ΔE，控制积分上下限
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for E_p = [0.75,0.80,0.85,0.90]
ETot = {Epp, Ep, E0, Es};
phiTot = {phi_pp, phi_p, phi_0};
sigmaTot = {sigma_pp,sigma_p,sigma_0,sigma_s};
tspan = [0, 5e2];
y0 = [0; n0/2; 0; n0];
options = odeset('RelTol',1e-5,'AbsTol',1e-5,'MaxStep', 1e5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,y] = draw_output(ETot, phiTot, sigmaTot, Ns, Tc, Erange, tspan, y0, options);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算主程序
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy_dt = odefun1(t, y, ETot, phiTot, sigmaTot, Ns, Tc, Erange)
k = 8.617333262145e-5; % [eV/K]
T = Tc + 273.15;       % [K]
nu0 = 1e13;            % [Hz]
me = 9.10956e-31;      % [kg]
m = me * 0.34;         % [kg]
S = 9e-20;             % [cm^2]
ve = sqrt((8 * k * T * 1.60217662e-19) / (pi * m)) * 1e2;   %[cm/s]
C_pp = ve * S;         % [cm^3/s]
C_p = 0.5 * C_pp;

[Epp, Ep, E0, Es] = ETot{:};
[phi_pp, phi_p, phi_0] = phiTot{:};
[sigma_pp, sigma_p, sigma_0, sigma_s] = sigmaTot{:};

% 四种类型载流子浓度：[一价氧空位，二价氧空位，零价*氧空位，电子]
Np  = y(1);
Npp = y(2);
N0  = y(3);
n   = y(4);

% 四种类型载流子高斯分布：[一价氧空位，二价氧空位，零价*氧空位，施主氧空位]
rho_p  = @(E) 1 / (sqrt(pi) * sigma_p) * exp(-(E - Ep).^2 / sigma_p^2);
rho_pp = @(E) 1 / (sqrt(pi) * sigma_pp) * exp(-(E - Epp).^2 / sigma_pp^2);
rho_0  = @(E) 1 / (sqrt(pi) * sigma_0) * exp(-(E - E0).^2 / sigma_0^2);
rho_s  = @(E) 1 / (sqrt(pi) * sigma_s) * exp(-(E - Es).^2 / sigma_s^2);

% 自适应积分
% 四种类型载流子贡献（向导带）：[一价氧空位，二价氧空位(x)，零价*氧空位，施主氧空位]
Np_up  = integral(@(E) Np  * rho_p(E)  * nu0 * exp(-Ep / k / T), Ep-Erange, Ep+Erange);
% Npp_up = integral(@(E) Npp * rho_pp(E) * nu0 * exp(-Epp / k / T), Epp-Erange, Epp+Erange);
N0_up  = integral(@(E) N0  * rho_0(E)  * nu0 * exp(-E0 / k / T), E0-Erange, E0+Erange);
Ns_up  = integral(@(E) Ns  * rho_s(E)  * nu0 * exp(-Es / k / T), Es-Erange, Es+Erange);

% 四种类型载流子贡献（向价带）：[一价氧空位，二价氧空位，零价*氧空位，施主氧空位(x)]
Np_dw  = integral(@(E) Np  * rho_p(E)  * n * C_p * exp(-phi_p / k / T) , Ep-Erange, Ep+Erange);
Npp_dw = integral(@(E) Npp * rho_pp(E) * n * C_pp * exp(-phi_pp / k / T) , Epp-Erange, Epp+Erange);
N0_dw  = integral(@(E) N0  * rho_0(E)  * nu0 * exp(-phi_0 / k / T) , E0-Erange, E0+Erange);
% Ns_dw  = integral(@(E) Ns  * rho_s(E)  * n * C_s * exp(-phi_s / k / T) , Es-Erange, Es+Erange);

% 四种类型载流子浓度随时间变化：[一价氧空位，二价氧空位，零价*氧空位，电子]
dNp_dt  = Npp_dw + N0_up - Np_up - Np_dw;
dNpp_dt = Np_up - Npp_dw;
dN0_dt  = Np_dw - N0_up - N0_dw;
dn_dt   = Np_up + N0_up + Ns_up - Npp_dw - Np_dw;

dy_dt = [dNp_dt; dNpp_dt; dN0_dt; dn_dt];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 绘图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,y] = draw_output(ETot, phiTot, sigmaTot, Ns, Tc, Erange, tspan, y0, options)
hlen = 1;
y = cell(hlen,1);
t = cell(hlen,1);
h1 = cell(hlen,1); h2 = cell(hlen,1); h3 = cell(hlen,1); h4 = cell(hlen,1);

[Epp, Ep, E0, Es] = ETot{:};
[phi_pp, phi_p, phi_0] = phiTot{:};

figure()
for hid = 1:hlen
    if hid == 1
        % y = [Np; Npp; N0; n]
        fun = @(t,y) odefun1(t, [y(1); y(2); y(3); y(1)+2*y(2)], ETot, phiTot, sigmaTot, Ns, Tc, Erange);
        [t{1}, y{1}] = ode45(fun, tspan, y0, options);
    end
    % subplot(1,1,hid)
    hold on
    ydata = y{hid};
    tdata = t{hid};
    delete("data1.xlsx");
    M = [tdata(:), ydata(:,4)];
    writematrix(M, ['data',num2str(hid),'.xlsx']);
    h1{hid} = plot(tdata, ydata(:,1), "LineWidth",1.2,"LineStyle","--","Color",[0.9,0.2,0.3]);
    h2{hid} = plot(tdata, ydata(:,2), "LineWidth",1.2,"LineStyle","--","Color",[0.9,0.5,0.3]);
    h3{hid} = plot(tdata, ydata(:,3), "LineWidth",1.2,"LineStyle","--","Color",[0.3,0.5,0.5]);
    h4{hid} = plot(tdata, ydata(:,4), "LineWidth",1.2,"LineStyle","-","Color",[0.3,0.5,0.9]);
    % maxIndices = islocalmax(ydata(:,4));
    % maxIndices = find(maxIndices > 0);
    % plot(tdata(maxIndices), ydata(maxIndices,4), "Marker","o","LineStyle","none","Color",[0.3,0.5,0.9])
    if hid == 1
        ylabel("载流子浓度n/cm^{-3}")
    end
    if hid == hlen
        legend([h1{hid},h2{hid},h3{hid},h4{hid}], ["N_{T}^{+}(t)","N_{T}^{++}(t)","N_0^{\ast}(t)","n(t)"])
    end
    if hid == 1
        title(sprintf("恒温(%.1f°C)",Tc))
    end
    % yyaxis("right")
    % t1 = ydata(2:end,4) - ydata(1:end-1,4);
    % plot(tdata(1:end-1), t1 > 0, "LineWidth",1.2,"LineStyle","-","Color",[0.5,0.7,0.5])
    % ylabel("是否后项大于前项")
    % ylim([-1,2])

    xlabel("时间t/s")
    text(0.1,0.7,sprintf("E_{+} = %.2feV \nE_{++} = %.2feV \nE_{0} = %.2feV \nE_{s} = %.2feV" + ...
        "\nφ(N_{0}^{*}->N_{0}) = %.2feV" + ...
        "\nφ(N^{+}->N_{0}^{*}) = %.2feV" + ...
        "\nφ(N^{++}->N^{+}) = %.2feV",...
        Ep,Epp,E0,Es,phi_0,phi_p,phi_pp),"FontWeight","bold","FontSize",11,"Units","normalized")
    set(gca, "LineWidth",1.2, "Box","on", "FontWeight","bold", "FontSize",15)
end
run("DLT2Energy.mlx")
drawnow()
end

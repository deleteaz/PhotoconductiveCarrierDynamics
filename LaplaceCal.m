%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 拉普拉斯变换寻峰算法
% !!!!!注意!!!!!
% 需手动调整最后的高斯拟合参数，以适应不同数量级导致的拟合不准的问题
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;clear;close all;warning off
% 常数
theta = 1e13; % 弛豫时间常数
N = 1e5+1;      % 插值数
kT = 0.026;
% 电流前因子
S = 1e-2 * 100e-9; % [m*m]
mu = 1e-4; % [m^2/(Vs)]
e = 1.6e-19; % [C]
epsilon = 20e3; % [V/m]
A = 1/(S * mu * e * epsilon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 数据导入
file = "data1.xlsx";
data = readmatrix(file);
t_data = data(:, 1); % Time data
I_data = data(:, 2); % Intensity data
% t_data = linspace(0,1000,N);
% tau1 = 10; tau2 = 200; %%
% I_data = 1 * exp(-t_data/tau1) + 1 * exp(-t_data/tau2);
t_data = t_data(:);
I_data = I_data(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 数据插值（提高稳健性，降低噪声干扰）
t_min = min(t_data);
t_max = max(t_data);
t_fine = linspace(t_min, t_max, N)';
I_fine = interp1(t_data, I_data, t_fine, 'spline');
% 统一向量取向
t_fine = t_fine(:);
I_fine = I_fine(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求解域定义
s_values = logspace(-5,2.05,N)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 拉普拉斯积分
parallel_num = 1; % 并行计算
integral_results = zeros(size(s_values));
for i = 1:parallel_num:length(s_values)
    if (length(s_values) - i) <= parallel_num
        s = s_values(i:end);
        integrand_fine = I_fine .* exp(-s' .* t_fine);
        integral_results(i:end) = cotes(t_fine, integrand_fine);
        break
    end
    s = s_values(i:i+parallel_num-1);
    integrand_fine = I_fine .* exp(-s' .* t_fine);
    integral_results(i:i+parallel_num-1) = cotes(t_fine, integrand_fine);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求导
s_hatn = s_values .* integral_results;    % s*hat{n}
lns = log(s_values);                      % lns
spl = spline(lns, s_hatn);
ds_hatn_dlns = ppval(fnder(spl, 1), lns); % d(s*hat{n})/d(lns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最终结果
g = (1 / kT) .* ds_hatn_dlns;
E = kT * log(theta ./ s_values);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 数据可视化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 原始数据
% fig1 = figure();
% hold on
% plot(t_data, I_data,"LineWidth",1.5);
% xlabel("t/s", "Interpreter","latex")
% ylabel("$n/cm^{-3}$", "Interpreter","latex")
% set(gca,"FontSize", 15)
% print(fig1, "原始数据", "-dpng");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 中间数据
% fig2 = figure();
% subplot(2,1,1)
% plot(s_values, integral_results,"LineWidth",1.5);
% xlabel("s", "Interpreter","latex")
% ylabel("拉普拉斯积分值")
% set(gca,"FontSize", 15)
% subplot(2,1,2)
% plot(s_values(1:end), s_hatn,"LineWidth",1.5);
% xlabel("s", "Interpreter","latex")
% ylabel("$s\times\hat{n}$", "Interpreter","latex")
% set(gca,"FontSize", 10)
% print(fig2, "拉普拉斯积分值与sn值", "-dpng");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求导结果
% fig3 = figure();
% hold on
% maxIndices = islocalmax(ds_hatn_dlns, 1);
% maxIndices = find(maxIndices > 0);
% plot(s_values, ds_hatn_dlns,"LineWidth",1.5)
% plot(s_values(maxIndices), ds_hatn_dlns(maxIndices), "Marker","o","LineStyle","none")
% for i = 1:length(maxIndices)
%     text(s_values(maxIndices(i)), ds_hatn_dlns(maxIndices(i)), sprintf("s = %.4f\ntau = %.2f", ...
%         s_values(maxIndices(i)), 1/s_values(maxIndices(i))),...
%         "FontSize",12,"FontWeight","bold");
% end
% xlabel("s", "Interpreter","latex")
% ylabel("$\frac{d(s\times\hat{n})}{d(lns)}$", "Interpreter","latex")
% set(gca,"FontSize", 15)
% print(fig3, "d(sn)d(lns)值", "-dpng");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 最终数据
fig4 = figure();
subplot(1,2,1)
hold on
maxIndices = islocalmax(g, 1);
maxIndices = find(maxIndices > 0);
h1 = plot(E, g,"LineWidth",2);
plot(E(maxIndices), g(maxIndices), "Marker","o","LineStyle","none")
for i = 1:length(maxIndices)
    text(E(maxIndices(i)), g(maxIndices(i)), sprintf("E = %.3f\ns = %.4f\ntau = %.2f", ...
        E(maxIndices(i)), 1/(exp(E(maxIndices(i))/kT)/theta), exp(E(maxIndices(i))/kT)/theta), ...
        "FontSize",12,"FontWeight","bold");
end
% 高斯峰拟合（非双峰需自行添加）
% 对比 g(E) 与 叠加峰，越贴近，拟合越好
LB = [min(g), 0.60, 1e-3, min(g), 0.60, 1e-3, min(g), 0.60, 1e-3];
UB = [max(g), 1.10, 1e0, max(g), 1.10, 1e0, max(g), 1.10, 1e0];
[gfit, pfit] = JADEfit(E, g, LB, UB, "off");
% h2 = plot(E, gfit,"LineWidth",1.5);
E1 = linspace(0,2,5e4)';
gmax = cell(length(LB)/3 ,1);
idmax = cell(length(LB)/3 ,1);
Area = cell(length(LB)/3 ,1);
hfit = cell(length(LB)/3 ,1);
for hfit_id = 1:3:length(LB)
    hid = (hfit_id+2)/3;
    [gmax{hid}, idmax{hid}] = max(gaussFun(E1, pfit(hfit_id), pfit(hfit_id+1), pfit(hfit_id+2)));
    Area{hid} = integral(@(E1) gaussFun(E1, pfit(hfit_id), pfit(hfit_id+1), pfit(hfit_id+2)), min(E1), max(E1));
    hfit{hid} = plot(E, gaussFun(E, pfit(hfit_id), pfit(hfit_id+1), pfit(hfit_id+2)),"LineStyle","--");
    text(E1(idmax{hid}),gmax{hid}/2,sprintf("Area%d = %.1e",hid,Area{hid}), "FontSize",12,"FontWeight","bold");
end
xlabel('E/eV', "Interpreter","latex");
ylabel('$g(E)/(m^{-3}\times eV^{-1})$', "Interpreter","latex")
% legend([h1,h2,h3,h4],["g(E)","叠加峰","拟合峰1","拟合峰2"])
% legend("g(E)")
set(gca,"FontSize", 15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 高斯峰拟合（非双峰需自行添加）
subplot(1,2,2)
hold on
h1 = plot(E, g,"LineWidth",2);
plot(E(maxIndices), g(maxIndices), "Marker","o","LineStyle","none")
% for i = 1:length(maxIndices)
%     text(E(maxIndices(i)), g(maxIndices(i)), sprintf("E = %.3f\ns = %.4f\ntau = %.2f", ...
%         E(maxIndices(i)), 1/(exp(E(maxIndices(i))/kT)/theta), exp(E(maxIndices(i))/kT)/theta), ...
%         "FontSize",12,"FontWeight","bold");
% end
h2 = plot(E, gfit,"LineWidth",1.5);
for hfit_id = 1:3:length(LB)
    hid = (hfit_id+2)/3;
    [gmax{hid}, idmax{hid}] = max(gaussFun(E1, pfit(hfit_id), pfit(hfit_id+1), pfit(hfit_id+2)));
    Area{hid} = integral(@(E1) gaussFun(E1, pfit(hfit_id), pfit(hfit_id+1), pfit(hfit_id+2)), min(E1), max(E1));
    hfit{hid} = plot(E, gaussFun(E, pfit(hfit_id), pfit(hfit_id+1), pfit(hfit_id+2)),"LineStyle","--");
    text(E1(idmax{hid}),gmax{hid}/2,sprintf("E%d = %.2e\ns%d = %.2e\ntau%d = %.2e",...
        hid,E1(idmax{hid}),...
        hid,1/exp(E1(idmax{hid})/kT)/theta,...
        hid,(exp(E1(idmax{hid})/kT)/theta)),...
        "FontSize",12,"FontWeight","bold");
end

xlabel('E/eV', "Interpreter","latex");
% ylabel('$g(E)/(m^{-3}\times eV^{-1})$', "Interpreter","latex")
% legend([h1,h2,h3,h4],["g(E)","叠加峰","拟合峰1","拟合峰2"])
set(gca,"FontSize", 15)
% print(fig4, "g(E)_E", "-dpng");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 原始数据拟合(非双峰需自行添加e指数)
% fig5 = figure();
% hold on
% n_data = I_data;
% n_fit = 0;
% for i = 1:length(LB)/3
%     n_fit = n_fit + Area{i} * exp(-t_data ./ ((exp(E1(idmax{i})/kT)/theta)));
% end
% % n_fit = 1 * exp(-t_data ./ (1/(exp(E(maxIndices(1))/kT)/theta))) + 1 * exp(-t_data ./ (1/(exp(E(maxIndices(2))/kT)/theta)));
% % n_fit = 1e19 * exp(-t_data ./ 30) + 1e19 * exp(-t_data ./ 1);
% h1 = plot(t_data, n_data, "LineWidth",1.5);
% h2 = plot(t_data, n_fit);
% legend([h1,h2],["原始数据","反推结果"])
% xlabel("t/s", "Interpreter","latex")
% ylabel("$n/cm^{-3}$", "Interpreter","latex")
% set(gca,"FontSize", 15, "XScale","linear", "YScale","linear")
% print(fig5, "能量逆推原始数据", "-dpng");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 数据导出
% str = {'拉普拉斯积分值','s*n','d(sn)/d(lns)','g（E）','E','逆推I'};
% type = {'double', 'double', 'double', 'double', 'double', 'double'};
% sz = [length(g),length(str)];
% T = table('Size',sz,'VariableTypes',type,'VariableNames',str);
% T{:,"d(sn)/d(lns)"} = ds_hatn_dlns;
% T{:,"s*n"} = s_hatn;
% T{:,"拉普拉斯积分值"} = integral_results;
% T{:,"g（E）"} = g;
% T{:,"E"} = E;
% % T{:,"逆推I"} = n_fit;
% writetable(T, "结果数据.xlsx");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 复化cotes积分法
function I = cotes(x, y)
% 单行计算
% n = size(x,1);
% if mod(n-1, 4) ~= 0
%     error('区间数量必须为4的倍数');
% end
% h = x(2) - x(1);
% sum_1_4 = sum(y(2:4:end-3, :), 1);
% sum_1_2 = sum(y(3:4:end-2, :), 1);
% sum_3_4 = sum(y(4:4:end-1, :), 1);
% sum_mid = sum(y(5:4:end-4, :), 1);
% I = (7*(y(1) + y(end)) + 32*sum_1_4 + 12*sum_1_2 + 32*sum_3_4 + 14*sum_mid);
% I = I * h * 4 / 90;

% 并行计算
n = length(x);  % 改用长度（假设x为向量）
if mod(n-1, 4) ~= 0  % 修正边界检查：子区间数为4的倍数
    error('子区间数量必须为4的倍数（节点数-1）');
end
h = x(2) - x(1);
m = size(y, 2);  % 被积函数数量（列数）
I = zeros(1, m);  % 初始化结果

I = I + 7 * (y(1,:) + y(end,:));

% 系数32的位置：x(2), x(6), x(10)...
idx = 2:4:n-3;
for k = idx
    I = I + 32 * y(k,:);
end
% 系数12的位置：x(3), x(7), x(11)...
idx = 3:4:n-2;
for k = idx
    I = I + 12 * y(k,:);
end
% 系数32的位置：x(4), x(8), x(12)...
idx = 4:4:n-1;
for k = idx
    I = I + 32 * y(k,:);
end
% 系数14的位置：x(5), x(9), x(13)...
idx = 5:4:n-4;
for k = idx
    I = I + 14 * y(k,:);
end
I = I * h * 4 / 90;  % 最终系数
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = LorentzFun(x, a, b, c)
y = a ./ ((x - b).^2 + c);
end

function y = gaussFun(x, a, b, c)
y = a .* exp(-(x - b).^2 / c^2);
end

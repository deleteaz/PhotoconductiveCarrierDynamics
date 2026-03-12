function [fitting_gbest, pop_gbest] = JADEfit(x, yTrue, LB, UB, display)
% 求解器初始化
POP_SIZE = 50;  % 求解器个数
DIM1 = length(LB);

% 交变参数
muCR = 0.5;
muSF = 0.5;
P = 0.05;
C = 0.1;
memoryA = [];

% 种群初始化
pop = LB + (UB - LB) .* rand(POP_SIZE, DIM1);
fitness = inf .* ones(POP_SIZE, 1);
fitting = zeros(length(yTrue), POP_SIZE);
fitness_gbest = inf;
record = [];
GER = 400;
IGER = 0;

while GER > IGER
    successCR = [];
    successSF = [];

    rand_CR = normrnd(muCR, 0.1, POP_SIZE, 1);
    if any(rand_CR > 1)
        rand_CR(rand_CR > 1) = 1;
    end
    if any(rand_CR < 0)
        rand_CR(rand_CR < 0) = 0;
    end

    rand_SF = cauchyrnd(muSF, 0.1, POP_SIZE, 1);
    while any(rand_SF <= 0)
        rand_SF(rand_SF <= 0) = cauchyrnd(muSF, 0.1, 1, 1);
    end
    if any(rand_SF > 1)
        rand_SF(rand_SF > 1) = 1;
    end

    pbestSize = P * POP_SIZE;

    r1 = randi([1, POP_SIZE], POP_SIZE, 1);
    r2 = randi([1, POP_SIZE], POP_SIZE, 1);

    for i = 1:POP_SIZE
        while r1(i) == i
            r1(i) = randi([1, POP_SIZE], 1);
        end
        while any([r2(i) == i, r2(i) == r1(i)])
            r2(i) = randi([1, POP_SIZE + size(memoryA, 1)], 1);
        end
    end

    popMemory = cat(1, pop, memoryA);
    pbest = ceil(rand(POP_SIZE, 1) .* pbestSize);
    pop_mutat = pop + ...
        rand_SF .* (pop(pbest, :) - pop) + ...
        rand_SF .* (pop(r1, :) - popMemory(r2, :));

    j = (1:DIM1) .* ones(POP_SIZE, DIM1);
    jrand = randi([0, DIM1], POP_SIZE, DIM1);
    rc = rand(POP_SIZE, DIM1);
    cross_id = (j == jrand) + (rc < muCR);
    pop_cross = cross_id .* pop_mutat + ~cross_id .* pop;

    LB_id = pop_cross < LB;
    UB_id = pop_cross > UB;
    LB_reset = (LB + pop) / 2;
    UB_reset = (UB + pop) / 2;
    pop_cross = LB_id .* LB_reset + ~LB_id .* pop_cross;
    pop_cross = UB_id .* UB_reset + ~UB_id .* pop_cross;

    for i = 1:POP_SIZE
        [fitting(:, i), fitness_cross] = fobj(pop_cross(i, :), x, yTrue);
        if fitness_cross < fitness(i)
            successCR = cat(1, successCR, rand_CR(i));
            successSF = cat(1, successSF, rand_SF(i));

            pop(i, :) = pop_cross(i, :);
            fitness(i) = fitness_cross;
            memoryA = cat(1, memoryA, pop(i, :));

            if fitness_cross < fitness_gbest
                fitness_gbest = fitness_cross;
                pop_gbest = pop(i, :);
                fitting_gbest = fitting(:, i);
            end
        end
    end

    while size(memoryA, 1) > POP_SIZE
        memoryA(randi([1, POP_SIZE], 1), :) = [];
    end

    if size(successCR, 1) > 0
        meanA = mean(successCR);
        meanL = sum(successSF .^ 2) ./ sum(successSF);
        muCR = (1 - C) * muCR + C * meanA;
        muSF = (1 - C) * muSF + C * meanL;
    end

    [fitness, best_id] = sort(fitness);
    pop = pop(best_id, :);

    record = [record, fitness_gbest];
    IGER = IGER + 1;
end

if display == "on"
    % 结果可视化
    % 迭代曲线
    figure()
    plot(record, "color", [0.1, 0.3, 0.7], "LineWidth", 1.1)
    yline(record(end), "LineWidth", 1.1, "LineStyle", "--", "Color", [0.1, 0.1, 0.1])
    text(GER-15, record(end), sprintf("最优MSE为:%.3e", record(end)), "color", [0.8, 0.4, 0.3])
    title("迭代最优曲线", "FontSize", 15, "FontWeight", "bold")
    xlabel("迭代次数")
    ylabel("MSE")
    set(gca, "FontSize", 12, "FontWeight", "bold", "LineWidth", 1.5);

    % 拟合结果
    figure()
    hold on
    plot(x, fitting_gbest, 'r--')
    plot(x, yTrue, 'b-')
    title("拟合结果", "FontSize", 15, "FontWeight", "bold")
    xlabel("x")
    ylabel("y")
    set(gca, "FontSize", 12, "FontWeight", "bold", "LineWidth", 1.5);
    hold off
end

end

function y = LorentzFun(x, a, b, c)
    y = a ./ ((x - b).^2 + c);
end

function y = gaussFun(x, a, b, c)
    y = a .* exp(-(x - b).^2 / c^2);
end

function [yTest, y] = fobj(param, x, yTrue)
% 求解函数
yTest = 0;
% for i = 1:3:length(param)
%     yTest = yTest + LorentzFun(x, param(i), param(i+1), param(i+2));
% end
for i = 1:3:length(param)
    yTest = yTest + gaussFun(x, param(i), param(i+1), param(i+2));
end
% 目标函数(R^2)
y = sqrt(sum((yTrue - yTest).^2) / length(yTrue));
end

function x = cauchyrnd(mu, sigma, varargin)
% 柯西随机数
% 输入: mu sigma m*n
% 输出: 随机数
x = mu + sigma .* tan(pi * (rand(varargin{:}) - 0.5));
end


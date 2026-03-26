clear; clc; close all;

%% 1. 系统基础参数设定
N = 16;
N_DFnT = 16;
L = 2;              
Block_Num = 2;    
C = 2;              
loop_Num = 100000;   % 循环次数，为了平滑可以设大一点（比如50000），跑得慢可调小
M = 2;              % 统一使用 2-DCSK 
Equal = 2;          % MMSE 均衡

bits_per_symbol = log2(M);

% 使用 k=9 的参数作为基准测试 (你也可以根据需要改成 k=5)
c1_dcsk = 5/(2*N_DFnT); 
c2_dcsk = 1.414; 

%% 2. 扫频参数设定 (对齐你的原版 csi.m)
rho_vec = logspace(-2, 0, 11); % 生成从 0.01 到 1 的 11 个对数分布点
EbN0_dB_vec = [6, 12, 18];     % 测试三个固定的信噪比

% 预分配矩阵保存结果 (3行 11列)
ber_dcsk_matrix = zeros(length(EbN0_dB_vec), length(rho_vec));

fprintf('开始自动化生成 DCSK-AFDM (k=9) 的 CSI 鲁棒性图表数据...\n');
fprintf('总计需要运行 %d 个组合。\n\n', length(EbN0_dB_vec) * length(rho_vec));

%% 3. 核心仿真双重循环
for i = 1:length(EbN0_dB_vec)
    dB = EbN0_dB_vec(i);
    EbN0_linear = 10^(dB/10);
    
    % 【公平性折算】计算 DCSK 的物理 SNR
    SNR_dcsk_linear = EbN0_linear * (bits_per_symbol / (2 * N));
    
    fprintf('=== 正在仿真 Eb/N0 = %d dB ===\n', dB);
    
    for j = 1:length(rho_vec)
        rou_current = rho_vec(j);
        err_dcsk = 0;
        
        for loop = 1:loop_Num
            % 发射端
            [Bits_d, Sym0_d] = Transmitter_MDCSK(N_DFnT, M, N, C, c1_dcsk, c2_dcsk);
            % 信道 (理想传输)
            [h_d, Sym1_d]    = Channel(Sym0_d, L, N, Block_Num, SNR_dcsk_linear);
            % 接收端 (注入当前 rou_current 的估计误差)
            Bitsre_d         = Receiver_MDCSK_rou(rou_current, N_DFnT, M, N, C, Equal, Sym1_d, h_d, SNR_dcsk_linear, c1_dcsk, c2_dcsk);
            
            % 统计错误比特
            err_dcsk = err_dcsk + sum(Bits_d ~= Bitsre_d);
        end
        
        % 计算当前 (SNR, rho) 组合的 BER
        ber_dcsk_matrix(i, j) = err_dcsk / (loop_Num * bits_per_symbol);
        fprintf('  rho = %.4f | BER = %.6f\n', rou_current, ber_dcsk_matrix(i, j));
    end
    fprintf('\n');
end

%% 4. 自动打印提取出的硬编码数组 (方便你复制到最终论文脚本里保存)
fprintf('====================================================\n');
fprintf('仿真完成！请复制以下代码以便日后直接画图：\n\n');
fprintf('x_rho = logspace(-2, 0, 11);\n');
fprintf('y6_dcsk = [%s];\n', num2str(ber_dcsk_matrix(1, :), '%g '));
fprintf('y12_dcsk = [%s];\n', num2str(ber_dcsk_matrix(2, :), '%g '));
fprintf('y18_dcsk = [%s];\n', num2str(ber_dcsk_matrix(3, :), '%g '));
fprintf('====================================================\n\n');

%% 5. 直接绘制图表
figure('Name', 'DCSK-AFDM BER vs CSI Error', 'Position', [200, 200, 700, 500]);
box on; hold on; grid on;

% 绘制三条曲线，使用不同的标记区分
plot(rho_vec, ber_dcsk_matrix(1,:), 'g-s', 'LineWidth', 1.5, 'MarkerSize', 7);
plot(rho_vec, ber_dcsk_matrix(2,:), 'm-s', 'LineWidth', 1.5, 'MarkerSize', 7);
plot(rho_vec, ber_dcsk_matrix(3,:), 'r-s', 'LineWidth', 1.5, 'MarkerSize', 7);

% 对数坐标轴设置
set(gca,'Yscale','log');
set(gca,'Xscale','log');

% 坐标轴范围与你原版 csi.m 完全一致
ylim([1e-5 1]); 
xlim([1e-2, 1e0]); 

% 刻度标签
xlabel('(\rho)', 'FontSize', 12);
ylabel('BER', 'FontSize', 12);
title('M-DCSK-AFDM Performance under Imperfect CSI', 'FontSize', 13);

legend('DCSK-AFDM (6 dB)', 'DCSK-AFDM (12 dB)', 'DCSK-AFDM (18 dB)', 'Location', 'southwest', 'FontSize', 11);
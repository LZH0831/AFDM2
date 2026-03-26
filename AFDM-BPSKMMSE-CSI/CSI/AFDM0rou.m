clear; clc; close all;

%% 1. 系统基础参数设定
N = 16;
N_DFnT = 16;
L = 2;              % 2径频率选择性信道 (论文假设)
Block_Num = 2;    % 建议加大块数以获得更平滑的曲线
C = 2;              % 循环前缀长度
P = N + C;
loop_Num = 10000;    % 蒙特卡洛仿真次数 (可根据电脑算力适当增减)
M = 2;              % BPSK 调制
Equal = 2;          % 2 = MMSE 均衡器

rou = 0.01;         % 论文图3设定的 CSI 不确定性条件 ρ=0.01
SNR_dB_vec = 0:2:18;% 遍历 SNR 从 0 到 18dB

%% 2. 定义三种对比方案的 DFnT 参数
% 方案 1: OFDM
c1_ofdm = 0;             c2_ofdm = 0;
% 方案 2: AFDM (k=5)
c1_afdm5 = 5/(2*N_DFnT); c2_afdm5 = 1.414; 
% 方案 3: AFDM (k=9)
c1_afdm9 = 9/(2*N_DFnT); c2_afdm9 = 0.4; 

%% 3. 初始化结果矩阵
ber_ofdm  = zeros(1, length(SNR_dB_vec));
ber_afdm5 = zeros(1, length(SNR_dB_vec));
ber_afdm9 = zeros(1, length(SNR_dB_vec));

fprintf('开始仿真复现 Fig.3 (rho = %.2f, MMSE, BPSK)...\n', rou);

%% 4. 核心仿真循环
for idx = 1:length(SNR_dB_vec)
    dB = SNR_dB_vec(idx);
    SNR = 10^(dB/10);
    fprintf('正在计算 SNR = %2d dB... ', dB);
    
    err_o = 0; err_a5 = 0; err_a9 = 0;
    
    for loop = 1:loop_Num
        % === 方案1: OFDM ===
        [Bits_o, Sym0_o] = Transmitter(N_DFnT, M, Block_Num, N, C, c1_ofdm, c2_ofdm);
        [h_o, Sym1_o]    = Channel(Sym0_o, L, N, Block_Num, SNR);
        Bitsre_o         = Receiver(rou, N_DFnT, M, Block_Num, N, C, Equal, Sym1_o, h_o, SNR, c1_ofdm, c2_ofdm);
        err_o = err_o + sum(Bits_o ~= Bitsre_o);
        
        % === 方案2: AFDM (k=5) ===
        [Bits_a5, Sym0_a5] = Transmitter(N_DFnT, M, Block_Num, N, C, c1_afdm5, c2_afdm5);
        [h_a5, Sym1_a5]    = Channel(Sym0_a5, L, N, Block_Num, SNR);
        Bitsre_a5          = Receiver(rou, N_DFnT, M, Block_Num, N, C, Equal, Sym1_a5, h_a5, SNR, c1_afdm5, c2_afdm5);
        err_a5 = err_a5 + sum(Bits_a5 ~= Bitsre_a5);
        
        % === 方案3: AFDM (k=9) ===
        [Bits_a9, Sym0_a9] = Transmitter(N_DFnT, M, Block_Num, N, C, c1_afdm9, c2_afdm9);
        [h_a9, Sym1_a9]    = Channel(Sym0_a9, L, N, Block_Num, SNR);
        Bitsre_a9          = Receiver(rou, N_DFnT, M, Block_Num, N, C, Equal, Sym1_a9, h_a9, SNR, c1_afdm9, c2_afdm9);
        err_a9 = err_a9 + sum(Bits_a9 ~= Bitsre_a9);
    end
    
    % 计算平均误码率
    total_bits = loop_Num * Block_Num * N * log2(M);
    ber_ofdm(idx)  = err_o / total_bits;
    ber_afdm5(idx) = err_a5 / total_bits;
    ber_afdm9(idx) = err_a9 / total_bits;
    
    fprintf('完成\n');
end

%% 5. 绘图 (还原论文图表样式)
figure('Name', 'Reproduction of Fig. 3');
semilogy(SNR_dB_vec, ber_ofdm,  'b--o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
semilogy(SNR_dB_vec, ber_afdm9, 'r-o',  'LineWidth', 1.5, 'MarkerSize', 8);
semilogy(SNR_dB_vec, ber_afdm5, 'y-x', 'LineWidth', 1.5, 'MarkerSize', 8);

grid on; box on;
set(gca, 'YScale', 'log');
ylim([1e-5 1]); 
xlim([0 18]);
xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('OFDM', 'AFDM (k=9)', 'AFDM (k=5)', 'Location', 'southwest');
title(sprintf('Imperfect CSI (\\rho=%.2f)', rou));

%% 6. 保存数据结果到本地文件
myfilename = sprintf('Fig3_Results_rho=%.2f.txt', rou);
fid1 = fopen(myfilename, 'w');
fprintf(fid1, 'SNR(dB)\tOFDM\t\tAFDM(k=9)\tAFDM(k=5)\n');
for idx = 1:length(SNR_dB_vec)
    fprintf(fid1, '%d\t%e\t%e\t%e\n', SNR_dB_vec(idx), ber_ofdm(idx), ber_afdm9(idx), ber_afdm5(idx));
end
fclose(fid1);
disp(['仿真结束，数据已保存至：', myfilename]);
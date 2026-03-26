clear; clc; close all;

%% 1. 系统基础参数设定
N = 16;
N_DFnT = 16;
L = 2;              
Block_Num = 2;    
C = 2;              
P = N + C;
loop_Num = 50000;    
M = 2;              % 统一使用 BPSK / 2-DCSK 对比
Equal = 2;          

rou = 0.1;          % 不完美 CSI
EbN0_dB_vec = 0:2:24;

%% 2. 定义 DFnT 参数 (继承你的设定)
c1_afdm5 = 5/(2*N_DFnT); c2_afdm5 = 1.414; 
c1_afdm9 = 9/(2*N_DFnT); c2_afdm9 = 0.4; 

%% 3. 初始化结果矩阵
ber_psk_k5  = zeros(1, length(EbN0_dB_vec));
ber_psk_k9  = zeros(1, length(EbN0_dB_vec));
ber_dcsk_k5 = zeros(1, length(EbN0_dB_vec));
ber_dcsk_k9 = zeros(1, length(EbN0_dB_vec));

fprintf('开始综合对比仿真 (rho=%.2f, MMSE)...\n', rou);

%% 4. 核心仿真循环
bits_per_symbol = log2(M);

for idx = 1:length(EbN0_dB_vec)
    dB = EbN0_dB_vec(idx);
    EbN0_linear = 10^(dB/10);
    
    % 物理 SNR 折算
    SNR_psk_linear  = EbN0_linear * bits_per_symbol;
    SNR_dcsk_linear = EbN0_linear * (bits_per_symbol / (2 * N));
    
    fprintf('正在计算 Eb/N0 = %2d dB... ', dB);
    err_psk5 = 0; err_psk9 = 0; err_dcsk5 = 0; err_dcsk9 = 0;
    
    for loop = 1:loop_Num
        % === 方案1: PSK-AFDM (k=5) ===
        [Bits_p5, Sym0_p5] = Transmitter(N_DFnT, M, Block_Num, N, C, c1_afdm5, c2_afdm5);
        [h_p5, Sym1_p5]    = Channel(Sym0_p5, L, N, Block_Num, SNR_psk_linear);
        Bitsre_p5          = Receiver(rou, N_DFnT, M, Block_Num, N, C, Equal, Sym1_p5, h_p5, SNR_psk_linear, c1_afdm5, c2_afdm5);
        err_psk5 = err_psk5 + sum(Bits_p5 ~= Bitsre_p5);
        
        % === 方案2: PSK-AFDM (k=9) ===
        [Bits_p9, Sym0_p9] = Transmitter(N_DFnT, M, Block_Num, N, C, c1_afdm9, c2_afdm9);
        [h_p9, Sym1_p9]    = Channel(Sym0_p9, L, N, Block_Num, SNR_psk_linear);
        Bitsre_p9          = Receiver(rou, N_DFnT, M, Block_Num, N, C, Equal, Sym1_p9, h_p9, SNR_psk_linear, c1_afdm9, c2_afdm9);
        err_psk9 = err_psk9 + sum(Bits_p9 ~= Bitsre_p9);
        
        % === 方案3: DCSK-AFDM (k=5) ===
        [Bits_d5, Sym0_d5] = Transmitter_MDCSK(N_DFnT, M, N, C, c1_afdm5, c2_afdm5);
        [h_d5, Sym1_d5]    = Channel(Sym0_d5, L, N, Block_Num, SNR_dcsk_linear);
        Bitsre_d5          = Receiver_MDCSK_rou(rou, N_DFnT, M, N, C, Equal, Sym1_d5, h_d5, SNR_dcsk_linear, c1_afdm5, c2_afdm5);
        err_dcsk5 = err_dcsk5 + sum(Bits_d5 ~= Bitsre_d5);
        
        % === 方案4: DCSK-AFDM (k=9) ===
        [Bits_d9, Sym0_d9] = Transmitter_MDCSK(N_DFnT, M, N, C, c1_afdm9, c2_afdm9);
        [h_d9, Sym1_d9]    = Channel(Sym0_d9, L, N, Block_Num, SNR_dcsk_linear);
        Bitsre_d9          = Receiver_MDCSK_rou(rou, N_DFnT, M, N, C, Equal, Sym1_d9, h_d9, SNR_dcsk_linear, c1_afdm9, c2_afdm9);
        err_dcsk9 = err_dcsk9 + sum(Bits_d9 ~= Bitsre_d9);
    end
    
    total_bits_psk  = loop_Num * Block_Num * N * bits_per_symbol;
    total_bits_dcsk = loop_Num * bits_per_symbol;
    
    ber_psk_k5(idx)  = err_psk5 / total_bits_psk;
    ber_psk_k9(idx)  = err_psk9 / total_bits_psk;
    ber_dcsk_k5(idx) = err_dcsk5 / total_bits_dcsk;
    ber_dcsk_k9(idx) = err_dcsk9 / total_bits_dcsk;
    fprintf('完成\n');
end

%% 5. 绘图对比
figure('Name', 'Comprehensive Comparison: Fig 3');
semilogy(EbN0_dB_vec, ber_psk_k9,  'm-^', 'LineWidth', 1.5, 'MarkerSize', 7); hold on;
semilogy(EbN0_dB_vec, ber_psk_k5,  'r-o', 'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbN0_dB_vec, ber_dcsk_k9, 'm--^', 'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbN0_dB_vec, ber_dcsk_k5, 'r--o', 'LineWidth', 1.5, 'MarkerSize', 7);

grid on; box on; set(gca, 'YScale', 'log');
ylim([1e-5 1]); xlim([0 25]);
xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('PSK-AFDM (k=9)', 'PSK-AFDM (k=5)', 'DCSK-AFDM (k=9)', 'DCSK-AFDM (k=5)', 'Location', 'southwest');
title(sprintf('Imperfect CSI (\\rho=%.2f)', rou));
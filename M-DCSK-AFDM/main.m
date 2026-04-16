clear; clc; close all;
N = 16;             % 子载波数 (DCSK 扩频因子 beta=16)
N_DFnT = 16;
L = 2;              
Block_Num = 2;      
C = 2;              % 循环前缀长度
Equal = 2;          % 2代表使用 MMSE 均衡
loop_Num = 100000;   % 循环次数

EbN0_dB = 0:2:30;   % 横坐标统一为 Eb/N0 
M_vec = [2, 4];     % 2对应 BPSK/2-DCSK，4对应 QPSK/4-DCSK

total_ber_psk = zeros(length(M_vec), length(EbN0_dB));
total_ber_dcsk = zeros(length(M_vec), length(EbN0_dB));

fprintf('开始统一仿真：M-PSK-AFDM 与 M-DCSK-AFDM 对比 (公平 Eb/N0 约束, 理想 CSI)...\n');

for m_idx = 1:length(M_vec)
    M = M_vec(m_idx);
    bits_per_symbol = log2(M);
    
    for snr_idx = 1:length(EbN0_dB)
        dB = EbN0_dB(snr_idx);
        EbN0_linear = 10^(dB/10); 
        
        SNR_psk_linear = EbN0_linear * log2(M);
        SNR_dcsk_linear = EbN0_linear * (log2(M) / (2 * N));
        
        err_bits_psk = 0;
        err_bits_dcsk = 0;
        
        for loop = 1:loop_Num
            %% ================= 1. 原版 M-PSK-AFDM =================
            [Bits_psk, Symbols0_psk] = Transmitter(N_DFnT, M, Block_Num, N, C);
            [H0_psk, Symbols1_psk] = Channel(Symbols0_psk, L, N, Block_Num, SNR_psk_linear);
            Bitsre_psk = Receiver(N_DFnT, M, Block_Num, N, C, Equal, Symbols1_psk, H0_psk, SNR_psk_linear);
            err_bits_psk = err_bits_psk + sum(Bits_psk ~= Bitsre_psk);
            
            %% ================= 2. 融合 M-DCSK-AFDM =================
            [Bits_dcsk, Symbols0_dcsk] = Transmitter_MDCSK(N_DFnT, M, N, C);
            [H0_dcsk, Symbols1_dcsk] = Channel(Symbols0_dcsk, L, N, Block_Num, SNR_dcsk_linear);
            Bitsre_dcsk = Receiver_MDCSK(N_DFnT, M, N, C, Equal, Symbols1_dcsk, H0_dcsk, SNR_dcsk_linear);
            err_bits_dcsk = err_bits_dcsk + sum(Bits_dcsk ~= Bitsre_dcsk);
        end
        
        % 计算 BER 
        total_ber_psk(m_idx, snr_idx) = err_bits_psk / (loop_Num * N * Block_Num * bits_per_symbol);
        total_ber_dcsk(m_idx, snr_idx) = err_bits_dcsk / (loop_Num * bits_per_symbol);
        
        fprintf('M=%d, Eb/N0=%2d dB | PSK-AFDM BER: %.2e | DCSK-AFDM BER: %.2e\n', ...
            M, dB, total_ber_psk(m_idx, snr_idx), total_ber_dcsk(m_idx, snr_idx));
    end
end

%% 画图对比
figure('Position', [150, 150, 800, 600]);
hold on; grid on; box on;

semilogy(EbN0_dB, total_ber_psk(1,:), 'r>-', 'LineWidth', 1.5,  'MarkerSize', 6);
semilogy(EbN0_dB, total_ber_psk(2,:), 'rs-', 'LineWidth', 1.5,  'MarkerSize', 6);
semilogy(EbN0_dB, total_ber_dcsk(1,:), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbN0_dB, total_ber_dcsk(2,:), 'bd-', 'LineWidth', 1.5, 'MarkerSize', 7);

set(gca, 'YScale', 'log', 'FontSize', 11);
ylim([1e-5 1]);
xlim([0 30]);
xlabel('E_b/N_0 (dB)', 'FontSize', 12);
ylabel('BER', 'FontSize', 12);
title('PSK-AFDM vs DCSK-AFDM (Ideal CSI)', 'FontSize', 13);

legend('AFDM(k=9) BPSK', ...
       'AFDM(k=9) QPSK', ...
       '2-DCSK-AFDM(k=9)', ...
       '4-DCSK-AFDM(k=9)', ...
       'Location', 'southwest', 'FontSize', 11);
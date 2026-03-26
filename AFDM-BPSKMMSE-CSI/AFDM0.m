clear; clc; close all;

%% 1. 系统基础参数设定
N = 16;
N_DFnT = 16;
L = 2;              
Block_Num = 2;    
C = 2;              
P = N + C;
loop_Num = 10000;   % 平滑曲线推荐 20000 次以上
Equal = 2;          % 2代表使用 MMSE 均衡

% Fig 2 横坐标: 0 到 30 dB (考虑到运行时间，这里设为步长 2)
EbN0_dB_vec = 0:2:30; 

%% 2. 方案参数设定
% OFDM 的参数设定 (k=0 时，T1和T2退化为单位阵)
c1_ofdm = 0; 
c2_ofdm = 0;

% AFDM 的参数设定 (k=9)
c1_afdm = 9/(2*N_DFnT); 
c2_afdm = 1.414;

%% 3. 初始化结果保存数组
ber_ofdm_bpsk = zeros(1, length(EbN0_dB_vec));
ber_afdm_bpsk = zeros(1, length(EbN0_dB_vec));
ber_ofdm_qpsk = zeros(1, length(EbN0_dB_vec));
ber_afdm_qpsk = zeros(1, length(EbN0_dB_vec));

fprintf('开始复现 Figure 2 (Perfect CSI, MMSE)...\n');

%% 4. 核心仿真循环
for idx = 1:length(EbN0_dB_vec)
    dB = EbN0_dB_vec(idx);
    EbN0_linear = 10^(dB/10);
    fprintf('正在计算 Eb/N0 = %2d dB... ', dB);
    
    err_o_b = 0; err_a_b = 0;
    err_o_q = 0; err_a_q = 0;
    
    for loop = 1:loop_Num
        %% ----- BPSK 组 (M=2) -----
        M_bpsk = 2;
        SNR_bpsk_linear = EbN0_linear * log2(M_bpsk); 
        
        % 1. OFDM (k=0) BPSK
        [Bits_ob, Sym0_ob] = Transmitter(N_DFnT, M_bpsk, Block_Num, N, C, c1_ofdm, c2_ofdm);
        [H0_ob, Sym1_ob]   = Channel(Sym0_ob, L, N, Block_Num, SNR_bpsk_linear);
        Bitsre_ob          = Receiver(N_DFnT, M_bpsk, Block_Num, N, C, Equal, Sym1_ob, H0_ob, SNR_bpsk_linear, c1_ofdm, c2_ofdm);
        err_o_b = err_o_b + sum(Bits_ob ~= Bitsre_ob);
        
        % 2. AFDM (k=9) BPSK
        [Bits_ab, Sym0_ab] = Transmitter(N_DFnT, M_bpsk, Block_Num, N, C, c1_afdm, c2_afdm);
        [H0_ab, Sym1_ab]   = Channel(Sym0_ab, L, N, Block_Num, SNR_bpsk_linear);
        Bitsre_ab          = Receiver(N_DFnT, M_bpsk, Block_Num, N, C, Equal, Sym1_ab, H0_ab, SNR_bpsk_linear, c1_afdm, c2_afdm);
        err_a_b = err_a_b + sum(Bits_ab ~= Bitsre_ab);

        %% ----- QPSK 组 (M=4) -----
        M_qpsk = 4;
        SNR_qpsk_linear = EbN0_linear * log2(M_qpsk); % QPSK 能量补偿
        
        % 3. OFDM (k=0) QPSK
        [Bits_oq, Sym0_oq] = Transmitter(N_DFnT, M_qpsk, Block_Num, N, C, c1_ofdm, c2_ofdm);
        [H0_oq, Sym1_oq]   = Channel(Sym0_oq, L, N, Block_Num, SNR_qpsk_linear);
        Bitsre_oq          = Receiver(N_DFnT, M_qpsk, Block_Num, N, C, Equal, Sym1_oq, H0_oq, SNR_qpsk_linear, c1_ofdm, c2_ofdm);
        err_o_q = err_o_q + sum(Bits_oq ~= Bitsre_oq);
        
        % 4. AFDM (k=9) QPSK
        [Bits_aq, Sym0_aq] = Transmitter(N_DFnT, M_qpsk, Block_Num, N, C, c1_afdm, c2_afdm);
        [H0_aq, Sym1_aq]   = Channel(Sym0_aq, L, N, Block_Num, SNR_qpsk_linear);
        Bitsre_aq          = Receiver(N_DFnT, M_qpsk, Block_Num, N, C, Equal, Sym1_aq, H0_aq, SNR_qpsk_linear, c1_afdm, c2_afdm);
        err_a_q = err_a_q + sum(Bits_aq ~= Bitsre_aq);
    end
    
    total_bits_bpsk = loop_Num * Block_Num * N * log2(2);
    total_bits_qpsk = loop_Num * Block_Num * N * log2(4);
    
    ber_ofdm_bpsk(idx) = err_o_b / total_bits_bpsk;
    ber_afdm_bpsk(idx) = err_a_b / total_bits_bpsk;
    ber_ofdm_qpsk(idx) = err_o_q / total_bits_qpsk;
    ber_afdm_qpsk(idx) = err_a_q / total_bits_qpsk;
    
    fprintf('完成\n');
end

%% 5. 绘图复现
figure('Name', 'Fig 2 Reproduction');
box on; hold on; grid on;

% 严格对齐论文中图 2 的线型和颜色
semilogy(EbN0_dB_vec, ber_ofdm_qpsk, 'b--o', 'LineWidth', 1.5, 'MarkerSize', 7); 
semilogy(EbN0_dB_vec, ber_afdm_qpsk, 'b--d', 'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbN0_dB_vec, ber_ofdm_bpsk, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbN0_dB_vec, ber_afdm_bpsk, 'r->', 'LineWidth', 1.5, 'MarkerSize', 7);

set(gca, 'YScale', 'log', 'FontSize', 12);
ylim([1e-6 1]); 
xlim([0 30]);

xlabel('E_b/N_0 (dB)', 'FontSize', 14);
ylabel('BER', 'FontSize', 14);
title('perfect CSI over 2-path channel', 'FontSize', 14);

legend('OFDM(k=0) QPSK', ...
       'AFDM(k=9) QPSK', ...
       'OFDM(k=0) BPSK', ...
       'AFDM(k=9) BPSK', ...
       'Location', 'southwest', 'FontSize', 12);
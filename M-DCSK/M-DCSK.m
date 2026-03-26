% 复现论文 Fig. 3: M-ary DCSK 在 AWGN 信道下的 BER 性能曲线 (\beta = 60)
clear; clc; close all;

beta = 60;                  
num_symbols = 50000;        
EbN0_dB_range = 0:2:24;     
M_vec = [2, 4, 8, 16, 32];  

markers = {'o', 's', '^', 'v', 'd'};
colors = {'b', 'r', 'k', 'm', 'g'};

figure('Position', [100, 100, 700, 550]);
hold on; grid on;
set(gca, 'YScale', 'log'); 

fprintf('开始仿真 M-ary DCSK (扩频因子 beta=%d, 符号数=%d)...\n', beta, num_symbols);

for m_idx = 1:length(M_vec)
    M = M_vec(m_idx);
    bits_per_symbol = log2(M);
    BER_list = zeros(1, length(EbN0_dB_range));
    
    fprintf('正在仿真 %d-ary DCSK...\n', M);
    
    for snr_idx = 1:length(EbN0_dB_range)
        EbN0_dB = EbN0_dB_range(snr_idx);
        EbN0 = 10^(EbN0_dB / 10);
        Eb = 1.0;                  
        Es = bits_per_symbol * Eb; 
        N0 = Eb / EbN0;            
        sigma = sqrt(N0 / 2);      
        
        data = randi([0, M-1], num_symbols, 1);
        
        phase_offset = 0; 
        if M == 4
            phase_offset = pi/4;
        end
        sym = pskmod(data, M, phase_offset, 'gray');
        a_s = real(sym);
        b_s = imag(sym);
        
        x = rand(num_symbols, 1) * 2 - 1;
        cx = zeros(num_symbols, beta);
        for k = 1:beta
            cx(:, k) = x;
            x = 1 - 2 * x.^2; 
        end
        
        cx = cx - mean(cx, 2);
        cx = cx ./ sqrt(sum(cx.^2, 2));
        
        cx_analytic = hilbert(cx')'; 
        cy = imag(cx_analytic);
        cy = cy ./ sqrt(sum(cy.^2, 2)); 
        
        m_s = a_s .* cx + b_s .* cy;
        
        tx_ref  = sqrt(Es / 2) * cx;
        tx_info = sqrt(Es / 2) * m_s;
        

        rx_ref  = tx_ref + sigma * randn(num_symbols, beta);
        rx_info = tx_info + sigma * randn(num_symbols, beta);
        

        rx_ref_analytic = hilbert(rx_ref')';
        rx_ref_y = imag(rx_ref_analytic);
        

        z_a = sum(rx_info .* rx_ref, 2);
        z_b = sum(rx_info .* rx_ref_y, 2);
        z = z_a + 1i * z_b; 
        
        data_hat = pskdemod(z, M, phase_offset, 'gray');
        
        [~, ber] = biterr(data, data_hat);
        BER_list(snr_idx) = ber;
    end
    
    semilogy(EbN0_dB_range, BER_list, ['-', markers{m_idx}], ...
        'Color', colors{m_idx}, 'LineWidth', 1.5, 'MarkerSize', 7, ...
        'DisplayName', sprintf('%d-DCSK', M));
        
    fprintf('  完成! 在 %ddB 处的 BER 为: %.5e\n', EbN0_dB_range(end), BER_list(end));
end

set(gca, 'YMinorGrid', 'on', 'XMinorGrid', 'off', 'FontSize', 11);
xlabel('E_b/N_0 (dB)', 'FontSize', 13, 'Interpreter', 'tex');
ylabel('BER', 'FontSize', 13);
title('BER Performance of M-ary DCSK (\beta = 60)', 'FontSize', 14);
legend('Location', 'southwest', 'FontSize', 11);
ylim([1e-5, 1]);
xlim([0, 26]);
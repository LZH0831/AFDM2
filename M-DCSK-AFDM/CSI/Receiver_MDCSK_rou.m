function Bitsre = Receiver_MDCSK_rou(rou, N_DFnT, M, N, C, Equal, Symbols1, h, SNR, c1, c2)
    P = N + C;
    L = length(h); % 动态获取多径数 (1x2 向量的长度 L=2)
    Block_Num = 2;
    R = [zeros(N, P-N), eye(N)];
    
    %% 1. 去除循环前缀
    Symbols2 = zeros(N, 1, Block_Num);
    for a = 1:Block_Num
        Symbols2(:,:,a) = R * Symbols1(:,:,a);
    end
    
    %% 2. DAFT 矩阵构建 (使用外部传入的 c1, c2)
    W_IDFT=zeros(N_DFnT);
    for a=1:N_DFnT
       for b=1:N_DFnT
           W_IDFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N_DFnT);
       end 
    end
    T_1=zeros(N_DFnT);
    for a=1:N_DFnT
       for b=1:N_DFnT
           if a==b, T_1(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c1); end
       end 
    end
    T_1 = T_1*exp(-1i*pi/4);
    T_1H = conj(T_1);

    T_2=zeros(N_DFnT);
    for a=1:N_DFnT
        for b=1:N_DFnT
            if a==b, T_2(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c2); end
        end 
    end
    T_2 = T_2*1/sqrt(N_DFnT);
    T_2H = conj(T_2);

    Q_IDFnT0 = T_2H*W_IDFT*T_1H;
    P_DFnT0 = conj(Q_IDFnT0);
      
    Symbols3 = zeros(size(Symbols2));
    for count = 1:Block_Num
        Symbols3(:,:,count) = P_DFnT0 * Symbols2(:,:,count);
    end

    %% 3. ================= 核心修复：构造不完美 CSI 矩阵 =================
    % 1. 在真实的 1xL 冲激响应上注入误差
    err_h = (1/sqrt(2)) * (randn(size(h)) + 1i*randn(size(h)));
    hi = sqrt(1 - rou^2) * h + rou * err_h; % 得到不完美的信道向量 hi

    % 2. 将 1xL 的 hi 扩充为 PxP 的循环卷积矩阵 H0_est
    H0_est = zeros(P); 
    for a = 1:P
        for b = 1:P                 
            if a-b < 0 || a-b > L-1     
                H0_est(a,b) = 0;
            else
                H0_est(a,b) = hi(a-b+1);
            end
        end
    end
    %% ====================================================================

    %% 4. 构造 MMSE 均衡矩阵
    S = eye(N);
    T = [S(2*N-P+1:N,:); S];
    H_est = R * H0_est * T;
    D_est = P_DFnT0 * H_est * Q_IDFnT0; 

    if Equal == 2
       G = pinv((D_est)' * D_est + eye(size(D_est)).*(1/SNR)) * (D_est)';
    else
       G = inv(D_est); 
    end
    
    %% 5. 均衡与 M-DCSK 解调
    Symbols4 = zeros(size(Symbols3));
    for count = 1:Block_Num
        Symbols4(:,:,count) = G * Symbols3(:,:,count);
    end

    rx_ref  = real(Symbols4(:,:,1));
    rx_info = real(Symbols4(:,:,2));
    
    rx_ref_analytic = hilbert(rx_ref);
    rx_ref_y = imag(rx_ref_analytic);
    
    z_a = sum(rx_info .* rx_ref);
    z_b = sum(rx_info .* rx_ref_y);
    z = z_a + 1i * z_b; 
    
    phase_offset = 0; 
    if M == 4, phase_offset = pi/4; end
    
    data_hat = pskdemod(z, M, phase_offset, 'gray');
    dec = dec2bin(data_hat, log2(M));
    Bitsre = zeros(1, log2(M));
    for n = 1:length(dec)
        Bitsre(n) = str2double(dec(n));
    end
end
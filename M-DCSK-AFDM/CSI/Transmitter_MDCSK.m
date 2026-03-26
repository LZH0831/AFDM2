function [Bits,Symbols0]=Transmitter_MDCSK(N_DFnT, M, N, C, c1, c2)
    P = N + C;                              
    Block_Num = 2; % 1个参考块，1个信息块
    bits_per_symbol = log2(M);
    
    %% 1. 生成比特并进行 M-PSK 映射
    Bits = randi([0 1], [1, bits_per_symbol]);
    dec_data = bin2dec(num2str(Bits));
    
    phase_offset = 0; 
    if M == 4
        phase_offset = pi/4;
    end
    sym = pskmod(dec_data, M, phase_offset, 'gray');
    a_s = real(sym);
    b_s = imag(sym);

    %% 2. 生成混沌序列 (扩频因子 beta = N)
    beta = N;
    x = rand() * 2 - 1;
    cx = zeros(beta, 1);
    for k = 1:beta
        cx(k) = x;
        x = 1 - 2 * x^2; 
    end
    
    cx = cx - mean(cx);
    cx = cx / norm(cx) * sqrt(beta); 
    
    cx_analytic = hilbert(cx);
    cy = imag(cx_analytic);
    cy = cy / norm(cy) * sqrt(beta); 
    
    tx_ref = cx;
    tx_info = a_s * cx + b_s * cy;

    %% 3. 发射机的离散 DA 变换 (使用外部传入的 c1, c2)
    W_IDFT = zeros(N_DFnT);
    for a=1:N_DFnT
       for b=1:N_DFnT
           W_IDFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N_DFnT);
       end 
    end
    T_1 = zeros(N_DFnT);
    for a=1:N_DFnT
       for b=1:N_DFnT
           if a==b
           T_1(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c1);
           end
       end 
    end
    T_1 = T_1*exp(-1i*pi/4);
    T_1H = conj(T_1);

    T_2 = zeros(N_DFnT);
    for a=1:N_DFnT
        for b=1:N_DFnT
            if a==b
               T_2(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c2);
            end
        end 
    end
    T_2 = T_2*1/sqrt(N_DFnT);
    T_2H = conj(T_2);

    Q_IDFnT0 = T_2H*W_IDFT*T_1H;
     
    %% IDFnT 与 循环前缀
    Symbols=zeros(N,1,Block_Num);
    Symbols(:,:,1) = Q_IDFnT0 * tx_ref;
    Symbols(:,:,2) = Q_IDFnT0 * tx_info;
    
    S=eye(N);
    T=[S(2*N-P+1:N,:);S];
    Symbols0=zeros(P,1,Block_Num);
    for a=1:Block_Num
        Symbols0(:,:,a)=T*Symbols(:,:,a);
    end
end
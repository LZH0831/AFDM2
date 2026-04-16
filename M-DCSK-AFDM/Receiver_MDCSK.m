function Bitsre=Receiver_MDCSK(N_DFnT,M,N,C,Equal,Symbols1,H0,SNR)
    P=N+C;
    Block_Num = 2;
    R=[zeros(N,P-N),eye(N)];
    Symbols2=zeros(N,1,Block_Num);
    
    %% 去除循环前缀
    for a=1:Block_Num
        Symbols2(:,:,a)=R*Symbols1(:,:,a);
    end
    
    %% 构造和应用 DFT 矩阵 
    c1=9/(2*N_DFnT);
    c2=1.414;
    W_IDFT=zeros(N_DFnT);
    for a=1:N_DFnT
       for b=1:N_DFnT
           W_IDFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N_DFnT);
       end 
    end
    T_1=zeros(N_DFnT);
    for a=1:N_DFnT
       for b=1:N_DFnT
           if a==b
           T_1(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c1);
           end
       end 
    end
    T_1=T_1*exp(-1i*pi/4);
    T_1H=conj(T_1);

    T_2=zeros(N_DFnT);
    for a=1:N_DFnT
        for b=1:N_DFnT
            if a==b
               T_2(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c2);
            end
        end 
    end
    T_2 =T_2*1/sqrt(N_DFnT);
    T_2H=conj(T_2);

    Q_IDFnT0=T_2H*W_IDFT*T_1H;
    P_DFnT0=conj(Q_IDFnT0);
      
    Symbols3=zeros(size(Symbols2));
    for count=1:Block_Num
        sample=P_DFnT0*Symbols2(:,:,count); 
        Symbols3(:,:,count)=sample;
    end

    %% 构造 MMSE 均衡矩阵
    S=eye(N);
    T=[S(2*N-P+1:N,:);S];
    H=R*H0*T;
    D=P_DFnT0*H*Q_IDFnT0; 

    if Equal==2
       G=pinv((D)'*D + eye(size(D)).*(1/SNR))*(D)';
    else
       G=inv(D); % ZF
    end
    
    %% 均衡
    Symbols4=zeros(size(Symbols3));
    for count=1:Block_Num
        Symbols4(:,:,count)=G*Symbols3(:,:,count);
    end

    %% M-DCSK 相关解调
    rx_ref  = real(Symbols4(:,:,1));
    rx_info = real(Symbols4(:,:,2));
    
    rx_ref_analytic = hilbert(rx_ref);
    rx_ref_y = imag(rx_ref_analytic);
    
    z_a = sum(rx_info .* rx_ref);
    z_b = sum(rx_info .* rx_ref_y);
    z = z_a + 1i * z_b; 
    
    phase_offset = 0; 
    if M == 4
        phase_offset = pi/4;
    end
    
    % 星座点判决
    data_hat = pskdemod(z, M, phase_offset, 'gray');
    
    % 转换为比特串
    dec = dec2bin(data_hat, log2(M));
    Bitsre = zeros(1, log2(M));
    for n=1:length(dec)
        Bitsre(n) = str2double(dec(n));
    end
end
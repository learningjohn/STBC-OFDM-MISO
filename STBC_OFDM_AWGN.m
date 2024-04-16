%AWGN下的STBC-OFDM系统
%%ZUST,2024/4/13
clc;clear;
N = 128;    %FFT长度，OFDM信号点数
Cp_len = 16; %循环前缀长度
Sym_num = 2000;%OFDM符号数
M = 16; %QAM调制阶数
SNR_dB = 0:1:20; %信噪比
NT = 2;NR = 1;%发送天线2，接收天线1

for i = 1:length(SNR_dB)
   biterr_count = 0;
   symerr_count = 0;
    for j=1:Sym_num/2 %一次传两个OFDM符号，传两回
       %%
       %发射端
       sym = randi([0,M-1],N*2,1);%符号信源
       sym_QAM  = qammod(sym,M,'UnitAveragePower',true);%QAM调制，功率归一化；
       sym_QAM1(:,1) = sym_QAM(1:N,:);
       sym_QAM2(:,1) = sym_QAM(N+1:end,:);%分为两组数据，并进行STBC编码
       STBC_code = cell(2);%元胞初始化编码矩阵

       %STBC编码
       STBC_code{1,1} = sym_QAM1 ; STBC_code{1,2} = sym_QAM2;
       STBC_code{2,1} = -conj(sym_QAM2);STBC_code{2,2} = conj(sym_QAM1);

       %IFFT,OFDM调制
       ofdm_mod = @(x) ifft(x,N).*sqrt(N);    %匿名函数,调制函数OFDM并归一化
       STBC_OFDM = cellfun(ofdm_mod,STBC_code,'UniformOutput',false);%ofdm调制

       %添加循环前缀
       cp_add = @(x) [x(end-Cp_len+1:end,:);x];%
       STBC_OFDM_addCP =  cellfun(cp_add,STBC_OFDM,'UniformOutput',false);
       %天线一的两个OFDM符号
       sym_OFDM_NT1 = cell2mat(STBC_OFDM_addCP(:,1)');
       %天线二的两个OFDM符号流
       sym_OFDM_NT2 = cell2mat(STBC_OFDM_addCP(:,2)');
       %%
       %信道
       %传输分两个时隙 
       h1 = 1;h2 = 1;  %信道响应，无衰落情况
       %天线一二同时发送经历不同的信道响应，但AWGN相同
       R1 = awgn(conv(h1,sym_OFDM_NT1(:,1))+conv(h2,sym_OFDM_NT2(:,1)),SNR_dB(i),'measured');%接收天线接收到的第一个信号
       R2 = awgn(conv(h1,sym_OFDM_NT1(:,2))+conv(h2,sym_OFDM_NT2(:,2)),SNR_dB(i),'measured');
       %R1 = sym_OFDM_NT1(:,1)+sym_OFDM_NT2(:,1);
       %R2 = sym_OFDM_NT1(:,2)+sym_OFDM_NT2(:,2);
       %%
       %接收端
       Y1 = fft(R1(Cp_len+1:end,:),N)./sqrt(N);%去循环前缀并FFT，OFDM解调
       Y2 = fft(R2(Cp_len+1:end,:),N)./sqrt(N);
       Y = [Y1;conj(Y2)];%第二个接收信号取共轭

       H1 = fft(h1,N);%计算信道响应频域，时域卷积频域点乘；
       H2 = fft(h2,N);
       D1 = diag(H1);D2 = diag(H2);
       D = [D1,D2;conj(D2),-conj(D1)];%构造频率响应矩阵

       %线性变换
       X_temp = D'*Y;
       %计算D1*D1+D2*D2*
       D_temp = diag(abs(H1).^2+abs(H2.^2));
       %克罗内尔积
       D_h = kron(eye(2),D_temp);
       %均衡
       X = inv(D_h)*X_temp;
       %解调
       sym_rec = qamdemod(X,M,'UnitAveragePower',true);
       %统计误符号数、误比特数
       symerr_count = symerr_count+sum(sym_rec~=sym);
       biterr_count = biterr_count+biterr(sym_rec,sym);
       
    end
    symerr_rate(i) = symerr_count/Sym_num/N;
    biterr_rate(i) = biterr_count/Sym_num/N/log2(M);
end
EbN0 = SNR_dB-10*log10(log2(M));
[ber_thoery,ser_thoery] = berawgn(EbN0,"qam",M);
figure()
semilogy(SNR_dB,symerr_rate);hold on
semilogy(SNR_dB,ser_thoery);
xlabel("信噪比（dB）");ylabel("误符号率");
title("STBC—OFDM,2发1收误符号率")
legend('STBC-OFDM','AWGN理论')
figure()
semilogy(SNR_dB,biterr_rate);
xlabel("信噪比（dB）");ylabel("误比特率");
title("STBC—OFDM,2发1收误比特率")
clear;
clc

%%% Performance comparison under different values of beta

Nr=4;
Nt=4;

ch_num=10^4;
SNR=(0:1:10)*3-10;

Cardinality=16;
sym_num=10^4; %% symbol number
%%%%%


SER_Huff_temp_beta3=zeros(1,length(SNR));
SER_Huff_temp_beta2=zeros(1,length(SNR));
SER_Huff_temp_beta1=zeros(1,length(SNR));
SER_Huff_temp_beta0=zeros(1,length(SNR));

for ich=1:1:ch_num  
    ich
    H=randn(Nr,Nt)+j*randn(Nr,Nt);
    
    source_int=randi([1 Cardinality], sym_num,1);
    noise=1/sqrt(2)*(randn(Nr,sym_num)+j*randn(Nr,sym_num));
    
    
for ix=1:1:length(SNR)
    snr=SNR(ix);
    gain=10^(snr/10);
%%%%% Codebook Generation of Different Modes
    %% Beta3
    [Huff_vec_beta3 Huff_codebook_beta3 Huff_metric_beta3]=ProVecSearchBeta(H,snr,3);
    %% Beta2
    [Huff_vec_beta2 Huff_codebook_beta2 Huff_metric_beta2]=ProVecSearchBeta(H,snr,2); 
    %% Beta1
    [Huff_vec_beta1 Huff_codebook_beta1 Huff_metric_beta1]=ProVecSearchBeta(H,snr,1); 
    %% Beta0
    [Huff_vec_beta0 Huff_codebook_beta0 Huff_metric_beta0]=ProVecSearchBeta(H,snr,0); 
    
%%%%%  ML Detection -- Find the codeword with minimum distance   %%%%%%%%%
    symbol_Huff_beta3 =Huff_codebook_beta3(:,source_int);
    symbol_Huff_beta2 =Huff_codebook_beta2(:,source_int);
    symbol_Huff_beta1 =Huff_codebook_beta1(:,source_int);
    symbol_Huff_beta0 =Huff_codebook_beta0(:,source_int);
%%%%% 
    R_Huff_beta3=sqrt(gain)*H*symbol_Huff_beta3+noise;   
    R_Huff_beta2=sqrt(gain)*H*symbol_Huff_beta2+noise;   
    R_Huff_beta1=sqrt(gain)*H*symbol_Huff_beta1+noise;   
    R_Huff_beta0=sqrt(gain)*H*symbol_Huff_beta0+noise;   
%%%%% 
        for ii=1:1:Cardinality
              temp_Huff_beta3=R_Huff_beta3-sqrt(gain)*H*Huff_codebook_beta3(:,ii)*ones(1,sym_num);
              temp_Huff_beta2=R_Huff_beta2-sqrt(gain)*H*Huff_codebook_beta2(:,ii)*ones(1,sym_num);
              temp_Huff_beta1=R_Huff_beta1-sqrt(gain)*H*Huff_codebook_beta1(:,ii)*ones(1,sym_num);
              temp_Huff_beta0=R_Huff_beta0-sqrt(gain)*H*Huff_codebook_beta0(:,ii)*ones(1,sym_num);
             if Nr>1
                norm_Huff_beta3(ii,:)=sum(conj(temp_Huff_beta3).*temp_Huff_beta3);   
                norm_Huff_beta2(ii,:)=sum(conj(temp_Huff_beta2).*temp_Huff_beta2); 
                norm_Huff_beta1(ii,:)=sum(conj(temp_Huff_beta1).*temp_Huff_beta1);   
                norm_Huff_beta0(ii,:)=sum(conj(temp_Huff_beta0).*temp_Huff_beta0); 
             else
                norm_Huff_beta3(ii,:)=abs(temp_Huff_beta3);    
                norm_Huff_beta2(ii,:)=abs(temp_Huff_beta2);  
                norm_Huff_beta1(ii,:)=abs(temp_Huff_beta1);  
                norm_Huff_beta0(ii,:)=abs(temp_Huff_beta0);  
             end
        end
        [val_beta3 decode_int_Huff_beta3]=min(norm_Huff_beta3);  
        [val_beta2 decode_int_Huff_beta2]=min(norm_Huff_beta2);  
        [val_beta1 decode_int_Huff_beta1]=min(norm_Huff_beta1);  
        [val_beta0 decode_int_Huff_beta0]=min(norm_Huff_beta0);  
%%%%%  ML Detection -- Find the codeword with minimum distance   %%%%%%%%%
        
        SER_Huff_snr_beta3(ix)=(1-sum(decode_int_Huff_beta3==source_int')/sym_num);  
        SER_Huff_snr_beta2(ix)=(1-sum(decode_int_Huff_beta2==source_int')/sym_num);  
        SER_Huff_snr_beta1(ix)=(1-sum(decode_int_Huff_beta1==source_int')/sym_num);  
        SER_Huff_snr_beta0(ix)=(1-sum(decode_int_Huff_beta0==source_int')/sym_num);
        
end
    SER_Huff_temp_beta3=SER_Huff_temp_beta3+SER_Huff_snr_beta3; 
    SER_Huff_temp_beta2=SER_Huff_temp_beta2+SER_Huff_snr_beta2; 
    SER_Huff_temp_beta1=SER_Huff_temp_beta1+SER_Huff_snr_beta1; 
    SER_Huff_temp_beta0=SER_Huff_temp_beta0+SER_Huff_snr_beta0; 
    
    SER_Huff_temp_beta3/ich
	SER_Huff_temp_beta2/ich
	SER_Huff_temp_beta1/ich
	SER_Huff_temp_beta0/ich
end
SER_Huff_array_beta3=SER_Huff_temp_beta3/ch_num;
SER_Huff_array_beta2=SER_Huff_temp_beta2/ch_num;
SER_Huff_array_beta1=SER_Huff_temp_beta1/ch_num;
SER_Huff_array_beta0=SER_Huff_temp_beta0/ch_num;


figure(1)
semilogy( SNR, SER_Huff_array_beta3,'-bo', SNR, SER_Huff_array_beta2,'-k^', SNR, SER_Huff_array_beta1,'-r*', SNR, SER_Huff_array_beta0,':bo','linewidth', 1.2 , 'MarkerSize', 7)
xlabel('SNR/dB','FontName','Times New Roman','FontSize',10)
ylabel('Capacity (Bits/Sec/Hz)','FontName','Times New Roman','FontSize',10,'Rotation',90)
legend( '\beta=3','\beta=2','\beta=1','\beta=0')
axis  auto
set(gca,'FontName','Times New Roman','FontSize',10)  
grid on

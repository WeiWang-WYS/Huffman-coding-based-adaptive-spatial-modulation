clear;
clc

%%% Performance comparison with spatial modulation and best antenna selection

Nr=4; % Receive Antenna Number
Nt=4; % Tramsmit Antenna Number

ch_num=10^4;
SNR=(0:1:10)*3-10;

Cardinality=16;
sym_num=10^4; %% symbol number
%%%%%

SER_spatial_temp=zeros(1,length(SNR)); %% Spatial Modulation
SER_Huff_temp=zeros(1,length(SNR)); %% Huffman coding based adaptive spatial modulation
SER_Ant_temp=zeros(1,length(SNR)); %% Antenna Selection
 
for ich=1:1:ch_num  
    ich
    H=randn(Nr,Nt)+j*randn(Nr,Nt);

    source_int=randi([1 Cardinality], sym_num,1);
    noise=1/sqrt(2)*(randn(Nr,sym_num)+j*randn(Nr,sym_num));

for ix=1:1:length(SNR)
    snr=SNR(ix);
    gain=10^(snr/10);
%%%%% Codebook Generation of Different Modes
    %% Spatial Modulation
    P_spatial=ones(1,Nt)/Nt; %% Probability vector
    [Spatial_codebook Spatial_metric]=CodeBookGen(P_spatial, snr, H);
    %% Huffman SM
    beta = 3;
    [Huff_vec Huff_codebook Huff_metric]= ProVecSearchBeta(H,snr,beta);
 
    %% Antenna Selection
    IndMat=diag(ones(Nt,1)); 
    [val index]=max(diag(H'*H));
    P_ant=IndMat(index,:); %% Probability vector
    [Ant_codebook Ant_metric]=CodeBookGen(P_ant, snr, H);
    
%%%%%  ML Detection -- Find the codeword with the minimum distance   %%%%%%%%%
    symbol_spatial=Spatial_codebook(:,source_int);
    symbol_Huff=Huff_codebook(:,source_int);
    symbol_Ant=Ant_codebook(:,source_int);
%%%%% 
    R_spatial=sqrt(gain)*H*symbol_spatial+noise;
    R_Huff=sqrt(gain)*H*symbol_Huff+noise;   
    R_Ant=sqrt(gain)*H*symbol_Ant+noise;           
%%%%% 
        for ii=1:1:Cardinality
             temp_spatial=R_spatial-sqrt(gain)*H*Spatial_codebook(:,ii)*ones(1,sym_num);
             temp_Huff=R_Huff-sqrt(gain)*H*Huff_codebook(:,ii)*ones(1,sym_num);
             temp_Ant=R_Ant-sqrt(gain)*H*Ant_codebook(:,ii)*ones(1,sym_num);
             if Nr>1
                norm_spatial(ii,:)=sum(conj(temp_spatial).*temp_spatial);    
                norm_Huff(ii,:)=sum(conj(temp_Huff).*temp_Huff);   
                norm_Ant(ii,:)=sum(conj(temp_Ant).*temp_Ant);   
            else
                norm_spatial(ii,:)=abs(temp_spatial);     
                norm_Huff(ii,:)=abs(temp_Huff);    
                norm_Ant(ii,:)=abs(temp_Ant); 
             end
        end
        [val decode_int_spatial]=min(norm_spatial);
        [val decode_int_Huff]=min(norm_Huff);  
        [val decode_int_Ant]=min(norm_Ant);  
%%%%%  ML Detection -- Find the codeword with the minimum distance   %%%%%%%%%

        SER_spatial_snr(ix)=(1-sum(decode_int_spatial==source_int')/sym_num);  
        SER_Huff_snr(ix)=(1-sum(decode_int_Huff==source_int')/sym_num);  
        SER_Ant_snr(ix)=(1-sum(decode_int_Ant==source_int')/sym_num); 
    end
    SER_spatial_temp=SER_spatial_temp+SER_spatial_snr;  
    SER_Huff_temp=SER_Huff_temp+SER_Huff_snr; 
    SER_Ant_temp=SER_Ant_temp+SER_Ant_snr; 
    SER_spatial_temp/ich
    SER_Huff_temp/ich
    SER_Ant_temp/ich
end

SER_spatial_array=SER_spatial_temp/ch_num;
SER_Huff_array=SER_Huff_temp/ch_num;
SER_Ant_array=SER_Ant_temp/ch_num;

figure(1)
semilogy(SNR,SER_spatial_array,'r-s',SNR,SER_Huff_array,'g-o',SNR,SER_Ant_array,'b-*','linewidth',1,'MarkerSize',5)
legend('SM','Huffman','AntSel','NorthWest'  )
xlabel('SNR/dB','FontName','Times New Roman','FontSize',11)
ylabel('SER','FontName','Times New Roman','FontSize',11,'Rotation',90)
axis auto
grid on
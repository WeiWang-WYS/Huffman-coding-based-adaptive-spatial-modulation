function [codebook metric]=CodeBookGen(P, SNR, H) 

%%% Generate the transmit codebook that corresponds to the probability vector P

BPSK=[-1 1];
QPSK=1/sqrt(2)*[1+j 1-j -1+j -1-j];
E8PSK=[exp(j*pi/8) exp(j*pi*3/8) exp(j*pi*5/8) exp(j*pi*7/8) ...
    exp(j*pi*9/8) exp(j*pi*11/8) exp(j*pi*13/8) exp(j*pi*15/8)];
E16QAM=1/sqrt(10)*[3+3*j 3+1*j 3-1*j 3-3*j 1+3*j 1+1*j 1-1*j 1-3*j -1+3*j -1+1*j -1-1*j -1-3*j ...
    -3+3*j -3+1*j -3-1*j -3-3*j];


Nt=length(P);

gain=10^(SNR/10);

cardinality=16;
const_size=P*cardinality;
const=[];

codebook=zeros(Nt,cardinality);
 

index_A=1;
index_B=const_size(1);
for i=1:1:Nt
    switch P(i)
        case 1
            const=E16QAM;
        case 1/2 
            const=E8PSK;
        case 1/4
            const=QPSK;
        case 1/8
            const=BPSK;
        case 0
            const=[];
    end
    codebook(i,index_A:index_B)=const;
    if i~=Nt
        index_A=index_A+const_size(i);
        index_B=index_B+const_size(i+1);
    end
end

 
cnt=0;
SER=0;
code_diff_square=zeros(Nt,Nt);
for ix=1:1:length(codebook)
    for jx=1:1:length(codebook)
       if ix~=jx
           cnt=cnt+1;
           code_diff(:,cnt)=codebook(:,ix)-codebook(:,jx);
           code_diff_square=code_diff_square+code_diff(:,cnt)*code_diff(:,cnt)';
           SER=SER+exp(-gain/4*code_diff(:,cnt)'*H'*H*code_diff(:,cnt));
       end
    end
end
SER=SER/(2*length(codebook));
metric=SER;

    
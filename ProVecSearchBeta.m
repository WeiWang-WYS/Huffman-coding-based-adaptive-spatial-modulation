function  [pro_vec codebook metric]=ProVecSearchBeta(H,SNR,beta)

switch beta
    case 3
        load('./ProbVecSpace/Nt4Beta3.mat');  %% Beta = 3, 
    case 2
        load('./ProbVecSpace/Nt4Beta2.mat');  %% Beta = 2 
    case 1
     	load('./ProbVecSpace/Nt4Beta1.mat');  %% Beta = 1 
    case 0
        load('./ProbVecSpace/Nt4Beta0.mat');  %% Beta = 0
end
    [Nr Nt]=size(H);
    metric=inf;
    codebook=[];
    pro_vec=[];
    for ix=1:1:length(vec_codebook) 
        [codebook_temp metric_temp]=CodeBookGen(vec_codebook(:,ix), SNR, H);
        if abs(metric_temp)<metric
            codebook=codebook_temp;
            metric=metric_temp; 
            pro_vec=vec_codebook(:,ix); 
        end
    end

 
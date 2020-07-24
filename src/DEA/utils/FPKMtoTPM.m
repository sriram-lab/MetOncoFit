function [TPM] = FPKMtoTPM(FPKM)   

%{ 
This function converts a vector of FPKM values to TPM values in order
to reduce bias caused by the length of the transcript

Inputs: 
FPKM: a vector of RNAseq values in FPKM

Outputs: 
TPM: a vector of size [n,1] of RNAseq values in TPM

%}

    TPM = zeros(length(FPKM),1);
    
    % cycle through each entry in the FPKM vector 
    % and convert to TPM 
    for a = 1:length(FPKM)
        TPM(a,1) = FPKM(a,1)./sum(FPKM).*10.^6;
    end

end
function [TPM] = FPKMtoTPM(FPKM)   

    TPM = zeros(length(FPKM),1);
    for a = 1:length(FPKM)
        TPM(a,1) = FPKM(a,1)./sum(FPKM).*10.^6;
    end

end
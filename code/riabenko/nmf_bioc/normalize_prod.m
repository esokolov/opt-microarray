function [AOut,COut] = normalize_prod(AIn, CIn)
norm = prod(AIn) ^ (1 / length(AIn));
AOut = AIn / norm;
COut = CIn * norm;
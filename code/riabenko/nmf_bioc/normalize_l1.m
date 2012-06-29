function [AOut,COut] = normalize_l1(AIn, CIn)
norm = sum(abs(AIn));
AOut = AIn / norm;
COut = CIn * norm;
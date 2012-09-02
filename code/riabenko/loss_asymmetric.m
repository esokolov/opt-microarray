function fit = loss_asymmetric(R)
fit = R;
fit(fit<0) = sqrt(-fit(fit<0));

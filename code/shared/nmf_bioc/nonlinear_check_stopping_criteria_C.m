function isConverged = nonlinear_check_stopping_criteria_C(C, C_prev, eps)
    fprintf('%f\n', mean(abs(C - C_prev)));
    if (mean(abs(C - C_prev)) < eps)
        isConverged = 1;
    else
        isConverged = 0;
    end
end
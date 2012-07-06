function isConverged = nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps)
    if (norm(I - Q, 'fro') < eps || ...
            (prevQuality > 0 && abs(currQuality - prevQuality) / currQuality < eps))
        isConverged = 1;
    else
        isConverged = 0;
    end
end
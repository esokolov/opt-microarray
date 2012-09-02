function isConverged = nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps)
    if (norm(I - A * C, 'fro') < eps || ...
            (prevQuality > 0 && abs(currQuality - prevQuality) / currQuality < eps))
        isConverged = 1;
    else
        isConverged = 0;
    end
end
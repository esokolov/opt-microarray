function isConverged = nmf_stopping(I, A, X, currQuality, prevQuality, eps)
    if (norm(I - A * X, 'fro') < eps || ...
            (prevQuality > 0 && abs(currQuality - prevQuality) / currQuality < eps))
        isConverged = 1;
    else
        isConverged = 0;
    end
end
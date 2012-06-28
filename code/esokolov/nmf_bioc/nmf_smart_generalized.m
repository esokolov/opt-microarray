function [A C] = nmf_smart_generalized(I, G, maxIterCnt, eps, errFuncName)
    if (strcmp(errFuncName, 'KL'))
        costFunc = @(Y, Q) sum(sum(bsxfun(@times, Q, log(bsxfun(@rdivide, Q + eps, Y + eps))) + Y - Q));
        errFunc = @(y, q) log(bsxfun(@rdivide, y, q));
    elseif (strcmp(errFuncName, 'BE1')) % Bose-Einstein divergence with alpha = 1
        costFunc = @(Y, Q) sum(sum(bsxfun(@times, Y, log(bsxfun(@rdivide, 2 * Y, Y + Q))) + ...
            bsxfun(@times, Q, log(bsxfun(@rdivide, 2 * Q, Y + Q)))));
        errFunc = @(y, q) log(bsxfun(@rdivide, y + q, 2 * q));
    elseif (strcmp(errFuncName, 'RJS'))
        costFunc = @(Y, Q) sum(sum(2 * bsxfun(@times, Y, log(bsxfun(@rdivide, 2 * Y, Y + Q))) + Q - Y));
        errFunc = @(y, q) log(bsxfun(@rdivide, y - q, y + q));
    else
        error('Unknown error function');
    end

    [A C] = nmf_init_als(I, G, eps);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        LearningRate_A = repmat(0.1 ./ sum(C, 2)', [size(A, 1) 1]);
        LearningRate_C = repmat(0.1 ./ sum(A, 1)', [1 size(C, 2)]);
        
        T = A' * errFunc(I, A * C);
        T = exp(bsxfun(@times, LearningRate_C, T));
        C = bsxfun(@times, C, T);
        
        T = errFunc(I, A * C) * C';
        T = exp(bsxfun(@times, LearningRate_A, T));
        A = bsxfun(@times, A, T);
        
        [A C] = nmf_normalize(A, C);
        
        currQuality = costFunc(I, A * C);
        if (nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end
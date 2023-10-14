function reI = tvgdIntensity(S, T, tau, eta, params)
%TVGDINTENSITY Estimate intensity using a total variation regularizer,
%optimized using subgradient method
%% Gradient descent
% parameters
Dx = make_convn_mat([1 -1], size(S), 'same', 'same');
Dy = make_convn_mat([1 -1]', size(S), 'same', 'same');
funcParams = struct('s', S(:), 'T', T, 'lambda', params.lambda, 'Dx', Dx, 'Dy', Dy);

maxIters = params.maxIters;
stepThresh = params.stepThresh;
maxStep = params.maxStep;

x = mleIntensity(S, T); % initialize as MLE
% x = zeros(size(S));
x(~isfinite(x)) = 0;
x(x<=0) = 0.2;
x = x(:);

% iterations
if params.display
    figure;
end
ftable = zeros(1, maxIters);
fbest = inf;
for i = 1:maxIters
    if params.display
        imshow(reshape(x/tau/eta, size(S)));
    end
    g = L(x, funcParams, 1);
    a = 0.1 /  sqrt(i);
    step = -a * g';
    step(step < -maxStep) = -maxStep;
    step(step > maxStep) = maxStep;
    
    x1 = x + step;
    x1(x1<=eps) = eps;
    step = x1 - x;
    
    stepNorm = norm(step);
    stepRatio = stepNorm / norm(x);
    x = x1;
    
    f = L(x, funcParams, 0);
    ftable(i) = f;
    if f < fbest
        fbest = f;
    end
    if params.display
        fprintf('Step %d:\t fval: %g,\t gnorm: %g,\t step: %g,\t stepRatio: %g\n',...
            i, f, norm(g), stepNorm, stepRatio);
    end
    if stepNorm < stepThresh
        break;
    end
    
    if any(~isfinite(f))
        fprintf('NaN or inf occured! Need to debug...\n');
        break;
    end
    
end

reI = reshape(x/tau/eta, size(S));
if params.display
    imshow(reI);
    figure; plot(ftable);
end

end

%% Function definition
function f = L(x, params, order)
if order == 0
    f = sum(-params.s/params.T.*log(1-exp(-x))+(params.T-params.s)/params.T.*x...
        +params.lambda*abs(params.Dx*x)...
        +params.lambda*abs(params.Dy*x));
else
    f = (-params.s/params.T./(1-exp(-x))+1)'...
        +params.lambda*sign(params.Dx*x)'*params.Dx...
        +params.lambda*sign(params.Dy*x)'*params.Dx;
end
end


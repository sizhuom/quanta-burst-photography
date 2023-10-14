function reI = tvgdRegularize(I, params)
%TVGDREGULARIZE Total variation regularization, implemented using
%subgradient method
%% Gradient descent
% parameters
Dx = make_convn_mat([1 -1], size(I), 'same', 'same');
Dy = make_convn_mat([1 -1]', size(I), 'same', 'same');
funcParams = struct('I', I(:), 'lambda', params.lambda, 'Dx', Dx, 'Dy', Dy);

maxIters = defaultField(params, 'maxIters', 300);
stepThresh = defaultField(params, 'stepThresh', 0);
maxStep = defaultField(params, 'maxStep', Inf);
display = defaultField(params, 'display', true);
displayFig = defaultField(params, 'displayFig', false);

x = I; % initialize as original field
% x = zeros(size(I));
x = x(:);

% iterations
if displayFig
    figure;
end
ftable = zeros(1, maxIters);
fbest = inf;
for i = 1:maxIters
    if displayFig
        imagesc(reshape(x, size(I))); colorbar; drawnow;
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
    if display
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

reI = reshape(x, size(I));
if display
    figure; plot(ftable);
end

end

%% Function definition
function f = L(x, params, order)
if order == 0
    f = sum((x-params.I).^2 ...
        +params.lambda*abs(params.Dx*x)...
        +params.lambda*abs(params.Dy*x));
else
    f = 2*(x-params.I)'...
        +params.lambda*sign(params.Dx*x)'*params.Dx...
        +params.lambda*sign(params.Dy*x)'*params.Dx;
end
end


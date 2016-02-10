function [JTJ_fd, gradient_fd, JTJ, grad] = scetch(a)

x = [1;2;3];

F = @(b)(exp(-b*x) - exp(-1.3*x));
JAC_fd = computeJAC(F, a, length(x), 1E-6);
JTJ_fd = JAC_fd'*JAC_fd;
gradient_fd = JAC_fd'*F(a);

JAC = -x.*exp(-a*x);
JTJ = JAC'*JAC;
grad = JAC'*F(a);

end

function JAC = computeJAC(fun, x, dataLen, tolX)

JAC = zeros(dataLen, length(x));
dx = max(0.25*tolX*abs(x), 1E-10);

for i = 1:length(x)
    xForw = x; xBackw = x;
    xForw(i) = x(i) + dx(i);
    xBackw(i) = x(i) - dx(i);
    
    F_forw = feval(fun, xForw);
    F_backw = feval(fun, xBackw);
    
    result = (F_forw - F_backw) / (2*dx(i));
    
    infInd = ~isfinite(result);
    if any(infInd)
        result(infInd) = mean(result(~infInd));
    end
    
    JAC(:,i) = result;
end

end
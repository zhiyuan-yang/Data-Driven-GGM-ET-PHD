function [y] = generateUniformMeasStarConvex(X,beta_D,p_D,Rot,sizeObject,R)

% Number of measurements
m_k = poissrnd(beta_D);

% Detections
dets = rand();
detections = (dets < p_D) & (m_k > 0);

% Allocate memory
y = struct('p',repmat({[]},1,1));

% Check if there is detection
if detections
    % Generate measurements
    y.p = StarConvexUniformsampling(X(1:2),Rot,sizeObject,R,m_k);
else
    % Return empty set of measurements
    y.p = [];
end

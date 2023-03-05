function pseudoMeasurement = f_meas_pseudo_LGM(x, z, phi, GMModel)
numberOfSigmaPoints = size(x,2);
pseudoMeasurement = zeros(2, numberOfSigmaPoints);
std = -pi:pi/8:7/8*pi;

for j = 1: numberOfSigmaPoints
    m = x(1:2,j);
    len = x(5,j);
    wid = x(6,j);
    yaw = atan2(x(4,j),x(3,j));
    theta = phi(j) - yaw  ;   %% theta = phi_SC - atan2(p_y,SC, p_x,SC)
    if theta > pi
        theta = theta - 2*pi;
    elseif theta < -pi
        theta = theta + 2*pi;
    end
    idx_theta = length(find( theta>= std));
    R = [cos(yaw), -sin(yaw)
     sin(yaw), cos(yaw)];
    S = [2/len   0
          0     2/wid];
    GMM = GMModel{idx_theta};
    % GMM.ComponentProportion: 1*3 vector
    % GMM.mu: 3*2 matrix
    % GMM.Sigma: 2*2*3 matrix
    %% previous method
    hat_z = S*(R'*(z - m));
    pdf = mvnpdf(hat_z',GMM.mu,GMM.Sigma);
    pdf = GMM.ComponentProportion .* pdf'; %pdf: 3*1vector
    [~ , ci] = max(pdf);
    mu = GMM.mu(ci,:)';
    v =  randn(1,2) * chol(GMM.Sigma(:,:,ci));
    pseudoMeasurement(:,j) =  m + R'*(inv(S)*( 1.1*mu ));
%    pseudoMeasurement(j) = norm(v) + norm(mu - hat_z) + 2*v*(mu - hat_z);
     
      %% new method
%      p = [1:GMM.NumComponents;GMM.ComponentProportion];
%      idx_mu = drnd(p,1);
%      mu = GMM.mu(idx_mu,:)';
%      v =  randn(1,2) * chol(GMM.Sigma(:,:,idx_mu));
%      pseudoMeasurement(:,j) =  m + R'*(inv(S)*( mu + v'));
     
end
end

function out = drnd(p,n)
    a = cumsum(p(2,:));
    b = rand(n,1);
    out = zeros(1,n);
    for k=1:n
        c = find(a<b(k));
        if isempty(c)
            out(k) = p(1,1);
        else
            out(k)=p(1,c(end)+1);
        end
    end
end

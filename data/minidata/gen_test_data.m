%testdata
load('GMModel.mat');
rng('default');
GMM = GMModel{9};
ego_car_m = [0;0];
ego_car_v = [12;0];
obs_car_m = [6;0];
obs_car_v = [12;0];
T=1:1:20;
obs_car_tra = obs_car_m + obs_car_v * T;
ego_car_tra = ego_car_m + ego_car_v * T;
w = 2.3;
l = 4.2;

lambda = 2;   %mean number of radar points
s = RandStream('mlfg6331_64');
theta = 0;
R = [cos(theta),-sin(theta);
     sin(theta),cos(theta)];
S =[l/2 0;
    0 w/2];
testdata = cell(1,20);
for i=1:1:20
    figure(1);
    n = poissrnd(lambda);
    c = randsample(s,[1,2,3],n,true,GMM.ComponentProportion);
    points = zeros(2,length(c));
    for j = 1:1:length(c)
        mu = GMM.mu(c(j),:);
        v =  randn(1,2) * chol(GMM.Sigma(:,:,c(j)));
        points(:,j) = obs_car_tra(:,i) + R*S*(mu + v)';
    end
    testdata{i} = points;
    plot(points(1,:),points(2,:),'go');
    hold on;
    rectangle('Position',[obs_car_tra(:,i)'-[l/2 w/2],l,w],'EdgeColor','b','LineWidth',1.5);
    rectangle('Position',[ego_car_tra(:,i)'-[l/2 w/2],l,w],'EdgeColor','r',LineWidth=1.5);
    plot(obs_car_tra(1,i),obs_car_tra(2,i),'bx');
    plot(ego_car_tra(1,i),ego_car_tra(2,i),'rx');
    axis([-10,260,-2,2]);
    xticks(-10:10:260);
end
hold off;
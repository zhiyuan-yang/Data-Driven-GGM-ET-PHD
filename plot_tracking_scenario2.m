%plot the automotive tracking scenario
clear;
clc;
addpath(genpath('data/scenario3'));
addpath('data');
load('output_ERHM_phd.mat');
load('output_GGIW_phd.mat');
load('output_LGM_phd.mat');
load("scenario3.mat");
meas = scenario3.meas;
truth = scenario3.truth;
ego = scenario3.ego;
len = 4.7;
wid = 1.8;

%% plot the center of target
figure(1);
hold on;
plot([98,99.5],[41,-21.7],'-k','HandleVisibility','off');
plot([43,84.8,91.3,93.4,94.9,96.1,96.8,97.6,98.2,98.2],[-1.64,-2,-1.8,-1.8,-1.8,-2,-3.3,-4.8,-8.5,-20.11],'-k','HandleVisibility','off');
start = scatter([43,98],[-1.64,41],'ok');
ending = scatter([99.7,98.2],[-21.7 -20.11],'xk');
legend([start,ending],{'起始位置','终止位置'});
xlabel('X/m');
ylabel('Y/m');
hold off;


%% plot the truth and track without extend
figure(2);
hold on;
for k=1:length(meas)
    s1 = plot(meas{k}(1,:),meas{k}(2,:),'gx');
end
for k=1:length(truth)
    if ~isempty(output_LGM_phd.state_X_hat{k})
        s2= scatter(output_LGM_phd.state_X_hat{k}(1,:),output_LGM_phd.state_X_hat{k}(2,:),'k.');
    end
end
s3 = plot([98,99.5],[41,-21.7],'-k','LineWidth',1.5);
plot([43,84.8,91.3,93.4,94.9,96.1,96.8,97.6,98.2,98.2],[-1.64,-2,-1.8,-1.8,-1.8,-2,-3.3,-4.8,-8.5,-20.11],'-k','HandleVisibility','off','LineWidth',1.5);
hold off;
xlabel('X/m');
ylabel('Y/m');
legend([s1,s2,s3],{'量测值','估计值','真实值'});


%% Plot the OSPA-Distance
figure(3)
load('data\montecarlo_ospa_GGIW_s2.mat');
s1 = plot(1:55,ave_ospa,'b.','LineStyle','-');
hold on;
load('data\montecarlo_ospa_LGM_s2.mat');
s2 = plot(1:55,ave_ospa,'rx','LineStyle','--');
load('data\montecarlo_ospa_ERHM_s2.mat');
s3 = plot(1:55,ave_ospa,'ks','LineStyle','--');
xlabel('TimeStep');
ylabel('OSPA/m');
legend([s1,s2,s3],{'GGIW','LGM','ERHM'});
hold off;

%% Plot the truth and track with extend
figure(3);
hold on;
for k=10:25:60
    rectangle('Position',[ego{k}.position(1)-len/2,ego{k}.position(2)-wid/2, ...
        len,wid],'EdgeColor','r','LineWidth',1.5,'LineStyle','--');
    for i=1:size(truth{k}.position,2)
        rectangle('Position',[truth{k}.position(1,i)-len/2,truth{k}.position(2,i)-wid/2, ...
            len,wid],'EdgeColor','b','LineWidth',1.5,'LineStyle','-');
    end
    plot(meas{k}(1,:),meas{k}(2,:),'g.');
    plot(output_LGM_phd.state_X_hat{k}(1,:),output_LGM_phd.state_X_hat{k}(2,:),'mx');
    plot(output_GGIW_phd.state_X_hat{k}(1,:),output_GGIW_phd.state_X_hat{k}(2,:),'cx');
    plot(output_ERHM_phd.state_X_hat{k}(1,:),output_ERHM_phd.state_X_hat{k}(2,:),'kx');
    %xlim([0 350]);
    ylim([-30 30]);


    
end
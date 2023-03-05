%plot the automotive tracking scenario
clear;
clc;
addpath(genpath('data/scenario1'));
addpath('data');
load('output_ERHM_phd.mat');
load('output_GGIW_phd.mat');
load('output_LGM_phd_s1.mat');
load("scenario1.mat");
meas = scenario1.meas;
truth = scenario1.truth;
ego = scenario1.ego;
len = 4.7;
wid = 1.8;

%% plot the center of target
figure(1);
hold on;
plot([59.35,351.35],[0.1044,0.1667],'-k','HandleVisibility','off');
plot([42.35,348.35],[4,4],'-k','HandleVisibility','off');
plot([33.35,144.35],[-4,-4],'-k','HandleVisibility','off');
start = scatter([59.35 42.35 33.35],[0 4 -4],'ok');
ending = scatter([351.35,348.35,144.35],[0.1667 4 -4],'xk');
legend([start,ending],{'起始位置','终止位置'});
xlabel('X/m');
ylabel('Y/m');
xlim([0 500]);
ylim([-10 10]);
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
s3 = plot([59.35,351.35],[0.1044,0.1667],'-k','HandleVisibility','off');
plot([42.35,348.35],[4,4],'-k','HandleVisibility','off');
plot([33.35,144.35],[-4,-4],'-k','HandleVisibility','off');
hold off;
xlabel('X/m');
ylabel('Y/m');
legend([s1,s2,s3],{'量测值','估计值','真实值'});


%% Plot the OSPA-Distance
figure(3)
load('data\montecarlo_ospa_GGIW_s1.mat');
s1 = plot(1:200,ave_ospa,'b.','LineStyle','-');
hold on;
load('data\montecarlo_ospa_LGM_s1.mat');
s2 = plot(1:200,ave_ospa,'rx','LineStyle','--');
load('data\montecarlo_ospa_ERHM_s1.mat');
s3 = plot(1:200,ave_ospa,'ks','LineStyle','--');
xlabel('TimeStep');
ylabel('OSPA/m');
legend([s1,s2,s3],{'GGIW','LGM','ERHM'});
hold off;



%% Plot the truth and track with extend
figure(4)
hold on;
load("data/scenario1/scenario1.mat");
truth = scenario1.truth;
meas = scenario1.meas;
len = 4.7;
wid = 2;
%% plot the true extend and meas
for k=51:5:71
        rectangle('Position',[truth{k}.position(1,2)-len/2,truth{k}.position(2,2)-wid/2, ...
            len,wid],'EdgeColor','b','LineWidth',1.5,'LineStyle','-');
    s1 = scatter(truth{k}.position(1,2),truth{k}.position(2,2),'bx');
    for i=1:size(meas{k},2)
        if meas{k}(2,i)>2
            scatter(meas{k}(1,i),meas{k}(2,i),'g.');   
        end
    end
end
%% plot GGIW
load('data\monte_carlo_GGIW_s1.mat');
GGIW = outcome{1}.output_GGIW_phd;
s2 = scatter(GGIW.state_X_hat{51}(1,1),GGIW.state_X_hat{51}(2,1),'rx');
L = GGIW.state_E_hat{51};
a = L(1,1,1);
b = L(1,2,1);
c = L(2,1,1);
d = L(2,2,1);
if a<d
    f = @(x,y) a*(x-GGIW.state_X_hat{51}(1,1))^2 + d*(y-GGIW.state_X_hat{51}(2,1))^2-1;
else
    f = @(x,y) d*(x-GGIW.state_X_hat{51}(1,1))^2 + a*(y-GGIW.state_X_hat{51}(2,1))^2-1;
end
fimplicit(f,[GGIW.state_X_hat{51}(1,1)-10, GGIW.state_X_hat{51}(1,1)+10,...
    GGIW.state_X_hat{51}(2,1)-10, GGIW.state_X_hat{51}(2,1)+10]);
for k=56:5:71
    s2 = scatter(GGIW.state_X_hat{k}(1,2),GGIW.state_X_hat{k}(2,2),'rx');
    L = GGIW.state_E_hat{k};
        a = L(1,1,2);
        b = L(1,2,2);
        c = L(2,1,2);
        d = L(2,2,2);
        if a<d
            f = @(x,y) a*(x-GGIW.state_X_hat{k}(1,2))^2 + d*(y-GGIW.state_X_hat{k}(2,2))^2-1;
        else
            f = @(x,y) d*(x-GGIW.state_X_hat{k}(1,2))^2 + a*(y-GGIW.state_X_hat{k}(2,2))^2-1;
        end
        fimplicit(f,[GGIW.state_X_hat{k}(1,2)-10, GGIW.state_X_hat{k}(1,2)+10,...
            GGIW.state_X_hat{k}(2,2)-10, GGIW.state_X_hat{k}(2,2)+10],'Color','r');
end
%% plot LGM
load('data\monte_carlo_LGM_s1.mat');
LGM = outcome{1}.output_LGM_phd;
for k=51:5:56
    s3 = scatter(LGM.state_X_hat{k}(1,1),LGM.state_X_hat{k}(2,1),'kx');
    L = LGM.state_E_hat{k};   
    rectangle('Position',[LGM.state_X_hat{k}(1,1)-1/2*LGM.state_E_hat{k}(1,1),LGM.state_X_hat{k}(2,1)-1/2*LGM.state_E_hat{k}(2,1), ...
            LGM.state_E_hat{k}(1,2),LGM.state_E_hat{k}(2,2)],'EdgeColor','k','LineWidth',1.5,'LineStyle','--');
end
for k=61:5:66
    s3 = scatter(LGM.state_X_hat{k}(1,2),LGM.state_X_hat{k}(2,2),'kx');
    L = LGM.state_E_hat{k};   
    rectangle('Position',[LGM.state_X_hat{k}(1,2)-1/2*LGM.state_E_hat{k}(1,2),LGM.state_X_hat{k}(2,2)-1/2*LGM.state_E_hat{k}(2,2), ...
            LGM.state_E_hat{k}(1,2),LGM.state_E_hat{k}(2,2)],'EdgeColor','k','LineWidth',1.5,'LineStyle','--');
end
k = 71;
s3 = scatter(LGM.state_X_hat{k}(1,3),LGM.state_X_hat{k}(2,3),'kx');
    L = LGM.state_E_hat{k};   
    rectangle('Position',[LGM.state_X_hat{k}(1,3)-1/2*LGM.state_E_hat{k}(1,3),LGM.state_X_hat{k}(2,3)-1/2*LGM.state_E_hat{k}(2,3), ...
            LGM.state_E_hat{k}(1,3),LGM.state_E_hat{k}(2,3)],'EdgeColor','k','LineWidth',1.5,'LineStyle','--');
ylim([-15 15]);
legend([s1,s2,s3],{'真实质心位置','GGIW估计质心位置','LGM估计质心位置'})
xlabel('X/m');
ylabel('Y/m')
hold off;

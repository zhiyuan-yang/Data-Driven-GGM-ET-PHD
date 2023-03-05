clear
clc
load('sensordata.mat');
N = length(sensordata);
num_obs_car = 2;

%% setting
% ego is 1*N cell each cell contains a struct
%     ego.position: 3*1 vector
%     ego.velocity: 3*1 vector
%     ego.ges: 3*1 vector[roll  pitch yaw]
% truth is n*N cell each cell contains a struct
%     truth.position: 3*1 vector
%     truth.velocity: 3*1 vector
%     truth.ges: 3*1 vector
% ego and truth are both in global system
% meas is in ego rectangle system
% note that the center of car in Matlab Driving toolbox is on the back rear

measurement = cell(1,N);
ego = cell(1,N);
truth = cell(1,N);

for i=1:1:N
    if i>=29 && i<52
        [ego{i}.position, truth{i}.position(1,:),truth{i}.position(2,:)...
                             ]=  sensordata(i).ActorPoses.Position;    
        [ego{i}.ges(1,1), truth{i}.ges(1,1),truth{i}.ges(1,2)]...
        = sensordata(i).ActorPoses.Roll;

        [ego{i}.ges(2,1), truth{i}.ges(2,1),truth{i}.ges(2,2)]...
        = sensordata(i).ActorPoses.Pitch;

        [ego{i}.ges(3,1), truth{i}.ges(3,1),truth{i}.ges(3,2)]...
        = sensordata(i).ActorPoses.Yaw;
    else
        [ego{i}.position,~, truth{i}.position(1,:)]= ...
                                                sensordata(i).ActorPoses.Position;
        [ego{i}.ges(1,1), ~,truth{i}.ges(1,1)]...
        = sensordata(i).ActorPoses.Roll;

        [ego{i}.ges(2,1), ~,truth{i}.ges(2,1)]...
        = sensordata(i).ActorPoses.Pitch;

        [ego{i}.ges(3,1), ~,truth{i}.ges(3,1)]...
        = sensordata(i).ActorPoses.Yaw;
    end
    ego{i}.position = ego{i}.position'+[1.35; 0; 0];
    truth{i}.position = truth{i}.position'+ [1.35;0;0];

    measurement{i} = zeros(2,length(sensordata(i).ObjectDetections));
    for k=1:1:length(sensordata(i).ObjectDetections)
        measurement{i}(:,k) = sensordata(i).ObjectDetections{k}.Measurement(1:2); %2*1 vector
        measurement{i}(1:2,k) = angle2rotmatrix(ego{i}.ges) * measurement{i}(1:2,k)+...
                    ego{i}.position(1:2,:) - [1.35;0];
    end
end
scenario3.ego = ego;
scenario3.truth = truth;
scenario3.meas = measurement;
save('scenario3.mat','scenario3');

function rotmatrix = angle2rotmatrix(ang)
    % this function only considers a 2-D scenario
    rotmatrix = [cos(ang(3)),-sin(ang(3));sin(ang(3)),cos(ang(3))];
end
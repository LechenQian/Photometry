function  all_RotEnc = get_aligned_RotEnc (ThisType,Data,posTrials,time,alignTo,protocol)

plotauxfigs = 0;

dbstop if error

switch alignTo
    case 'CS'
        alignstate = 'Foreperiod';
    case 'US'
        alignstate = 'Trace';
end

distance = nan(length(posTrials), length(time));
velocity = distance;
acceleration = distance;
ZeroEncoderTimes = [];
for i = 1:length(posTrials)   % Loop across trials of this type
    
    aligning_time = Data.States(posTrials(i)).(alignstate)(2) + 0.001; % align to the sample 1ms after the end of the foreperiod/trace period (depending on CS or US alignment)        
    
    RE_angle = Data.RotEncoder(posTrials(i)).Positions;
    RE_time = Data.RotEncoder(posTrials(i)).Times;
    
    drvposicoes = [0; diff(RE_angle')];  % To find discontinuities
    ind = find(abs(drvposicoes) > 5);   % positions of discontinuities
    
    % This loop is to eliminate rare single/double samples that are out of place
    % due to some error in the rotary encoder acquisition, creating 2
    % consecutive discontinuities that are not real.
    k = 1;
    aux = length(ind);
    while k < aux
        if (ind(k+1) - ind(k) <= 2) && RE_angle(ind(k)) ~= 180 && RE_angle(ind(k)) ~= -180  % Because there can be some jumping at the edges (+-180deg), these disontinuities must be corrected instead of ignored
            ind = [ind(1:k-1); ind(k+1:end)];
            RE_angle(ind(k):ind(k)+1) = mean([RE_angle(ind(k)-1) RE_angle(ind(k)+2)]);
        end
        k = k+1;
        aux = length(ind);
    end
    
    if ~isempty(ZeroEncoderTimes)
        [~, RE_zero] = min(abs(RE_time(ind)-ZeroEncoderTimes));
    else
        RE_zero = length(ind)+1;
    end
    
    RE_angle_corr = RE_angle;    % Arrange the angle vector to get rid of the discontinuities (when forced to zero or wen reaching + or - 150 deg)
    for j = 1:length(ind)
        if j ==  RE_zero
            RE_angle_corr(1:ind(j)-1) = RE_angle_corr(1:ind(j)-1) - RE_angle_corr(ind(j)-1);
            RE_angle_corr(ind(j):end) = RE_angle_corr(ind(j):end) - RE_angle_corr(ind(j));
        else
            if RE_angle_corr(ind(j)) < RE_angle_corr(ind(j)-1)      % if the signal goes down at this dicontinuity
                RE_angle_corr(1:ind(j)-1) = RE_angle_corr(1:ind(j)-1);
                RE_angle_corr(ind(j):end) = RE_angle_corr(ind(j):end) + abs(RE_angle_corr(ind(j)) - RE_angle_corr(ind(j)-1));
            else                                                      % if the signal goes up at this dicontinuity
                RE_angle_corr(1:ind(j)-1) = RE_angle_corr(1:ind(j)-1);
                RE_angle_corr(ind(j):end) = RE_angle_corr(ind(j):end) - abs(RE_angle_corr(ind(j)) - RE_angle_corr(ind(j)-1));
            end
        end
        %         length(RE_angle_corr)
        %         figure; plot(RE_time,RE_angle_corr)
        %         pause
    end
    if strcmp(protocol, 'EBC1CS')
        RE_distance = -(RE_angle_corr/360)*2*pi*0.0635;  % 0.0635 m is the radius of the treadmill (to calculate distance travelled)
    elseif strcmp(protocol, 'EBC2CS')
        RE_distance = (RE_angle_corr/360)*2*pi*0.102;  % 0.102 m is the radius of the treadmill (to calculate distance travelled)
    elseif strcmp(protocol(end-3:end), 'MSFP') || strcmp(protocol,'Es163O') || strcmp(protocol,'Eshl16')
        RE_distance = -(RE_angle_corr/360)*2*pi*0.102;  % 0.102 m is the radius of the treadmill (to calculate distance travelled)
    end
    RE_time_aligned = RE_time - aligning_time;
    
    RE_velocity = diff(RE_distance)./diff(RE_time);  % m/s
    
    if length(RE_distance) > 1
        RE_velocity = [RE_velocity(1) RE_velocity];
        RE_velocity(isinf(abs(RE_velocity))) = NaN;
        RE_acceleration = diff(RE_velocity)./diff(RE_time);  % m/s2
        RE_acceleration = [RE_acceleration(1) RE_acceleration];
        RE_acceleration(isinf(abs(RE_acceleration))) = NaN;
    else
        RE_velocity = nan(size(RE_distance));
        RE_acceleration = nan(size(RE_distance));
    end
    
    % Resample to have the same samples as the rest of the data:
    pos1 = find(time-RE_time_aligned(1) >= 0, 1,'first');
    pos2 = find(time-RE_time_aligned(end) <= 0, 1,'last');
    
    if pos2 > pos1
        serietempo_dist = timeseries(RE_distance,RE_time_aligned,'Name','Distance');
        serietempo_vel = timeseries(RE_velocity,RE_time_aligned,'Name','Velocity');
        serietempo_acc = timeseries(RE_acceleration,RE_time_aligned,'Name','Acceleration');
        res_serietempo_dist = resample(serietempo_dist,time(pos1:pos2));
        res_serietempo_vel = resample(serietempo_vel,time(pos1:pos2));
        res_serietempo_acc = resample(serietempo_acc,time(pos1:pos2));

        distance(i,pos1:pos2) = squeeze(res_serietempo_dist.data);
        velocity(i,pos1:pos2) = squeeze(res_serietempo_vel.data);
        acceleration(i,pos1:pos2) = squeeze(res_serietempo_acc.data);
    end
    
    if plotauxfigs 
        figure; hold on
        subplot(2,2,1); plot(RE_time_aligned,RE_angle_corr); title('Angle'); ylabel('Deg'); xlabel('Time (s)')
        subplot(2,2,2); plot(RE_time_aligned,RE_distance, 'linewidth', 2); hold on; plot(time(pos1:pos2),distance(i,pos1:pos2)); title('Distance'); ylabel('m'); xlabel('Time (s)')
        subplot(2,2,3); plot(RE_time_aligned,RE_velocity, 'linewidth', 2); hold on; plot(time(pos1:pos2),velocity(i,pos1:pos2)); title('Velocity'); ylabel('m/s'); xlabel('Time (s)')
        subplot(2,2,4); plot(RE_time_aligned,RE_acceleration, 'linewidth', 2); hold on; plot(time(pos1:pos2),acceleration(i,pos1:pos2)); title('Acceleration'); ylabel('m/(s2)'); xlabel('Time (s)')
    end
    
all_RotEnc = cat(3,distance', velocity', acceleration');
end
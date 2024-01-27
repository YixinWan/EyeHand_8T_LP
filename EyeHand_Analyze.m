%% for the latest version of eye-hand experiment
% hand file: 1 excel
%   parameters: target_x  target_y	baseline rotation aftereffect  RT
%                movement_frames	MoveOnsetTime1	MoveEndTime1 HandTrialStart
% eye file: annotations  gaze_positions_on_surface_screen
%  parameters: 'timestamp','label'
%              'gaze_timestamp','x_norm','y_norm','on_surf','confidence' 

%% choose which part to run
hand_prepro = 1;
eye_prepro = 1;
plot_hand = 1;
plot_eye = 1;


%% upload hand file
% Define the path to the hand data file
hand_location = 'D:/PSYCHOPY/eyehand_baseline/data';
hand_export = fullfile(hand_location, '001_new_8T_2024-01-27_10h04.34.736.csv');
% Read the CSV file into a table
hand_table = readtable(hand_export);
% Keep the specified columns
hand_table = hand_table(:, {'target_x','target_y','baseline','rotation','aftereffect','RT','movement_frames',...
    'MoveOnsetTime1','MoveEndTime1','HandTrialStart',...
    'baseline_trials_thisN','rotation_trials_thisN','aftereffect_trials_thisN','group'});
% Remove the first row, drop NaN values in 'target_x', and reset the index
hand_table = hand_table(2:end, :);
hand_table = rmmissing(hand_table, 'DataVariables', 'target_x');
%% upload eye files
% Define the path to the eye data file
eye_location = 'D:/pupils/recordings/2024_01_27/001/exports/000';
eye_export_annotation = fullfile(eye_location, 'annotations.csv');
eye_export_gaze = fullfile(eye_location, 'surfaces', 'gaze_positions_on_surface_screen.csv');
% Read the CSV file into a table
eye_annotation = readtable(eye_export_annotation);
eye_gaze = readtable(eye_export_gaze);
% Keep the specified columns
eye_annotation = eye_annotation(:, {'timestamp','label'});
eye_gaze = eye_gaze(:, {'gaze_timestamp','x_norm','y_norm','on_surf','confidence'});
% Remove the row that on_surf=FALSE and confidence < 0.9
% Define the conditions
% % Convert 'on_surf' column to logical
% logical_on_surf = strcmp(eye_gaze.on_surf, 'TRUE');
% Define the conditions
% condition_on_surf = eye_gaze.confidence >= 0.6; %logical_on_surf & (eye_gaze.confidence >= 0.9);
% eye_gaze = eye_gaze(condition_on_surf, :);


%% Define scales for data
    % Cut hand data from 1.5 second before movement onset to 2.5 seconds after
    % movmement onset. This should be enough to capture both the origin hold
    % and the feedback (I hope).
    tMin = -1.5;
    tMax = 2.5;
    dt = 0.01; % We'll resample everything to 10 ms sampling
    tVals = (tMin:dt:tMax)';
    numTVals = length(tVals);
    numCoord = 2;
    numTrials = nansum(hand_table.baseline) + nansum(hand_table.rotation) + nansum(hand_table.aftereffect); %80;
    
    % Units are image width/height in a normalized coordinate system 
    % with 0,0 origin in the bottom left and 1,1 at top right (pupil labs docs)
    % screen is 53*30 cm 
    % https://docs.pupil-labs.com/core/terminology/#coordinate-system
    screenWidth = 61; %35.9; %53;
    screenHeight = 34.3; %26.7; %30;
    
    % The eye starts from the bottom left of the screen and has normalized
    % coordinates over the screen. 
    x0Eye = 0.5; %原点坐标
    y0Eye = 0.5;
    xUnitEye = screenWidth;
    yUnitEye = screenHeight;
    
    % Determine a confidence threshold for the gaze data
    confidenceThreshold = 0.9;  % 0.97


%% eye
if eye_prepro
    %% 对齐手眼的时间
    timeGazeStart = NaN*zeros(numTrials, 1);
    timeGazeEnd = NaN*zeros(numTrials, 1);
    gaze = NaN*zeros(numTVals, numTrials, numCoord);
    
    trialLabel = eye_annotation.label;  %annotations表里的'label'项
    timestamp = eye_annotation.timestamp;  %annotations表里的'timestamp'项
    trialStartLabels = find(startsWith(trialLabel, "start_trial"));  
    trialEndLabels = find(startsWith(trialLabel, "end_trial"));   
    movementStartLabels = find(startsWith(trialLabel, "move_start"));
    movementEndLabels = find(startsWith(trialLabel, "move_end"));    
    numTrialStartLabels = length(trialStartLabels);
    
    expStart = timestamp(startsWith(trialLabel, "start_experiment")); 
    trialStarts = zeros(numTrialStartLabels,1);  %定义每次trial开始时的timestamp
    trialEnds = zeros(numTrialStartLabels,1);    
    movementStarts = zeros(numTrialStartLabels,1); %定义每次eye开始动时的timestamp
    movementEnds = zeros(numTrialStartLabels,1); 
    
    for trialNum = 1:numTrialStartLabels            
        trialStarts(trialNum) = timestamp(trialStartLabels(trialNum)); %每次trial开始时的timestamp
        trialEnds(trialNum) = timestamp(trialEndLabels(trialNum)); %每次trial结束时的timestamp
        
        movementStarts(trialNum) = timestamp(movementStartLabels(trialNum));
        movementEnds(trialNum) = timestamp(movementEndLabels(trialNum));
    end
    
    timeGazeExpStart = expStart;  %实验开始的时间 %from anotations file
    timeGazeStart = trialStarts;  %每次trial开始时的时间 %from anotations file
    timeGazeEnd = trialEnds;      %每次trial结束时的时间 %from anotations file
    timeMoveStart = movementStarts;  %每次eye开始动时的timestamp %from anotations file
    timeMoveEnd = movementEnds;      %每次eye结束动时的timestamp %from anotations file
    handTrialTimeStarts = hand_table.HandTrialStart; %每次trial开始的时刻（in psychopy）
    handTrialTime0 = hand_table.MoveOnsetTime1; %手开始动的时刻（in psychopy）
    
    for trialNum = 1:numTrials
        thisStartTime = trialStarts(trialNum);
        thisEndTime = trialEnds(trialNum);
                    
        theseGazeIndexes = find( ...    %找到gaze_timestamp的时间大于实验开始的时间且小于试验结束的部分
                            eye_gaze.gaze_timestamp > thisStartTime & ...  %筛掉不是实验过程中的数据
                            eye_gaze.gaze_timestamp < thisEndTime & ...
                            eye_gaze.confidence > confidenceThreshold);
        theseGazeT = eye_gaze.gaze_timestamp(theseGazeIndexes); %真正的实验过程中获取gaze的时间点
    
        theseGazeX = eye_gaze.x_norm(theseGazeIndexes);
        theseGazeY = eye_gaze.y_norm(theseGazeIndexes);
    
        badSample = find(diff(theseGazeT) <= 0);
        while ~isempty(badSample)
            allSample = 1:length(theseGazeT);
            keepSample = setdiff(allSample, badSample+1);
            theseGazeT = theseGazeT(keepSample);
            theseGazeX = theseGazeX(keepSample);
            theseGazeY = theseGazeY(keepSample);
            badSample = find(diff(theseGazeT) <= 0);
        end
                        
        if length(theseGazeIndexes) > 1
        %theseGazeTZeroed 让手动和眼动的时间起点一致
            theseGazeTZeroed = (theseGazeT - thisStartTime) - (handTrialTime0(trialNum) - handTrialTimeStarts(trialNum));
            t = tVals;
            x = interp1(theseGazeTZeroed, theseGazeX, tVals, 'makima',nan); %makima
            y = interp1(theseGazeTZeroed, theseGazeY, tVals, 'makima',nan);
                      
            gaze(:,trialNum,1) = (x-x0Eye)*xUnitEye;  %gaze x （屏幕上实际尺寸距离）
            gaze(:,trialNum,2) = (y-y0Eye)*yUnitEye;  %gaze y
        end
    end

    %% Zero-correct the eye origin for origin location
    
    originMin = -1.5; % -1.5
    originMax = -1.0; % -1.0
    method = 'rloess'; %用method参数指定平滑数据的方法
    span = 25;  %用span参数指定移动平均滤波器的窗宽,span为奇数，默认为5
    
    tOr = find(originMin < tVals & tVals < originMax); %盯着原点看的时间段
    g = gaze;
    % Get the median position of the eyes at origin fixation
    gOr = squeeze(median(g(tOr,:,:), 1));
    % And then smooth it across trials
    sgOr = cat(2, ...
        smooth(gOr(:,1), span, method), ...
        smooth(gOr(:,2), span, method));
                
    % Then shift the data so the smoothed 0 is 0
    gZeroed = g - shiftdim(sgOr, -1);    % 纠正眼位置偏移 %每次trial都依照自己看着原点的位置纠偏
       
    gShift = squeeze(median(gZeroed(tOr,:,:), 1));
   
             
    %%% plot对比了纠偏前后，在盯着原点的阶段
    figure;
    scatter3(gOr(:,1), gOr(:,2), 1:numTrials,'.')    
    hold on
%     scatter3(sgOr(:,1), sgOr(:,2), 1:numTrials)
    scatter3(gShift(:,1), gShift(:,2), 1:numTrials,'.')
    plot3(0,0, 1:numTrials,'k.',LineWidth=1.5)
%     ylim([-10 10]);
%     xlim([-6 6]);
    xlabel('X(cm)');
    ylabel('Y(cm)');
    zlabel('Trials');
    legend('gaze origin', 'gaze shift');
    title("Eye position at origin - check Zero-Correct");

    %% Rotate eye movements to vector origin to target 把所有眼动转到同一个target的方向
    target = [hand_table.target_x, hand_table.target_y];
    or = [0, 0];
        
    tgtDir = target - or;
    tgtDir = tgtDir ./ vecnorm(tgtDir,2,2);
    perpDir = [tgtDir(:,2) -tgtDir(:,1)];
        
    tarRotated = [
            dot(target, perpDir, 2) ...
            dot(target, tgtDir, 2) ...
            ];                            % 旋转后的target位置在大约（0,10)的位置
    % 旋转后的眼位置  
    gRotated = cat(3,...
                    dot(gZeroed, repmat(shiftdim(perpDir,-1), [numTVals 1 1]),3), ...
                    dot(gZeroed, repmat(shiftdim(tgtDir,-1), [numTVals 1 1]),3) ...
                    );
    
%     %%% 另一种实现方法
%     % 旋转eye
%     % 步骤 1
%     theta = atan2(hand_table.target_y, hand_table.target_x);
%     target_rotate_angle = rad2deg(theta);
%     % 以原点为圆心，计算与 y 轴正方向的夹角
%     target_rotate_angle = mod(90 - target_rotate_angle, 360); % 360度取模确保角度在[0, 360)范围内
%     % 步骤 2: 对每个 trial 的眼睛位置进行旋转
%     gazeRotated = cell(length(numTrials),1);
%     
%     for i = 1:numTrials
%         x_scaled = gaze(:, i, 1);
%         y_scaled = gaze(:, i, 2);    
%         % 逆时针旋转眼睛位置
%         rotated_gaze_trials = rotate_coordinates([x_scaled, y_scaled], target_rotate_angle(i));        
%         % 存储旋转后的坐标
%         gazeRotated{i} = rotated_gaze_trials;
%     end


end


if plot_eye
    %% 画出选定 trial 的旋转前后眼坐标轨迹
    selected_trial = 14;
    figure;
    hold on;
    % 旋转前眼坐标轨迹
    plot(gZeroed(:, selected_trial, 1), gZeroed(:, selected_trial, 2), 'bo', 'DisplayName', 'Raw Eye');
    % 旋转后眼坐标轨迹
    plot(gRotated(:, selected_trial, 1), gRotated(:, selected_trial, 2), 'r-', 'DisplayName', 'Rotated Eye');
%     plot(gazeRotated{selected_trial}(:, 1), gazeRotated{selected_trial}(:, 2), 'r-', 'DisplayName', 'Rotated Eye');
    hold off;
    % 添加标签和标题
    xlabel('X轴');
    ylabel('Y轴');
    title('选定 trial 的旋转前后眼坐标轨迹');
    legend('show');
   
end


if hand_prepro
    %% 旋转hand
    % 步骤 1
    theta = atan2(hand_table.target_y, hand_table.target_x);
    target_rotate_angle = rad2deg(theta);
    % 以原点为圆心，计算与 y 轴正方向的夹角
    target_rotate_angle = mod(90 - target_rotate_angle, 360); % 360度取模确保角度在[0, 360)范围内
    % 步骤 2
    % 提取出 movement_frames 代表的手位置的坐标点
    raw_hand = cellfun(@(str) extract_hand_position(str), hand_table.movement_frames, 'UniformOutput', false);
    % 步骤 3
    % 逆时针旋转 raw_hand 中的每个坐标
    rotated_hand = cellfun(@(hand, angle) rotate_coordinates(hand, angle), raw_hand, num2cell(target_rotate_angle), 'UniformOutput', false);
end


if plot_hand
    %% 画出所有 trials 的旋转后手坐标轨迹
    figure;
    hold on;
    % 初始化图例标记
    baseline_legend_added = false;
    rotation_legend_added = false;
    aftereffect_legend_added = false;
    
    % 循环绘制每个 trial 的旋转后手坐标轨迹
    for i = 1:size(rotated_hand, 1)
        if hand_table.baseline(i) == 1
            if ~baseline_legend_added
                baseline_hand = plot(rotated_hand{i}(:, 1), rotated_hand{i}(:, 2), '-', 'Color', 'r');
                baseline_legend_added = true;
            else
                plot(rotated_hand{i}(:, 1), rotated_hand{i}(:, 2), '-', 'Color', 'r');
            end
        elseif hand_table.rotation(i) == 1
            if ~rotation_legend_added
                rotation_hand = plot(rotated_hand{i}(:, 1), rotated_hand{i}(:, 2), '-', 'Color', 'b');
                rotation_legend_added = true;
            else
                plot(rotated_hand{i}(:, 1), rotated_hand{i}(:, 2), '-', 'Color', 'b');
            end
        elseif hand_table.aftereffect(i) == 1
            if ~aftereffect_legend_added
                aftereffect_hand = plot(rotated_hand{i}(:, 1), rotated_hand{i}(:, 2), '-', 'Color', 'k');
                aftereffect_legend_added = true;
            else
                plot(rotated_hand{i}(:, 1), rotated_hand{i}(:, 2), '-', 'Color', 'k');
            end
        end
    end
    hold off;
    % 设置坐标轴范围
    axis([-15, 15, -5, 25]);
    % 添加标签和标题
    xlabel('X');
    ylabel('Y');
    title('hand trace');
    legend([baseline_hand, rotation_hand, aftereffect_hand],'Baseline','Rotation', 'Aftereffect')
    
    %% 画出选定 trial 的旋转前后手坐标轨迹
    % selected_trial = 3;
    % figure;
    % hold on;
    % % 旋转前手坐标轨迹
    % plot(raw_hand{selected_trial}(:, 1), raw_hand{selected_trial}(:, 2), '-', 'DisplayName', 'Raw Hand');
    % % 旋转后手坐标轨迹
    % plot(rotated_hand{selected_trial}(:, 1), rotated_hand{selected_trial}(:, 2), '-', 'DisplayName', 'Rotated Hand');
    % hold off;
    % % 添加标签和标题
    % xlabel('X轴');
    % ylabel('Y轴');
    % title('选定 trial 的旋转前后手坐标轨迹');
    % legend('show');
    
    %% Trace plot
    numShowTrialsB = 16;
    numShowTrialsEA = 26;
    numShowTrialsLA = 56;
    numShowTrialsAE = 66;
    numShowTrialsAEL = 80;
    
    selectedTrialsB = 1:2:numShowTrialsB;
    selectedTrialsEA = 17:1:numShowTrialsEA;
    selectedTrialsLA = 27:3:numShowTrialsLA;
    selectedTrialsAE = 57:1:numShowTrialsAE;
    selectedTrialsAEL = 67:2:numShowTrialsAEL;
   
    figure();
    set(gcf, 'position', [0 0 1200 800]);

    %% BASELINE
    ax1 = axes;
    set(ax1, 'position', [0.13 0.11 0.775 0.815]);
    for i = 1:length(selectedTrialsB)
        thisTrialNum = selectedTrialsB(i);
        hx = rotated_hand{thisTrialNum}(:, 1);
        hy = rotated_hand{thisTrialNum}(:, 2);
        
        B = plot(hx, hy, 'color','#808080', 'LineWidth', 1); %#ffff4d
        hold on;
    end
    axis([-15 5 -2 20]);
    %% Aftereffect  
    ax2 = axes;
    set(ax2, 'position', [0.13 0.11 0.775 0.815]);
    % Initialize an array to store the custom colors for Aftereffect trials
    customColorsAE = zeros(numShowTrialsAE - 57 + 1, 3);

    for i = 1:length(selectedTrialsAE)
        thisTrialNum = selectedTrialsAE(i);
        hx = rotated_hand{thisTrialNum}(:, 1);
        hy = rotated_hand{thisTrialNum}(:, 2);
        
        % caculate the value of color of current trial
        % map color value to (0,1)
        value = (thisTrialNum - 57) / (numShowTrialsAE - 57);
        % R、G、B from 0 to 1
        thisColor = [1, value, 0];
        % Store the custom color in the array
        customColorsAE(i, :) = thisColor;

        AE = plot(hx, hy, 'color',thisColor, 'LineWidth', 1); %#77AC30
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;

    %% ADAPTATION (early)
    ax3 = axes;
    set(ax3, 'position', [0.13 0.11 0.775 0.815]);
    % Initialize an array to store the custom colors for Aftereffect trials
    customColorsEA = zeros(numShowTrialsEA - 17 + 1, 3);
    for i = 1:length(selectedTrialsEA)
        thisTrialNum = selectedTrialsEA(i);
        hx = rotated_hand{thisTrialNum}(:, 1);
        hy = rotated_hand{thisTrialNum}(:, 2);

        % caculate the value of color of current trial
        % map color value to (0,1)
        value = (thisTrialNum - 17) / (numShowTrialsEA - 17);
        % R、G、B from 0 to 1
        thisColor = [0, value, 1];
        % Store the custom color in the array
        customColorsEA(i, :) = thisColor;

        EA = plot(hx, hy, 'color',thisColor, 'LineWidth', 1); % Use a thicker blue line for late adaptation trials
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;
        %% ADAPTATION (late)
    ax4 = axes;
    set(ax4, 'position', [0.13 0.11 0.773 0.813]);
    for i = 1:length(selectedTrialsLA) % Consider the last 40 trials as early adaptation
        thisTrialNum = selectedTrialsLA(i);
        hx = rotated_hand{thisTrialNum}(:, 1);
        hy = rotated_hand{thisTrialNum}(:, 2);
        LA = plot(hx, hy, 'color',"#696969", 'LineWidth', 1,'LineStyle','--'); %#00FFFF
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;
        %% aftereffect (late 10 trials)
    ax5 = axes;
    set(ax5, 'position', [0.1 0.11 0.773 0.813]);
    for i = 1:length(selectedTrialsAEL) % Consider the last 40 trials as early adaptation
        thisTrialNum = selectedTrialsAEL(i);
        hx = rotated_hand{thisTrialNum}(:, 1);
        hy = rotated_hand{thisTrialNum}(:, 2);
        AEL = plot(hx, hy, 'color',"#696969", 'LineWidth', 1,'LineStyle',':'); %#00FFFF
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;
    scatter(0, 12, 100, 'k','filled')
    scatter(0, 0, 100, 'k','filled')
     
    % Set the custom colormap for earlyadapt
    colormap(ax3, customColorsEA);
    % Add colorbar for Aftereffect
    AColorbar = colorbar(ax3, 'Location', 'west');  
    % Modify colorbar ticks and labels
    tickValues = linspace(selectedTrialsEA(1), selectedTrialsEA(end), 10); % You can adjust the number of ticks as needed
    tickLabels = arrayfun(@(x) sprintf('%d', round(x)), tickValues, 'UniformOutput', false);
    AColorbar.Ticks = (tickValues - selectedTrialsEA(1)) / (selectedTrialsEA(end) - selectedTrialsEA(1));
    AColorbar.TickLabels = tickLabels;
    AColorbar.Label.String = 'trialNum-Early adapt';

    % Set the custom colormap for Aftereffect
    colormap(ax2, customColorsAE);
    % Add colorbar for Aftereffect
    hColorbar = colorbar(ax2, 'Location', 'east');  
    % Modify colorbar ticks and labels
    tickValues = linspace(selectedTrialsAE(1), selectedTrialsAE(end), 20); % You can adjust the number of ticks as needed
    tickLabels = arrayfun(@(x) sprintf('%d', round(x)), tickValues, 'UniformOutput', false);
    hColorbar.Ticks = (tickValues - selectedTrialsAE(1)) / (selectedTrialsAE(end) - selectedTrialsAE(1));
    hColorbar.TickLabels = tickLabels;
    hColorbar.Label.String = 'trialNum-aftereffect';
    
    hold off;
    % axis equal;
    xlim([-15 4]);
    ylim([-2 20]);
    legend([B, LA, AEL], 'Baseline', 'Late adaptation', 'Late aftereffect',Location='south');
    title('Hand trace');
    xlabel("x (cm)");
    ylabel("y (cm)");
end



% 辅助函数1：提取手位置字符串中的坐标点
function hand_position = extract_hand_position(str)
    num_str = regexp(str, '[\d.-]+', 'match');
    num_array = str2double(num_str);
    hand_position = reshape(num_array, 2, []).';
end

% 辅助函数2：逆时针旋转坐标点
function rotated_coords = rotate_coordinates(coords, angle)
    angle_rad = deg2rad(angle);
    rotation_matrix = [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)];
    rotated_coords = (rotation_matrix * coords.').';
end



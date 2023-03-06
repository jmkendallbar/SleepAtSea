%% SLEEP SPIRAL ANALYSIS
% Author: Jessica Kendall-Bar, 2022
% Script to identify and describe the occurrence of sleep spirals in
% EEG/3D Motion Data. 

% Read in hypnotrack (with 3D motion data [pitch, roll, heading] and sleep data at regular
% intervals to examine spiral behavior).
% Set data directory; change as necessary.
cd .. 
cd 'Data'

filename = '01_01_Hypnotrack_1Hz_ALL_ANIMALS.csv'; % Insert name of your Sleep & Motion 1Hz dataset 
T_ALL = readtable(filename);

SealIDs = string({'test31_FatiguedFiona', 'test33_HypoactiveHeidi', 'test35_JauntingJuliette'})
s = 1; % Going through seals one at a time
idx = contains(T_ALL.SealID, SealIDs(1), 'IgnoreCase', true); 
T = T_ALL(idx,:);

%% IDENTIFY SLEEP SPIRALS

% DEFINITIONS:
% Nap: consecutive segment of sleep
% Right Turn: turning right (diff(heading) between  0 and POS. pi or jumping to NEG. ~2pi
% Left Turn: turning left (diff(heading) between  0 and NEG. pi or jumping to POS. ~2pi
% Right Spin: turning right past 180 causing jump of NEG. 2pi
% Left Spin: turning left past 180 causing jump of POS. 2pi
% Spiral: two (or more) consecutive spins past 180 in the same direction
% Loop: a single loop of a spiral (between two same-direction spins past 180)

T.is_sleep = contains(T.Simple_Sleep_Code, ({'SWS', 'REM'}), 'IgnoreCase', true); % if SWS or REM
T.is_REM   = contains(T.Simple_Sleep_Code, ({'REM'}), 'IgnoreCase', true); % if REM, assign 1

% Calculate heading change (1st derivative of heading)
T.headdiff  = [0; diff(T.heading)]; 
T.RightSpin = T.headdiff < -pi ; % a right spin past 180 South is a diff of ~ -2*pi
T.LeftSpin  = T.headdiff > pi ;  % a left spin past 180 South is a diff of ~ +2*pi

% Creating heading column that does not jump between 180 and -180 
% 1. Create cumulative sum of right turns and left turns to keep track of
% overall turning (past 180) to the left or right.
T.CumulTurns_Rpos = cumsum(T.RightSpin - T.LeftSpin); 
% 2. Correct heading to add 2*pi to the heading for every right spin past 180 (and
% subtract 2*pi from the heading for every left spin past 180).
T.headcorr_CumulTurns = T.heading + 2*pi*(T.CumulTurns_Rpos);
% 3. Check that that gets rid of large jumps in heading for smooth animations.
T.headcorrdiff = [diff(T.headcorr_CumulTurns); 0];
max(T.headcorrdiff)

% Defining turning left and turning right
% LEFT TURNS: either slow (< pi/s or < 180 degrees/second) turns to the left (negative)
% OR sudden jumps (> 180 deg/s) "to the right" (positive)
T.TurningLeftCriteria = ((T.headdiff <= 0 & T.headdiff >= -pi) | ...
                        (T.headdiff > pi ));
% RIGHT TURNS: either slow (< pi/s or < 180 degrees/second) turns to the right (positive)
% OR sudden jumps (> 180 deg/s) "to the left" (negative)
T.TurningRightCriteria = ((T.headdiff <= pi & T.headdiff >= 0) | ...
                        (T.headdiff < -pi ));                        

% Create table with each consecutive (OCEAN) nap     
notonland           = ~contains(T.Water_Code, ({'LAND'}), 'IgnoreCase', true);
Nap_Criteria        = T.is_sleep & notonland; % Looking for naps NOT on land
Naps                = table(yt_setones(Nap_Criteria),'VariableNames',{'Indices'});
Naps.Duration_s     = (Naps.Indices(:,2)-Naps.Indices(:,1));
Naps                = Naps(find(Naps.Duration_s~=0),:);

T.SleepSpinNum(:) = nan;
for d = 1:height(Naps)
    startix = Naps.Indices(d,1);
    endix = Naps.Indices(d,2);
    duration = (endix-startix)+1;
    
    [uniqueSleeps, ~, idx] = unique(T.Simple_Sleep_Code(startix:endix));
    counts = accumarray(idx, 1);
    [~, maxIndex] = max(counts);
    modeSleep = uniqueSleeps(maxIndex);

    Naps.MostlySimpleSleepNum(d) = modeSleep; % Most common sleep stage for the nap period
    Naps.REM_seconds(d) = sum(T.is_REM(startix:endix)); % seconds spent in REM in the nap
    Naps.REM_percentage(d) = sum(T.is_REM(startix:endix))/sum(T.is_sleep(startix:endix));  % percent of nap spent in REM
    T.standardsleepxposition(startix:endix) = T.x(startix:endix)-T.x(startix) + d*20; % zero'd x position plus offset
    T.standardsleepyposition(startix:endix) = T.y(startix:endix)-T.y(startix); % zero'd y position 
    T.standardsleepzposition(startix:endix) = T.z(startix:endix)-T.z(startix); % zero'd z position
    Naps.RightSpins(d) = sum(T.RightSpin(startix:endix)); % count number of right spins past 180 during nap
    Naps.LeftSpins(d)  = sum(T.LeftSpin(startix:endix));  % count number of left spins past 180 during nap
    Naps.OverallSpins(d) = Naps.RightSpins(d) - Naps.LeftSpins(d); % overall spins right or left
    % Creating a new column in T where 0 is no turn, positive is a right spin past 180
    % and negative is a left spin past 180
    T.SleepSpinNum(startix:endix) = cumsum(T.RightSpin(startix:endix))-cumsum(T.LeftSpin(startix:endix));
end

T.diffSleepSpinNum = [0; diff(T.SleepSpinNum)];

% FIND SPIRALS
% Criteria: Find consecutive chunks where there are no spins or spins (past 180) to
% the left. Filter results by showing only spirals with at least 2 left spins.
LeftSpirals = table(yt_setones(T.diffSleepSpinNum <= 0 & T.diffSleepSpinNum >= -1),'VariableNames',{'Indices'});
LeftSpirals.Duration_s   = (LeftSpirals.Indices(:,2)-LeftSpirals.Indices(:,1));
LeftSpirals.direction(:) = {'left'};

RightSpirals = table(yt_setones(T.diffSleepSpinNum >= 0 & T.diffSleepSpinNum <= 1),'VariableNames',{'Indices'});
RightSpirals.Duration_s   = (RightSpirals.Indices(:,2)-RightSpirals.Indices(:,1));
RightSpirals.direction(:) = {'right'};

Spirals = vertcat(LeftSpirals,RightSpirals);
Spirals.Start_Turns    = T.SleepSpinNum(Spirals.Indices(:,1));
Spirals.End_Turns      = T.SleepSpinNum(Spirals.Indices(:,2));
Spirals.totalturns     = Spirals.End_Turns - Spirals.Start_Turns;
Spirals                = Spirals(find(abs(Spirals.totalturns) > 1),:);
Spirals.Start_Depth    = T.Depth(Spirals.Indices(:,1));
Spirals.End_Depth      = T.Depth(Spirals.Indices(:,2));
Spirals.Start_SleepCode    = T.Simple_Sleep_Code(Spirals.Indices(:,1));
Spirals.End_SleepCode      = T.Simple_Sleep_Code(Spirals.Indices(:,2));

scatter3(T.standardsleepxposition, T.standardsleepyposition, T.standardsleepzposition,...
    [],T.Simple_Sleep_Num,'filled')

T.standardspiralxposition(:) = nan;
T.standardspiralyposition(:) = nan;
T.standardspiralzposition(:) = nan;

for d = 1:height(Spirals)
    startix = Spirals.Indices(d,1);
    endix = Spirals.Indices(d,2);
    duration = (endix-startix)+1;
    Spirals.MostlySimpleSleepNum(d) = mode(T.Simple_Sleep_Code(startix:endix), 'all');
    Spirals.REM_seconds(d) = sum(T.is_REM(startix:endix));
    Spirals.REM_percentage(d) = sum(T.is_REM(startix:endix))/sum(T.is_sleep(startix:endix));
    Spirals.sleep(d) = sum(T.is_sleep(startix:endix));
    T.is_SleepSpiral(startix:endix) = 1;
    T.SleepSpiralDirection(startix:endix) = {Spirals.direction(d)};
    T.SleepSpiralTurns(startix:endix) = Spirals.totalturns(d);
    T.SleepSpiralDuration_s(startix:endix) = Spirals.Duration_s(d);
    T.standardspiralxposition(startix:endix) = T.x(startix:endix)-T.x(startix) + d*20;
    T.standardspiralyposition(startix:endix) = T.y(startix:endix)-T.y(startix);
    T.standardspiralzposition(startix:endix) = T.z(startix:endix)-T.z(startix);
end

scatter3(T.standardspiralxposition, T.standardspiralyposition, T.standardspiralzposition,...
    [],T.Simple_Sleep_Num,'filled')

Loop_Criteria        = T.is_SleepSpiral & T.diffSleepSpinNum == 0; % Looking for curls in sleep spirals (between areas where a spin is detected (diffSleepSpinNum==1))
Loops                = table(yt_setones(Loop_Criteria),'VariableNames',{'Indices'});
Loops.Duration_s     = (Loops.Indices(:,2)-Loops.Indices(:,1));
Loops                = Loops(find(Loops.Duration_s~=0),:);

T.standardloopxposition(:) = nan;
T.standardloopyposition(:) = nan;
T.standardloopzposition(:) = nan;

for d = 1:height(Loops)
    startix = Loops.Indices(d,1);
    endix = Loops.Indices(d,2);
    duration = (endix-startix)+1;
    Loops.MostlySimpleSleepNum(d) = mode(T.Simple_Sleep_Num(startix:endix));
    Loops.REM_seconds(d) = sum(T.is_REM(startix:endix));
    Loops.REM_percentage(d) = sum(T.is_REM(startix:endix))/sum(T.is_sleep(startix:endix));
    Loops.sleep(d) = sum(T.is_sleep(startix:endix));
    Loops.mean_speed(d) = mean(T.speed(startix:endix));
    Loops.diameter(d) = (Loops.mean_speed(d) * Loops.Duration_s(d))/pi;
    T.LoopNum(startix:endix) = d;
    T.LoopDur(startix:endix) = Loops.Duration_s(d);
    T.LoopModeSleepCode(startix:endix) = Loops.MostlySimpleSleepNum(d);
    T.standardloopxposition(startix:endix) = T.x(startix:endix)-T.x(startix);
    T.standardloopyposition(startix:endix) = T.y(startix:endix)-T.y(startix);
    T.standardloopzposition(startix:endix) = T.z(startix:endix)-T.z(startix) + d*20;
end

Loops_Only = T(find(~isnan(T.standardloopxposition)),:);
Naps_Only = T(find(~isnan(T.standardsleepxposition)),:);
Spirals_Only = T(find(~isnan(T.standardspiralxposition)),:);

scatter3(T.standardloopxposition, T.standardloopyposition, T.standardloopzposition,...
    [],T.Simple_Sleep_Num,'filled')

% Save data if necessary, compile T across seals to get 01_01_Hypnotrack_1Hz_ALL_ANIMALS.csv
% writetable(Spirals,strcat(SealIDs(s),'_09_Spirals_Stats.csv'));
% writetable(Loops,strcat(SealIDs(s),'_09_Loops_Stats.csv'));
% writetable(Loops_Only,strcat(SealIDs(s),'_09_Hypnotrack_1Hz_All_Loops.csv'));
% writetable(Spirals_Only,strcat(SealIDs(s),'_09_Hypnotrack_1Hz_All_Sleep_Spirals.csv'));
% writetable(Naps_Only,strcat(SealIDs(s),'_09_Hypnotrack_1Hz_All_Naps.csv'));
% writetable(T,strcat(SealIDs(s),'_09_Hypnotrack_1Hz_withSpirals.csv'));

%% READ IN SPIRAL DATA FOR ALL SEALS TO SUMMARIZE SLEEP SPIRAL DATA ACROSS ANIMALS

T_ALL = readtable('01_01_Hypnotrack_1Hz_ALL_ANIMALS.csv');

% Find dives
% Criteria: dives must be continuous segments deeper than dive_threshold.
T_ALL.is_dive          = abs(T_ALL.Depth) > 2;
Dives                   = table(yt_setones(T_ALL.is_dive),'VariableNames',{'Indices'}); 

Dives.Duration_s        = (Dives.Indices(:,2)-Dives.Indices(:,1))*1;
% Eliminate dives with duration 0 or longer than 500 min (false IDs)
Dives                   = Dives(find(Dives.Duration_s~=0),:);

for d = 1:height(Dives)
        startix = Dives.Indices(d,1);
        endix = Dives.Indices(d,2);
        duration = (endix-startix)+1;

        T_ALL.is_DiveNum(startix:endix) = d;
end

% Get dive indices of all dives with sleep spirals
spiraldives = T_ALL(find(T_ALL.is_SleepSpiral),:);
spiraldivenums = table(unique(spiraldives.is_DiveNum));

% Get dive indices of all REM sleep dives
REM_OO = T_ALL(find(T_ALL.Water_Num==3 & T_ALL.is_REM),:);
REM_OO_Divenums = table(unique(REM_OO.is_DiveNum));

% Get dive indices of all REM episodes that don't overlap with sleep spirals
nonSpiralREM = T_ALL(find(T_ALL.Water_Num==3 & T_ALL.is_REM & T_ALL.is_SleepSpiral==0),:);
nonSpiralREMdivenums = table(unique(nonSpiralREM.is_DiveNum));
for d = 1:height(nonSpiralREMdivenums)
        nonSpiralREMdivenums.is_DiveNum(startix:endix) = d;
end
T_ALL.is_SleepSpiralDive = ismember(T_ALL.is_DiveNum, spiraldivenums.Var1);

isPresent = any(ismember(T_ALL, [x y], 'rows'));

% Code generates total seconds in each state for a certain animal 
% Because frequency is 1Hz, adding together # of rows / 3600 returns hours

% Total recording time in Open Ocean
OO_h = length(find(T_ALL.Water_Num==3))/3600;
% OPEN OCEAN sleep in minutes
Sleep_OO_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_sleep==1))/60;
% Non-sleep spirals
NonSleep_Spirals_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_sleep==0 & T_ALL.is_SleepSpiral))/60;

% OPEN OCEAN SWS in minutes
SWS_OO_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_sleep==1 & T_ALL.is_REM==0))/60;
% OPEN OCEAN SWS within SLEEP SPIRALS in minutes
SWS_OO_SS_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_sleep==1 & T_ALL.is_REM==0 & T_ALL.is_SleepSpiral))/60;
% OPEN OCEAN SWS NOT within SLEEP SPIRALS, but in dives with sleep spirals in minutes
SWS_OO_SSD_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_sleep & T_ALL.is_REM==0 & T_ALL.is_SleepSpiral==0 & T_ALL.is_SleepSpiralDive))/60;
% OPEN OCEAN SWS while upside down in minutes
SWS_OO_SU_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_sleep==1 & T_ALL.is_REM==0 & abs(T_ALL.roll) > 2))/60;

% OPEN OCEAN REM in minutes
REM_OO_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_REM))/60;
% OPEN OCEAN REM within SLEEP SPIRALS in minutes
REM_OO_SS_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_REM & T_ALL.is_SleepSpiral))/60;
% OPEN OCEAN REM NOT within SLEEP SPIRALS, but in dives with sleep spirals in minutes
REM_OO_SSD_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_REM & T_ALL.is_SleepSpiral==0 & T_ALL.is_SleepSpiralDive))/60;

% OPEN OCEAN REM in minutes
REM_OO_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_REM))/60;
% OPEN OCEAN REM within SLEEP SPIRALS in minutes
REM_OO_SS_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_REM & T_ALL.is_SleepSpiral))/60;
% OPEN OCEAN REM while upside down in minutes
REM_OO_SU_min = length(find(T_ALL.Water_Num==3 & T_ALL.is_REM & abs(T_ALL.roll)>2))/60;

CS_h = length(find(T_ALL.Water_Num==2))/3600;
% CONTINENTAL SHELF sleep in minutes
Sleep_CS_min = length(find(T_ALL.Water_Num==2 & T_ALL.is_sleep==1))/60;
% CONTINENTAL SHELF SWS in minutes
SWS_CS_min = length(find(T_ALL.Water_Num==2 & T_ALL.is_sleep==1 & T_ALL.is_REM==0))/60;
% CONTINENTAL SHELF SWS while upside down in minutes
SWS_CS_SU_min = length(find(T_ALL.Water_Num==2 & T_ALL.is_sleep==1 & T_ALL.is_REM==0 & abs(T_ALL.roll) > 2))/60;

% CONTINENTAL SHELF REM in minutes
REM_CS_min = length(find(T_ALL.Water_Num==2 & T_ALL.is_REM))/60;
% CONTINENTAL SHELF REM while upside down in minutes
REM_CS_SU_min = length(find(T_ALL.Water_Num==2 & T_ALL.is_REM & abs(T_ALL.roll)>2))/60;

sleepspirals_min = length(find(T_ALL.is_SleepSpiral))/60;
sleepspiralSWS_min = length(find(T_ALL.is_SleepSpiral & T_ALL.is_sleep & T_ALL.is_REM==0))/60;

T_ALL.DriftRate = [0 ; diff(T_ALL.Depth)]; 
REM_Drifts = T_ALL.DriftRate(find(T_ALL.Water_Num==3 & T_ALL.is_REM));
SWS_Drifts = T_ALL.DriftRate(find(T_ALL.Water_Num==3 & T_ALL.is_sleep & T_ALL.is_REM==0));

REM_CS = length(find(T_ALL.Water_Num==2 & T_ALL.is_REM & abs(T_ALL.DriftRate) < 0.1));
SWS_CS = length(find(T_ALL.Water_Num==2 & T_ALL.is_sleep & T_ALL.is_REM==0 & abs(T_ALL.DriftRate) < 0.1));

sleep_CS_flats = (SWS_CS + REM_CS) / 3600 % HOURS OF SLEEP on the continental shelf
sleep_CS = (REM_CS_min + SWS_CS_min) / 60 % HOURS OF SLEEP on or above the continental shelf
sleep_CS_flats/sleep_CS

percentSWS_onflats = SWS_CS / (SWS_CS + REM_CS)
percentREM_onflats = REM_CS / (SWS_CS + REM_CS)

mean(REM_Drifts)
std(REM_Drifts)

mean(SWS_Drifts)
std(SWS_Drifts)

figure
histogram(SWS_Drifts,-1.3:.01:1.3,'facecolor','blue','facealpha',.5,'edgecolor','none')
hold on
histogram(REM_Drifts,-1.3:.01:1.3,'facecolor','yellow','facealpha',.5,'edgecolor','none')
box off
axis tight
legend('SWS Drift Rates','REM Drit Rates','location','northwest')
legend boxoff

strcat('% SWS Continental Shelf:  ',...
    num2str((SWS_CS_min/Sleep_OO_min)*100),' %.')
strcat('% REM Continental Shelf:  ',...
    num2str((REM_CS_min/Sleep_OO_min)*100),' %.')
strcat('% SWS Open Ocean:  ',...
    num2str((SWS_OO_min/Sleep_OO_min)*100),' %.')
strcat('% REM Open Ocean:  ',...
    num2str((REM_OO_min/Sleep_OO_min)*100),' %.')

strcat('% SWS in Sleep Spirals in Open Ocean: ',...
    num2str((SWS_OO_SS_min/SWS_OO_min)*100),...
    ' %.')
strcat('% REM in Sleep Spirals in Open Ocean: ',...
    num2str((REM_OO_SS_min/REM_OO_min)*100),...
    ' %.')
strcat('% SWS in Sleep Spirals or Dives with Sleep Spirals in Open Ocean: ',...
    num2str(((SWS_OO_SS_min+SWS_OO_SSD_min)/SWS_OO_min)*100),...
    ' %.')
strcat('% REM in Sleep Spirals or Dives with Sleep Spirals in Open Ocean: ',...
    num2str(((REM_OO_SS_min+REM_OO_SSD_min)/REM_OO_min)*100),...
    ' %.')

strcat('% upside down during SWS in Open ocean: ',...
    num2str((SWS_OO_SU_min/SWS_OO_min)*100),...
    ' %.')
strcat('% upside down during REM in Open ocean: ',...
    num2str((REM_OO_SU_min/REM_OO_min)*100),...
    ' %.')

strcat('% upside down during SWS in Continental Shelf: ',...
    num2str((SWS_CS_SU_min/SWS_CS_min)*100),...
    ' %.')
strcat('% upside down during REM in Continental Shelf: ',...
    num2str((REM_CS_SU_min/REM_CS_min)*100),...
    ' %.')

strcat('Total open ocean spiral REM time: ',...
    num2str(length(find(T_ALL.Water_Num==3 & T_ALL.is_SleepSpiral & T_ALL.is_REM))/60),...
    ' minutes.')
strcat('Total open ocean non-spiral REM time: ',...
    num2str(length(find(T_ALL.Water_Num==3 & T_ALL.is_SleepSpiral==0 & T_ALL.is_REM))/60),...
    ' minutes.')
strcat('Total open ocean non-spiral REM time (excluding dives with sleep spirals): ',...
    num2str(length(find(T_ALL.Water_Num==3 & T_ALL.is_SleepSpiral==0 & T_ALL.is_REM & T_ALL.is_SleepSpiralDive==0))/60),...
    ' minutes.')
length(find(T_ALL.Water_Num==3 & T_ALL.is_SleepSpiral==0 & T_ALL.is_REM & T_ALL.is_SleepSpiralDive==0))/60
nonSpiralREM = T_ALL(find(T_ALL.Water_Num==3 & T_ALL.is_REM & T_ALL.is_SleepSpiral==0),:);

strcat(SealIDs(s), ' Total sleep: ',...
    num2str(length(find(T.is_sleep==1))/3600),...
    ' hours.')
length(find(T.is_sleep==0))


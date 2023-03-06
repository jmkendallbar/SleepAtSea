%% CUSTOM SPEED ESTIMATION for 3D Visualization
% Author: Jessica Kendall-Bar

%% Read in data and define assumptions and constants

% Set data directory; change as necessary.
Data_path='..\Data';
cd(Data_path);

% Get Stroke Rate and Guestimate speed
% These are assumed values based on previous studies where velocity was measured.
minstrokerate = 10; % Assume minimum stroke rate of 10 strokes per minute for swimming to contribute to forward speed.
maxstrokerate = 80; % Assume maximum stroke rate of 80 strokes per minute (associate this and above with max speed).
maxswimspeed = 2; % Assume that max swimming speed is 2 m/s when stroking at >= maxstrokerate.
minswimspeed = 1; % Assume that min swimming speed is 1 m/s when stroking at <= minstrokerate.
driftspeed = 0.2; % Assume forward speed of 0.2 m/s during drift segments.
bottomspeed = 0; % Speed when animal is on the bottom, not moving.
standardlength = 150; % REPLACE with your animal's standard length (in centimeters)

% Use a **1Hz** dataset that has: 
%       - Depth (in meters)
%       - Pitch (in degrees)
%       - Roll (in degrees)
%       - Stroke_Rate (column with stroke rate in strokes per minute)
rates_file = 'YOUR FILENAME HERE'; 

% Read in data
Rates = readtable(rates_file,opts);
Rates.Sec(:) = round(linspace(0,height(Rates),height(Rates)));
Rates.Stroke_Rate(find(isnan(Rates.Stroke_Rate))) = 0; % replacing NaNs with 0s

Rates.smoothDepth = smoothdata(Rates.Depth,'gaussian',20); % Moving average across 20 depth samples
Rates.smoothDiffDepth = [-diff(Rates.smoothDepth); 0]; % slope; meters per sec
vertspeed = Rates.smoothDiffDepth; % vertical speed
diagspeed = vertspeed(:) ./ sind(Rates.pitch(:)); % diagonal speed
Rates.speed0 = NaN(height(Rates),1);

% IF Swimming:
swimming = Rates.Stroke_Rate <= maxstrokerate & Rates.Stroke_Rate >= minstrokerate;
% Set maxstrokerate (and above) to maxswimspeed
Rates.speed0(find(Rates.Stroke_Rate >= maxstrokerate)) = maxswimspeed;
% Map maxstrokerate to maxswimspeed and minstrokerate to minswimspeed.
Rates.speed0(find(swimming)) = ((Rates.Stroke_Rate(find(swimming)) - minstrokerate) ...
             * (maxswimspeed - minswimspeed) / (maxstrokerate - minstrokerate)) + minswimspeed;
% Mapping maxstrokerate of 80 to 2 m/s and minstrokerate of 1 to 1 m/s
% ((oldValue - oldMin) * newRange / oldRange) + newMin
         
% IF Gliding:
gliding = Rates.Stroke_Rate <= minstrokerate;
depth_threshold = abs(Rates.smoothDepth) >= 5; % Can't be on the surface
% If gliding up or down upon ascent or descent (with high pitch), make speed = 0.8 m/s
Rates.speed0(find(gliding & abs(Rates.pitch)>40 & abs(Rates.roll)<150 & depth_threshold)) = 0.8;
Rates.speed0(find(gliding & abs(Rates.pitch)<40 & abs(Rates.roll)<150 & depth_threshold)) = 0.8;
% If drifting upside down (causing pitch to be small and positive), make speed 0.2 m/s.
drift_threshold = abs(Rates.smoothDiffDepth) <= 3 & abs(Rates.smoothDiffDepth) >= 0.1; % change from 0.5 for Heidi
bottom_threshold = abs(Rates.smoothDiffDepth) <= 0.1;
Rates.speed0(find(gliding & abs(Rates.pitch)<40 & abs(Rates.roll)>150 & drift_threshold & depth_threshold)) = 0.2;
Rates.speed0(find(gliding & abs(Rates.pitch)<20 & bottom_threshold & depth_threshold)) = 0;

% IF on land:
% Use shallow depth and low pitch to find periods on land (usually high
% pitch in the water at surface)
land_threshold = abs(Rates.smoothDepth) <= 5 & abs(Rates.pitch)<20;
% If no galumphs are detected, set speed to zero.
Rates.speed0(find(gliding & land_threshold)) = 0;
% IF GALUMPHING set max speed:
galumphing = swimming & land_threshold;

% 0.12–0.71 body lengths s−1 - We are using average standard (straight) body length of 172 
bodylength = standardlength/100; % Body length in meters
maxgalumphspeed = 0.12*bodylength; % VALUES FROM https://journals.biologists.com/jeb/article/221/18/jeb180117/19448/Terrestrial-locomotion-of-the-northern-elephant
mingalumphspeed = 0.71*bodylength; % VALUES FROM https://journals.biologists.com/jeb/article/221/18/jeb180117/19448/Terrestrial-locomotion-of-the-northern-elephant

% Set maxstrokerate (and above) to maxswimspeed
Rates.speed0(find(land_threshold & Rates.Stroke_Rate >= maxstrokerate)) = maxgalumphspeed;
% Map maxstrokerate to maxswimspeed and minstrokerate to minswimspeed.
Rates.speed0(find(galumphing)) = ((Rates.Stroke_Rate(find(galumphing)) - minstrokerate) ...
             * (maxgalumphspeed - mingalumphspeed) / (maxstrokerate - minstrokerate)) + mingalumphspeed;

% Adds interpolated speed estimate where no good estimate parameter exists.
Rates.speed0 = double(fixgaps(Rates.speed0));
Rates.smoothspeed0 = smoothdata(Rates.speed0,'gaussian',50);
% Rates.smoothspeed0 = sgolayfilt(Rates.speed0,5,21); Introduces artifacts
% for estimates with abrupt shifts

ax1=subplot(5,1,[1:2]);
plot(ax1,Rates.Sec, Rates.Depth);
ax1.YDir='reverse';
title([SealIDs(s) 'Speed Plots']);
ylabel('Depth (m)');
xlabel('Seconds');
hold on

ax2=subplot(5,1,3);
plot(ax2,Rates.Sec, Rates.speed0);
ylabel('Speed0');
xlabel('Seconds');
hold on

ax2=subplot(5,1,3);
plot(ax2,Rates.Sec, Rates.smoothspeed0);
ylabel('Speed0');
xlabel('Seconds');
hold on

ax3=subplot(5,1,4);
plot(ax3,Rates.Sec, Rates.pitch, 'Color', [1 0 0 0.3])
ylabel('Degrees');
hold on
plot(ax3,Rates.Sec, Rates.roll, 'Color', [0 1 0 0.3])

ax4=subplot(5,1,5);
plot(ax4,Rates.Sec, Rates.Stroke_Rate, 'Color', [0 0.7 0.25 0.3])
ylabel('Stroke Rate (spm)');

legend(ax2,'Speed Estimate','Smoothed Speed Estimate')
legend(ax3,'Pitch','Roll')

ylim(ax2, [0 2.5])
ylim(ax4, [-10 max(Rates.Stroke_Rate)])
linkaxes([ax1,ax2,ax3,ax4],'x');

%% LINE UP DATA WITH CATS TOOLBOX DATA 
% BEFORE LINKING THESE UP - make sure the beginning and end of your Rates file is consistent with
% ON.ANIMAL and OFF.ANIMAL tagon times (in CATS Toolbox) - will use these to match data back to prh file
fs = 5; % sampling rate of CATS Toolbox data
speed_5hz = interp(Rates.smoothspeed0,fs);
manualspeed = NaN(height(tagon),1);
manualspeed(tagon) = speed_5hz(1:height(tagon(tagon)));

if ~exist('speed','var') 
    speed = table(manualspeed);
else
    speed.manualspeed(:) = manualspeed;
end

% Return to CATS Toolbox & run step 13 (once you also have GPS data).
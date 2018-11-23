%% Projectile/ Motion Visualiser
%% Written by Benjamin Langford
%% To Do

% Take Vectors in for a projectile in 3 dimensions
% Add them together to get a resultant vector
%Plot the motion along a 3d graph
%

%% Housekeeping
clc;
clear;
close all;

%% Inputs

disp('For each vector acting on the projectile, please input the Magnitude of Velocity, Bearing and Angle');
exit = 0;
n = 1;

inputs = struct;
while exit == 0
    %Vector type
    tt = 0;
    while tt == 0
        type = input(sprintf("Type of Vector %.0f, (a)cceleration or (i)nitial: ",n),'s');
        if(isscalar(type) & type == 'a' | type == 'i') %#ok<OR2,AND2>
            tt = 1;
            inputs.types(n) = type;
        else
            disp('Invalid input, please enter i for an intial vector and a for an acceleration vector');
        end
    end
    
    %Magnitude
    mt = 0;
    while mt == 0
        magnitude = input(sprintf("Enter magnitude for Vector %.0f: ", n), 's');
        if isnan(str2double(magnitude))
            disp("Invalid input, please enter a number");
        else
            inputs.magnitudes(n) = str2double(magnitude);
            mt = 1;
        end
    end
    
    %Bearing
    bt = 0;
    while bt == 0
        bearing = input(sprintf("Enter bearing for Vector %.0f (Horizontal angle): ", n), 's');
        if isnan(str2double(bearing))
            disp("Invalid input, please enter a number");
        else
            inputs.bearings(n) = str2double(bearing);
            bt = 1;
        end
    end
    %Angle
    at = 0;
    while at == 0
        angle = input(sprintf("Enter an angle for Vector %.0f (Vertical angle): ", n), 's');
        if isnan(str2double(angle))
            disp("Invalid input, please enter a number");
        else
            inputs.angles(n) = str2double(angle);
            at = 1;
        end
    end
    
    %Done?
    yt = 0;
    while yt == 0
        finished = input("Any more vectors? (y or n): ", 's');
        switch upper(finished)
            case 'N'
                exit = 1;
                yt = 1;
            case 'Y'
                n = n + 1;
                yt = 1;
            otherwise
                disp("Invalid input, please enter y for Yes or n for No");
        end
    end
end

%% Processing

disp("Processing...");
initials = struct;
accelerations = struct;

%Gravity
accelerations(1).magnitude = 9.81;
accelerations(1).bearing = 0;
accelerations(1).angle = 270;

% Sorting
in = 1;
an = 2;
for i = 1:n
    if inputs.types(i) == 'i'
        initials(in).magnitude = inputs.magnitudes(i);
        initials(in).bearing = inputs.bearings(i);
        initials(in).angle = inputs.angles(i);
        in = in + 1;
    elseif inputs.types(i) == 'a'
        accelerations(an).magnitude = inputs.magnitudes(i);
        accelerations(an).bearing = inputs.bearings(i);
        accelerations(an).angle = inputs.angles(i);
        an = an + 1;
    end
end
%Initial Components
i_Components = struct;

if in ~= 1
    for i = 1:(in-1)
        i_Components(i).x = initials(i).magnitude * cosd(initials(i).angle) * sind(initials(i).bearing);
        i_Components(i).y = initials(i).magnitude * cosd(initials(i).angle) * cosd(initials(i).bearing);
        i_Components(i).z = initials(i).magnitude * sind(initials(i).angle);
    end
end

%Acceleration Components
a_Components = struct;

for i = 1:(an-1)
    a_Components(i).x = accelerations(i).magnitude * cosd(accelerations(i).angle) * sind(accelerations(i).bearing);
    a_Components(i).y = accelerations(i).magnitude * cosd(accelerations(i).angle) * cosd(accelerations(i).bearing);
    a_Components(i).z = accelerations(i).magnitude * sind(accelerations(i).angle);
end

% Resultant calculation
i_resultant.x = 0;
i_resultant.y = 0;
i_resultant.z = 0;
a_resultant.x = 0;
a_resultant.y = 0;
a_resultant.z = 0;
% Initial Resultant
if in ~= 1
    for i = 1:(in-1)
        i_resultant.x = i_resultant.x + i_Components(i).x;
        i_resultant.y = i_resultant.y + i_Components(i).y;
        i_resultant.z = i_resultant.z + i_Components(i).z;
    end
end
%Acceleration Resultant
for i = 1:(an-1)
    a_resultant.x = a_resultant.x + a_Components(i).x;
    a_resultant.y = a_resultant.y + a_Components(i).y;
    a_resultant.z = a_resultant.z + a_Components(i).z;
end

% Maximum time calculation


% Displacement Calculation

disp("Done!");
%% Displaying
% Graph initialisation
%hold all;
%grid on;
%view(3);
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

disp('For each vector acting on the projectile, please input the Magnitude, Bearing and Angle');
exit = 0;
n = 1;

inputs = struct;
while exit == 0
    %Vector type
    tt = 0;
    while tt == 0
        type = input(sprintf("Type of Vector %.0f, (c)onstant or (i)nitial: ",n),'s');
        if(type == 'c' || type == 'i')
            tt = 1;
            inputs.types(n) = type;
        else
            disp('Invalid input, please enter i for an intial vector and c for a constant vector');
        end
    end
    %Magnitude
    mt = 0;
    while mt == 0;
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

%% Displaying
% Graph initialisation
hold all;
grid on;
view(3);
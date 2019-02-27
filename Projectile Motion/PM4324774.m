%% Projectile/ Motion Visualiser
%% Written by Benjamin Langford
% This program is designed to take input from the user and solve a set of
% linear motion problems using SUVAT equations, specifically for 4 cases
% that are typical SUVAT problems.

%% Housekeeping
clc; % This command clears the command window of all text
clear; % This command clears all variables stored in memory. Not entirely necessary for reasons I will discuss later, but still good practice.
close all; %This command closes any previously initiated figures. I do this since any new plots I create will draw themselves on previous figures.
%% Case Selection
% I start by asking the user which case they want to calculate and then
% detailing the cases that this program can support. It's useful to note
% that I use the disp function rather than fprintf since disp finishes its
% message with a new line, so that I don't have to use \n in each message
disp("Please enter a number for which case you want to calculate.");
disp("1. Plot the path of a ball given a set of initial velocities");
disp("2. Calculate the horizontal displacement of a ball released from a height, given the height, angle, and initial speed");
disp("3. Calculate the initial height of a ball and the initial speed given the time and displacement");
disp("4. Calculcate the initial speed and angle of a ball given the time and displacement");

% I then take the input for which case the user wants to calculate. The
% while loop below is the typical error checking layout that I will use for
% most of my inputs.
ct = 0; % I initialise a variable to act as the exit condition for the while loop.
while ct == 0 % I then to initialise the while loop, telling it to run while the exit condition is a certain value.
    c = input("Which case?: ", 's');
    %         I then take the input from the user. In the input() function, I
    %         tell it to take the input as a string using the 's' identifier.
    %         This is so that I can error check it in the line below. Taking it
    %         in as a string means the user can input whatever they want
    %         without the input function comlaining. The standard error check
    %         that I run is the following: Is the input a number and is the
    %         input a real number. I do this by checking if the ouput if the
    %         str2double function is NaN (Not a Number), if true then I know
    %         it's invalid. I check for complex numbers by checking if the
    %         complex part of the input is greater than zero, if so, then I
    %         know it's invalid. I then check if the input is within the
    %         allowed boundaries. I surround all these checks with OR
    %         operators, since if any of these checks return true, then the
    %         entire bracket will return true and I know the input is invalid.
    
    % NOTE: I know the method I use below isn't the most efficient with
    % converting it to a double each time I check, but MATLAB is so efficient
    % and this script is so small, the extra time taken is inconsequential.
    if (isnan(str2double(c))  || str2double(c) > 4 || abs(imag(str2double(c))) > 0 || str2double(c) < 1)
        disp("Invalid input, please input 1, 2, 3, or 4 depending on which case you want to calculate from above");
        % If the inputis invalid, I inform the user and let them know what an
        % appropriate input would be.
    else
        c = str2double(c); % If the input is valid then I convert the string input to it's numerical counterpart and store it.
        ct = 1; % I then change the condition variable  to exit the while loop.
    end
end

% I then navigate to the appropriate case using a switch statement, which
% is basically a condensed if statement.
switch c
    case 1
        caseOne() % I actually call each case as a function. While nesting the case code within an if statement
%         would achieve the same effect, I instead use functions out of
%         preference, since in the MATLAB IDE I can collapse the contents
%         of a function and ignore it while I work on another.
    case 2
        caseTwo()
    case 3
        caseThree()
    case 4
        caseFour()
end

% NOTE: One advantage of me doing the cases like this is that I can easily
% add more cases by creating new ones within a new function and calling
% that function here.
%% Cases
% Below are the functions that hold the code for each case. I mentioned
% earlier that using the clear command was unnecessary and this is why. The
% variables that are created within a function get deleted once a function
% is closed, so the only variables that are left when the program is done
% are the case variables. Technically the program ends on line 70 and the
% code below is accessed separately.
function caseOne
%% Inputs
% In this case, the program takes in a set of velocity and acceleration
% vectors and plots the motion of a ball given these vectors. As a
% challenge, I decided to plot this motion in three dimensions.
disp('For each vector acting on the ball, please input the Magnitude of Velocity, Bearing and Angle');

exit = 0; %I've nested each while loop for each input in a larger while loop so that the user can input multiple vectors. This is the condition variable for this loop.
n = 1; %This variable is to count the total number of vectors
inputs = struct; 
%I initiate a struct to hold the input data for each vector. A struct or
%structure array is a data structure that groups data under fields. This is
%different from a Hashmap or matrix in that each field under the struct can
%be a different data type.

while exit == 0
    %%% Vector type %%%
    % I chose to allow the user to input acceleration vectors as well as
    % velocity vectors, since they'll have an extra dimension to play with.
    
    tt = 0;
    while tt == 0
        type = input(sprintf("Type of Vector %.0f, (a)cceleration or (v)elocity: ",n),'s'); 
        % Since this input is a string, I only need to check if the input
        % is the right length(in this case length one) and check if it's an
        % allowed character
        
        % I use | instead of || and & instead of && since if the input is
        % longer than one character, we must compare a character vector
        % rather than an individual character.
        if(isscalar(type) & type == 'a' | type == 'v') %#ok<AND2,OR2> <- This is here to tell MATLAB to shut up about the operators I'm using during code inspection
            tt = 1;
            inputs.types(n) = type; %If the input is valid, I store it under the "types" field, which will be used to sort the data later.
        else
            disp('Invalid input, please enter i for an intial vector and a for an acceleration vector');
        end
    end
    
    %%% Magnitude %%%
    %Standard error check.
    mt = 0;
    while mt == 0
        magnitude = input(sprintf("Enter magnitude for Vector %.0f: ", n), 's');
        if (isnan(str2double(magnitude)) || abs(imag(str2double(magnitude))) > 0)
            disp("Invalid input, please enter a number");
        else
            inputs.magnitudes(n) = str2double(magnitude); %I store the input here under the "Magnitudes" field, which will be used to calculate components later.
            mt = 1;
        end
    end
    
    %%% Bearing %%%
    %Another standard error check, it doesn't really matter what the user
    %inputs due to how the sine function works.
    bt = 0;
    while bt == 0
        bearing = input(sprintf("Enter bearing for Vector %.0f (Horizontal angle): ", n), 's');
        if (isnan(str2double(bearing)) || abs(imag(str2double(bearing))) > 0)
            disp("Invalid input, please enter a number");
        else
            inputs.bearings(n) = str2double(bearing); % I store this input under "bearings" which will also be used to calculate components.
            bt = 1;
        end
    end
    
    %%% Angle %%%
    %Basically the same input as bearing but for a different axis.
    at = 0;
    while at == 0
        angle = input(sprintf("Enter an angle for Vector %.0f (Vertical angle): ", n), 's');
        if (isnan(str2double(angle)) || abs(imag(str2double(angle))) > 0)
            disp("Invalid input, please enter a number");
        else
            inputs.angles(n) = str2double(angle);
            at = 1;
        end
    end
    
    %%% Done? %%%
    % I ask the user if they want to input any more vectors.
    yt = 0;
    while yt == 0
        finished = input("Any more vectors? (y or n): ", 's');
        switch upper(finished) % I convert the input here to upper case so that the user can enter Y or y for yes and it would still work.
            case 'N' 
                exit = 1; % If the user says no, then the outer while loop condition is changed to break the loop.
                yt = 1;
            case 'Y'
                n = n + 1; % If the user says yes, then the vector number variable is incremented and the input process is allowed to repeat.
                yt = 1; % The while loop for this input still needs to be broken.
            otherwise
                disp("Invalid input, please enter y for Yes or n for No"); % If the input doesn't satisfy any of the cases above, it defaults to this and the loop repeats.
        end
    end
end

%% Processing

disp("Processing..."); % I let the user know the program is thinking.

velocities = struct; % I then initiate structs for the velocity vectors and acceleration vectors.
accelerations = struct;

%%% Gravity %%%
% While the user can input as many or as little accelerations as they want,
% there will always be the constant acceleration due to gravity.
accelerations(1).magnitude = 9.81;
accelerations(1).bearing = 0;
accelerations(1).angle = 270; % Structs can be vectors, which is what I do here to hold every input.
% Having it as a vector means that the user can input as many or as little
% vectors as they want and the program will be able to process it.

%%% Sorting %%%
% Here I sort the inputs into their respective groups based on the value of
% their type.

vn = 1; % These two variables keep track of how many of each vector type there are.
an = 2; % This one starts at two since gravity is already in the struct.
for i = 1:n % This for loop will increment up to the total number of vectors counted earlier.
    if inputs.types(i) == 'v'
        velocities(vn).magnitude = inputs.magnitudes(i); %For each type of input, I store its magnitude, angle and bearing in similar fields in its respective struct.
        velocities(vn).bearing = inputs.bearings(i);
        velocities(vn).angle = inputs.angles(i);
        vn = vn + 1;
    elseif inputs.types(i) == 'a'
        accelerations(an).magnitude = inputs.magnitudes(i);
        accelerations(an).bearing = inputs.bearings(i);
        accelerations(an).angle = inputs.angles(i);
        an = an + 1;
    end
end
%%% Initial Components %%%
% Here is where I calculate the x, y, and z components for each vector. I
% use the spherical co-ordinate system to calculate the components.

% Sperical co-ordinates are defined as a magnitude, an angle from the x
% axis in the x-y plane and an agle from the z axis in the x-z plane. My
% method, to make it more in line with mechanics, takes an angle from the y
% axis in the x-y plane and an angle from the x axis in the x-z plane.
% Because of this, I replaced sin with cos and vice verse to make sure the
% maths checks out.

v_Components = struct;
% I store the components of the velocities in their own struct and
% calculate each component based off of the spherical co-ordinates system.
if vn ~= 1 % Arrays in MATLAB are 1 based, so if I try to access an array at value 0, MATLAB crashes, so I check for this here.
    for i = 1:(vn-1)
        v_Components(i).x = velocities(i).magnitude * cosd(velocities(i).angle) * sind(velocities(i).bearing);
        v_Components(i).y = velocities(i).magnitude * cosd(velocities(i).angle) * cosd(velocities(i).bearing);
        v_Components(i).z = velocities(i).magnitude * sind(velocities(i).angle);
    end
end

%%% Acceleration Components %%%
% I do the same as above but with vectors of acceleration.
a_Components = struct;

for i = 1:(an-1)
    a_Components(i).x = accelerations(i).magnitude * cosd(accelerations(i).angle) * sind(accelerations(i).bearing);
    a_Components(i).y = accelerations(i).magnitude * cosd(accelerations(i).angle) * cosd(accelerations(i).bearing);
    a_Components(i).z = accelerations(i).magnitude * sind(accelerations(i).angle);
end

%% Resultant calculation

v_resultant = struct('x',0,'y',0,'z',0); % I initiate more structs to hold the resultant vectors.
a_resultant = struct('x',0,'y',0,'z',0); % I must initiate the fields here so that I can use the field in calculations without them being empty.

%%% Initial Resultant %%%
if vn ~= 1 % I check if any velocity vectors are present before trying to calculate, otherwise MATLAB would crash
    for i = 1:(vn-1)
        v_resultant.x = v_resultant.x + v_Components(i).x; % Here I sum up the values for each component into a resultant component.
        v_resultant.y = v_resultant.y + v_Components(i).y;  %I store each resultant component in its own struct.
        v_resultant.z = v_resultant.z + v_Components(i).z; % This is for the velocity vectors
    end
end

%%% Acceleration Resultant %%%`
for i = 1:(an-1)
    a_resultant.x = a_resultant.x + a_Components(i).x; % I do the same for acceleration vectors
    a_resultant.y = a_resultant.y + a_Components(i).y;
    a_resultant.z = a_resultant.z + a_Components(i).z;
end

%%% Maximum time calculation %%%
% Here I calculate the time that the ball will take to fly. I only
% need to consider z components here sinze when the z displacement is 0 the
% ball is on the ground. I do this so that I can A) only calculate
% what I need to calculate and B) calculate displacement in reasonable
% steps of time, to avoid a memory overflow for very large time frames.
T_max = 2*(-v_resultant.z/a_resultant.z);  % This equation is a rearrangement of a = (v - u)/t. Doing this for the z axis can find the time up, so  to find the time 
                                                                   % we multiply it by two
t = 0:T_max/1000:T_max; % I then create an array that displacement can be calculated with.
                                           % I chose steps of 1000 arbitrarily and it gave a favourable result, so I stuck with it.

%%% Displacement Calculation %%%
% I then calculate the displacement for each component using S = u*t + 1/2 * a * t^2;
S = struct('x',0,'y',0,'z',0); % I store the displacement in one final struct.
S.x = v_resultant.x .* t + 0.5 * a_resultant.x .* t.^2;
S.y = v_resultant.y .* t + 0.5 * a_resultant.y .* t.^2;
S.z = v_resultant.z .* t + 0.5 * a_resultant.z .* t.^2;

v_vector = sqrt(v_resultant.x^2 + v_resultant.y^2 + v_resultant.z^2); % I calculate the magnitude of the resultant velocity
a_vector = sqrt(a_resultant.x^2 + a_resultant.y^2 + a_resultant.z^2); % I calculate the magnitude of the resultant acceleration
x_vector = S.x(1:100:end); % These three lines every hundredth value from the displacement arrays to create a subset of ten values.
y_vector = S.y(1:100:end);
z_vector = S.z(1:100:end);
maxDist = sqrt(S.x(end)^2 + S.y(end)^2); % I then calculate the maximum displacement in the x-y plane, since the final z displacement will be zero.

disp("Done!");
%% Displaying
%%% Graph initialisation %%%
hold all; % This command prevents matlab from erasing previous plots when drawing new ones.
grid on; % This command makes the grid visible
view(3); % This command tells MATLAB I want to view a 3d graph.
plot3(S.x,S.y,S.z); % I then plot the motion of the ball using the displacement values calculated earlier.
rotate3d('on'); % This tells MATLAB to rotate the plot when I click and drag on the plot.
plot(S.x(end),S.y(end), 'hr'); % I plot the maximum decision and display it as a red star.
text(S.x(end),S.y(end), sprintf("Maximum distance: %.3fm", maxDist)); % I then annotate this star with the value of the final displacement.

quiver3(0,0,0,v_resultant.x/2,v_resultant.y/2,v_resultant.z/2);
% This command plots an arrow at the initial position of the ball
% showing the direction of the resultant velocity. I half the values of the
% resultant components to ensure that the arrow isn't too long relative to
% the graph.
text(v_resultant.x/2,v_resultant.y/2,v_resultant.z/2, sprintf("Initial speed: %.4fm/s", v_vector)); % I then annotate this vector with the value of the magnitude of the initial velocity.
% I intentionally chose to not draw the resultant angles since the text
% would overlap and clutter the plot.

half = ones([1 length(x_vector)]); % This command creates an array and populates it with a set of ones
half = half ./10; % I then half each of these ones to create an array of halfs.

%I do this so that I can periodically show the acceleration vector on the
%curve with a set of arrows. This is why I also created a set of ten values
%from the displacement vectors. This command will draw ten arrows emerging
%from the graph in the direction of the resultant vector.
quiver3(x_vector,y_vector,z_vector,half.*a_resultant.x,half.*a_resultant.y,half.*a_resultant.z,0,'b');
text(0,0,0,sprintf("Acceleration: %.4fm/s^2", a_vector)); % I then annotate the graph to display the magnitude of the acceleration.
pbaspect([1 1 1]); % This command locks the aspect ratio of the graph so it doesn't warp as the user rotates it.
end

function caseTwo
% For this case the user is asked to enter an initial height, a launch
% angle and a launch speed. The program then calculates the maximum
% horizontal displacement and displays it. I chose not to calculate this in
% three dimensions as I thought it unnecessary.
%% Inputs
% Standard error check, here I also check if the height is greater than
% zero, since if it's less than zero then the ball would be underground
ht = 0;
while ht == 0
    h = input("Please enter the launch height from the ground: ", 's');
    if(abs(imag(str2double(h))) > 0 || isnan(str2double(h)) || str2double(h) < 0)
        disp("Invalid input, please enter a number that is greater than 0.");
    else
        h = str2double(h);
        ht = 1;
    end
end
% Another standard error check. I check if the angle is less that 180
% degress otherwise it would launch into the ground.
at = 0;
while at == 0
    a = input("Please enter the launch angle in degrees: ", 's');
    if(isnan(str2double(a)) || str2double(a) > 180 && str2double(a) > 0 || abs(imag(str2double(a))) > 0)
        disp("Invalid input, please enter an angle that is less than 180 degrees");
    else
        a = str2double(a);
        at = 1;
    end
end
% Another Standard error check, making sure the speed is greater than zero
% since negative speed won't make sense.
ut = 0;
while ut == 0
    u = input("Please enter the launch speed in metres per second: ", 's');
    if(isnan(str2double(u)) || str2double(u) < 0 || abs(imag(str2double(u))) > 0)
        disp("Invalid input, please enter a positive value for speed: ", 's');
    else
        u = str2double(u);
        ut = 1;
    end
end

%% Processing
% The motion of the ball will be split into 2 parts: a parabola up to
% when it reaches its initial height and a period when it's falling from
% that height.
T_parabola = 2*(-u*sind(a))/-9.81; % I can work out the time for the parabola using a rearranged SUVAT equation
u_down = u*sind(a) + -9.81 * T_parabola; %  I calculate the initial speed for the time down
[t1, t2] = quadraticEquation(-4.905, u_down, h); % I then calculate the time down by solving the quadratic formula s-u*t + 0.5 * a * t ^2 = 0
if t1 > 0
    T_down = t1; % The equation will have a positive and negative answer, the positive answer is the true value for time.
else
    T_down = t2; % Here I select the positive time.
end

T_max = T_parabola + T_down; % I add up the calculated times to get the total time.

t = 0:(T_max/1000):T_max; % I create an array I can use in calculations.

S_x = u*cosd(a)*t; %Calculating horizontal displacement.
S_y = u*sind(a)*t + 0.5*-9.81*(t.^2); % Calculating vertical displacement.
max_x = max(S_x); % I store the value of the maximum displacement for later

y_max = max(S_y) + h; % I store the value of the maximum height for later
x_max = u*cosd(a) * T_parabola / 2; % I calculate the x value of the max height


%% Displaying
hold on;
grid on;

plot(S_x,S_y+h); % I plot the motion of the ball, offset by the maximum height
plot(x_max, y_max,'hr'); % I mark the maximum height the ball reached
text(x_max, y_max + (y_max/10), sprintf("Max height: %.2fm",y_max)); % and annotate it
quiver(0, 0, max_x,0,0); %I draw a vector to show the maximum displacement
quiver(0,h, u*cosd(a),u*sind(a)); %and another vector to show the initial velocity.
text(max_x/2,-(y_max/4),sprintf("Max distance: %.2fm", max_x)); % I then annotate the displacement on the graph
bw = xlim; % I then create a small matrix that holds the values of the graph limits
bar(0,h,bw(2)/5); % and create a bar at the initial height of the ball and of a width  proportional to the screen.
% This bar represents a platform that the ball could be launched from


end

function caseThree
%% Inputs
% Standard error check, also check if the time is greater than zero since
% negative time doesn't make sense
tt = 0;
while tt == 0
    t = input("Please enter the time of flight in seconds: ", 's');
    if( isnan(str2double(t)) || str2double(t)<0 || abs(imag(str2double(t))) > 0)
        disp("Invalid input, please enter a positive number for time");
    else
        t = str2double(t);
        tt = 1;
    end
end

% Standard error check.
st = 0;
while st == 0
    s = input("Please enter the displacement in metres: ", 's');
    if (isnan(str2double(s)) || abs(imag(str2double(s))) > 0)
        disp("Invalid input, please enter a number as the value for speed");
    else
        s = str2double(s);
        st = 1;
    end
end

%% Processing
T = 0:t/1000:t; % I create a matrix for time that I can use for calculations
S = 0:s/1000:s; % and a similar matrix for displacement, 
h = 0.5 * -9.81 * T.^2; %I then calculate the height that a ball would travel in this time
h = h + abs(min(h)); % And then offset it so it lands at zero

u_x = s/t; % I then calculate the initial speed
H = max(h); % And store the value for maximum height for later calculations.

%% Displaying
hold on;
grid on;
plot(S,h); % I plot the motion of the ball
bw_t = xlim;
bw = bw_t(2)/5; % and create a bar similar to the one in the previous case
bar(0-bw/2,H,bw);

q1 = quiver(0+bw/4, 0, 0, H, 0,'g'); % I then create two arrow objects to represent the height
q1.MaxHeadSize = 0.01; % I initialise them as objects so I can control how large the heads get
% as in some cases they got too large for the plot.
q2 = quiver(0+bw/4, H, 0, -H, 0,'g');
q2.MaxHeadSize = 0.01;
text(0+bw/2, H/2, sprintf("Launch height: %.2fm", H)); % I then annotate these arrows.

quiver(0+bw/2,H,u_x,0,0,'r'); % I then create another arrow to display the initial speed.
text(0+bw/2+u_x, H, sprintf("Initial speed: %.2fm/s", u_x));
end

function caseFour
% In this case the user is asked for the time of flight and the total
% horizontal displacement of the ball and the program will calculate the
% initial launch speed and the initial launch angle.
%% Inputs
% Standard error check, I also check if the time is less than or equal to
% zero, since negative time doesn't make sense and zero time means the ball
% hasn't moved (also wouldn't make sense if there's any displacement).
tt = 0;
while tt == 0
    ts = input("Please enter the time of flight in seconds: ", 's');
    t = str2double(ts); 
    if (isnan(t) || t<=0 || abs(imag(t)) > 0)
        disp("Invalid input, please enter a positive number for time that is greater than zero");
    else
        tt = 1;
    end
end
% Another standard error check. Here I only check if the displacement is
% equal to zero since I allow for negative displacement (just means firing
% to the left).
st = 0;
while st == 0
    ss = input("Please enter the displacement in metres: ", 's');
    s = str2double(ss);
    if (isnan(s) || abs(imag(s)) > 0 || s == 0)
        disp("Invalid input, please enter a number as the value for displacement that is greater than zero");
    else
        st = 1;
    end
end
%% Processing
u_x = s/t; % Calculates the x component of the initial speed
u_y = (9.81 * t)/2; %Calculates the Y component (this equation is a rearranged
% form of equations for time previously used)

% With the components of the velocity calculated, I can now calculate:
u = sqrt(u_x^2 + u_y^2); % The magnitude of the initial velocity
theta = atan(u_y/u_x); % And the angle of the initial velocity (left in radians for later use).

S_x = 0:s/1000:s; % I then create some arrays for time and displacement which I can
% use for calculations.
T = 0:t/1000:t;

S_y = u_y .* T + 0.5 * -9.81 .* T.^2; % I calculate the vertical displacement of the ball
% with the calculated initial speed.
S_max = max(S_y); % I then store the maximum displacement for later use.
%% Displaying
hold on;
grid on;
daspect([1 1 1]); % I lock the aspect ratio of the data and the graph to be square/
pbaspect([1 1 1]); % I do this due to how the arc function works.


plot(S_x, S_y); % I then plot the motion of the ball
plot(s/2, S_max, 'h'); % I mark the position of the maximum height with a star
text(s/2, S_max + (S_max / 10), sprintf("Max height: %.2fm", S_max)); % and annotate it.

q1 = quiver(0,0,u_x,u_y,0); % I then initialise an arrow object so that I can control the scaling of it
q1.MaxHeadSize = 0.1; % This arrow will represent the initial velocity.
text(u_x, u_y, sprintf("Initial speed: %.2fm/s", u)); % I then annotate the arrow

ang([0 0],u/2, [0 theta], 'k'); % Here I draw an arc from the X axis to the arrow to represent the angle.
text(u_x/2, u_y/2, sprintf("Launch angle: %.2f degrees", rad2deg(theta))); % I then annotate inside this arc, converting the angle to degrees due to convention.
 axis manual; % I then lock the axis so that they don't change.
plot([-500,500],[0 0],'k-'); % I then plot a line to represent the ground.
% This line would cause the plot to zoom out a long way, which is why I
% locked the axis.
end

%% Extra functions
% This function applies the quadratic formula and returns both values for X
function [x1, x2] = quadraticEquation(a,b,c)
x1 = (-b + sqrt(b^2 - 4 * a * c))/(2*a);
x2 = (-b - sqrt(b^2 - 4 * a * c))/(2*a);
end

% This function draws an arc given a centre point, radius and angle, as
% well as support for graph styling. I copied this from the MATLAB file
% exchange since MATLAB doesn't have a function for drawing arc sections
% (even though other languages do :thinking:)

% I found this function here:
% https://uk.mathworks.com/matlabcentral/fileexchange/21940-plot-the-angle-arc-for-some-angle/
function h = ang(centre,radius,span,style)
% Credit goes to Husam Aldahiyat, October, 2008.
theta = linspace(span(1),span(2),100);
rho = ones(1,100) * radius;
[x,y] = pol2cart(theta,rho);
x = x + centre(1);
y = y + centre(2);
h = plot(x,y,style);
end
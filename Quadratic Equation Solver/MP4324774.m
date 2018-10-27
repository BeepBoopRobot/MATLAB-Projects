%The code below was written by Benjamin Langford
% The function of this script is to calculate and display the roots of a
% quadratic using the quadratic formula,

clc;    %clc removes any previous commands or messages from the command window.
% It is good practice to do this so that any new error messages aren't
% confused with previous ones and to see if the program is outputting
% correctly
clear; %This command clears any previously initialised variables from memory
%  We use this command to ensure that memory is free, since variables
%  persist between scripts, thus if I made any changes after making a large
%  matrix, the matrix would still be in memory and could cause a memory
%  overflow.
close all; % This command closes all figures that have graphs on them.
hold on; % This command makes it so that everything plotted is persistent, meaning 
% that when plot() is called more than once, it doesn't overwrite the previous plot

% TO DO: 
% Make it so intersection points are displayed
% IF TIME ALLOWS:
% Upgrade it so it finds the intersection between a quadratic and a line
% Upgrade itso it finds the intersection between two quadratics
% Make it so that all the lines are displayed on one graph
% Make it so that the quadratic solver parses a quadratic equation instead
% of asking for a b and c separately
% DEADLINE: NOVEMBER 7TH
%

disp('Please enter a, b, and c for a quadratic from the form ax^2 + bx + c, as well as the bounds for the graph');
%This command prints to the command window only. I do this to prompt the
%user for input.


at = 0; % This variable is the condition variable for the while loop below.
while at == 0 %The while loop is set to repeat until a condition is met. In this case I set it to a value of the condition variable above
    a = input('a: ', 's'); %the input() function defaults to a number input and as such has its own error checking
    %     I intend to bypass this by having the program print out a more
    %     user-friendly message to the user. To do this, I must first take the
    %     input as a string rather than a number.
    if isnan(str2double(a)) % I then convert the string to a double. If the variable is anything but a number, this function returns NaN
        % NaN thankfully has a function that can check for it and return a status
        % value which I can use in an if statement.
        disp('Invalid input, expected a number'); % This is the message to prompt the user to enter a numeric input.
        % The if statement skips the else part and the loop starts again
    else
        at = 1; % if the input is numeric, then the condition variable is set to break the loop and continue
        a = str2double(a); % I then convert the string input to a double and store it. I use str2double to include fractions
    end
end
% I then do the same for the other variables.
bt = 0;
while bt == 0
    b = input('b: ', 's');
    if isnan(str2double(b))
        disp('Invalid input, expected a number');
    else
        bt = 1;
        b = str2double(b);
    end
end
ct = 0;
while ct == 0
    c = input('c: ', 's');
    if isnan(str2double(c))
        disp('Invalid input, expected a number');
    else
        ct = 1;
        c = str2double(c);
    end
end

lt = 0;
while lt == 0
    L = input('Lower Bound of graph (default -5): ', 's');
    if isnan(str2double(L))
        L = -5;
        lt = 1;
    else
        lt = 1;
        L = str2double(L);
    end
end

rt = 0;
while rt == 0
    R = input('Upper bound of graph (default 5): ', 's');
    if isnan(str2double(R))
        R = 5;
        rt = 1;
    else
        rt = 1;
        R = str2double(R);
    end
end
fprintf('You have entered %d, %d, and %d \n', a,b,c);

x = L:0.01:R;
y = a*x.^2 + b.*x + c;

plot(x,y);

% I then check the discriminant of the equation to check the number of
% roots the equation will have. This means that I can create a function for
% each case.
discriminant = b^2 - 4 * a * c; % Now that I have calculated this, I can use it in later calculations.

if discriminant > 0 %If the discriminant is greater than 0, the the equation will have two roots
    [x1,x2] = twoRoots(a,b,discriminant); % I call the function that handles two roots and 
%     pass the values for a, b and the discriminant to it. I don't need to
%     pass c since it was already used in the calculation for the
%     discriminant

    plot(x1,0, 'p');
    plot(x2,0, 'p');
    fprintf('x1 = %.3f, x2 = %.3f \n', x1, x2); % I then print out 
%     the values for the roots
elseif discriminant == 0 % If the discriminant equals 0, then there will be one root.
    disp(oneRoot(a,b)); % I only need to pass a and b since the discriminant is 0
else
    disp('Equation will have no x axis intersects'); % I let the user know that there will be no interesect on the graph
    [xi1, xi2] = twoRoots(a,b,discriminant); % I use the same function in the first case since there will still be 2 roots,
%     however they will not be real roots, that is to say the function will
%     return complex numbers rather than floats
    fprintf('x1 = %.3f%+.3fi ,  x2 = %.3f+%.3fi \n',real(xi1),imag(xi1),real(xi2),imag(xi2));
%     Since fprintf doesn't have a format for complex numbers, I need to
%     split them into their real and imaginary parts and concatenate them
%     separately
end

function [x1, x2] = twoRoots(a,b,d)
x1 = (-b + sqrt(d))/(2*a);
x2 = (-b - sqrt(d))/(2*a);
end

function x = oneRoot(a,b)
x = -b/2*a;
end
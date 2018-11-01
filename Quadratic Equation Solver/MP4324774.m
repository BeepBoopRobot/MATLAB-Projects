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
a = [0 0 0 0 0 0 0 0 0 0];
b = [0 0 0 0 0 0 0 0 0 0];
c = [0 0 0 0 0 0 0 0 0 0];
maxt = 0;
while maxt == 0
    max = input('How many equations do you want to plot? (Max 10) ', 's');
    if isnan(str2double(max))
        disp('Invalid input, expected a number');
    elseif str2double(max) >10
        disp('Input out of range, maximum input is 10');
    else
        maxt = 1;
        max = str2double(max);
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
for i = 1:max
    fprintf('Please enter a, b, and c for quadratic %.0f \n', i);
    %This command prints to the command window only. I do this to prompt the
    %user for input. I type %.0f to tell MATLAB I want the value i to be
    %substituted there. The .0 tells MATLAB the precision of the number in
    %decimal places. Since it will be an integer, I do not want any
    %decimalplaces, hence the .0. I also use put a \n in there to indicate
    %a new line.
    
    
    at = 0; % This variable is the condition variable for the while loop below.
    while at == 0 %The while loop is set to repeat until a condition is met. In this case I set it to a value of the condition variable above
        as = input('a: ', 's'); %the input() function defaults to a number input and as such has its own error checking
        %     I intend to bypass this by having the program print out a more
        %     user-friendly message to the user. To do this, I must first take the
        %     input as a string rather than a number.
        if isnan(str2double(as)) % I then convert the string to a double. If the variable is anything but a number, this function returns NaN
            % NaN thankfully has a function that can check for it and return a status
            % value which I can use in an if statement.
            disp('Invalid input, expected a number'); % This is the message to prompt the user to enter a numeric input.
            % The if statement skips the else part and the loop starts again
        else
            at = 1; % if the input is numeric, then the condition variable is set to break the loop and continue
            a(i) = str2double(as); % I then convert the string input to a double and store it. I use str2double to include fractions
        end
    end
    % I then do the same for the other variables.
    bt = 0;
    while bt == 0
        bs = input('b: ', 's');
        if isnan(str2double(bs))
            disp('Invalid input, expected a number');
        else
            bt = 1;
            b(i) = str2double(bs);
        end
    end
    ct = 0;
    while ct == 0
        cs = input('c: ', 's');
        if isnan(str2double(cs))
            disp('Invalid input, expected a number');
        else
            ct = 1;
            c(i) = str2double(cs);
        end
    end
    
    fprintf('You have entered %d, %d, and %d \n', a(i),b(i),c(i));
end
hold on; % This command tells matlab not to overwrite the previous graph when I plot a new one
% (for some reason matlab counts plotting individual points as separate
% graphs worth clearing the graph for
grid on; % This makes the grid visible on the graph.

title(sprintf('Graph of %.0f equations',max));
% This command sets the title for the graph. I use sprintf to allow me to
% concatenate the values for a, b, and c into the title.


plot([L R],[0 0],'k-');
plot([0 0], [-500 500], 'k-');
axis([L R 10*L 10*R]);

hold on;
for i = 1:max
    x = L:0.01:R;
    if ~(a(i) == 0)
        y = a(i)*x.^2 + b(i)*x + c(i);
        plot(x,y,'-.');
    else
        y = b(i).*x + c(i);
        plot(x,y,'-.');
    end
    
    % I then check the discriminant of the equation to check the number of
    % roots the equation will have. This means that I can create a function for
    % each case.
    discriminant = b(i)^2 - 4 * a(i) * c(i); % Now that I have calculated this, I can use it in later calculations.
    if discriminant > 0 %If the discriminant is greater than 0, the the equation will have two roots
        [x1, x2] = twoRoots(a(i),b(i),discriminant); % I call the function that handles two roots and
        %     pass the values for a, b and the discriminant to it. I don't need to
        %     pass c since it was already used in the calculation for the
        %     discriminant.
        plot(x1,0, 'pg'); %With the roots calculated, I then plot them on the graph.
        plot(x2,0, 'pg'); %I use 'p' as a constructor to tell the function to plot these points as green stars
        
        fprintf('Equation %.0f: x1 = %.3f, x2 = %.3f \n',i , x1, x2); % I then print out
        %     the values for the roots
    elseif discriminant == 0 % If the discriminant equals 0, then there will be one root.
        x = oneRoot(a(i),b(i)); % I only need to pass a and b since the discriminant is 0
        plot(x,0,'pr'); %This command tells MATLAB to plot a red star at the x intercept
        fprintf('Equation %.0f: x = %.3f \n',i,x);
    else
        fprintf('Equation %.0f will have no x axis intersects \n', i); % I let the user know that there will be no interesect on the graph
        [xi1, xi2] = twoRoots(a(i),b(i),discriminant); % I use the same function in the first case since there will still be 2 roots,
        %     however they will not be real roots, that is to say the function will
        %     return complex numbers rather than floats
        fprintf('Equation %.0f: x1 = %.3f%+.3fi ,  x2 = %.3f-%.3fi \n',i,real(xi1),imag(xi1),real(xi2),imag(xi2));
        %     Since fprintf doesn't have a format for complex numbers, I need to
        %     split them into their real and imaginary parts and concatenate them
        %     separately. The first root will have a positive imaginary
        %     part while the second will have negative part. imag() only
        %     returns the magnitude so I have to put the signs of each part
        %     in manually.
    end
end

function [x1, x2] = twoRoots(a,b,d) % This function executes the
%     quadratic equation. It returns each root in a matrix, which is how
%     MATLAB handles multiple inputs
x1 = (-b + sqrt(d))/(2*a); % Since the quadtratic equation has a ± in it, I can split it into 2 equations.
x2 = (-b - sqrt(d))/(2*a);
end

function x = oneRoot(a,b)
x = (-b)/(2*a); % Since the discriminant will equal zero, then the equation will have
% just one root, the two equations defined by the quadrataic formula will
% have the same value so I do not need to pass the discriminant into the
% equation.
end

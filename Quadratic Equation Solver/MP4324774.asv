%The code below was written by Benjamin Langford

%We first start by using the commands clc and clear
clc;    %clc removes any previous commands or messages from the command window.
% It is good practice to do this so that any new error messages aren't
% confused with previous ones and to see if the program is outputting
% correctly
clear; %This command clears any previously initialised variables from memory
%  We use this command to ensure that memory is free, since variables
%  persist between scripts, thus if I made any changes after making a large
%  matrix, the matrix would still be in memory and could cause a memory
%  overflow.

% TO DO: Take input
% declare and validate variables
% allow user to input, a,b,c as well as domain and range for graph
% Program part to calculate roots, allowing for: no roots, one root, two
% roots
% Make it so intersection points are displayed
% IF TIME ALLOWS:
% Upgrade it so it finds the intersection between a quadratic and a line
% Upgrade itso it finds the intersection between two quadratics
% Make it so that all the lines are displayed on one graph
% Make it so that the quadratic solver parses a quadratic equation instead
% of asking for a b and c separately
% DEADLINE: NOVEMBER 7TH
%

disp('Please enter a, b, and c for a quadratic from the form ax^2 + bx + c');
%This command prints to the command window only. I do this to prompt the
%user for input.

%The three lines below initialises a variable for a, b, and c and takes
%their value from the user. The input() command reads a the value from the
%command window after the prompt when the user presses enter
at = 0;
while at == 0
    a = input('a: ', 's');
    if isnan(str2double(a)) && ~isscalar(a)
        disp('invalid');
    else
        at = 1;
        a = str2double(a);
    end
end
bt = 0;
while bt == 0
    b = input('b: ', 's');
    if isnan(str2double(b)) && ~isscalar(b)
        disp('invalid');
    else
        bt = 1;
        b = str2double(b);
    end
end
ct = 0;
while ct == 0
    c = input('c: ', 's');
    if isnan(str2double(c)) && ~isscalar(c)
        disp('invalid');
    else
        ct = 1;
        c = str2double(c);
    end
end

fprintf('You have entered %d, %d, and %d \n', a,b,c);
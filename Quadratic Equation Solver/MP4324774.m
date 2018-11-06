%% QUADRATIC EQUATION SOLVER
%% Written by Benjamin Langford
% This program is designed to take the values for a, b, and c from a
% quadratic equation in the form y = ax^2 + bx + c (or straight lines in the
% form y = bx + c
%% Housekeeping
clc;    %This command clears the console. We to this to make the command window neater and easier yo read
clear;      %This command deletes any variable stored in the workspace. We do this to A) there are no conflicts
% with previously initialised variables and B) so that there's enough space in memory to run the program.
close all;  %This command closes all previously opened graphs. I do this to
%% Initial values

<<<<<<< HEAD
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
=======
disp('For a quadratic, enter a,b, and c for a quadratic in the form y = ax^2 + bx + c\n.');
disp('For straight lines leave a blank and enter b and c for equations in the form y = bx + c.');
% I display a prompt to let the user know what values to enter. I use two
% disp() commands to display it on to different lines in the console, using
% "\n" as a line break.
a = [0 0 0 0 0 0 0 0 0 0]; %These three matrices are what will allow me to take multiple inputs for each variable.
b = [0 0 0 0 0 0 0 0 0 0]; %While I could expand the matrixes to be much larger, I chose to limit them to 10 values
c = [0 0 0 0 0 0 0 0 0 0]; %since any more graphs than that would be hard to look at.

%% Graph Inputs

% I surround each input prompt in a while loop. I do this so that the
% program prompts the user for an input until it gets an input it can
% accept.
maxt = 0;   % I initialise this variable to act as the condition for the while loop below
while maxt == 0                                                                                              %The while loop is set to run while the condition variable is zero
    max = input('How many equations do you want to plot? (Max 10) ', 's');      % I take the user input in as a string, since it means the user
                                                                                                                         % can enter whatever they want the input() won't throw an error
    if isnan(str2double(max))   %I then use the str2double() function to convert the string into a double variable.
        disp('Invalid input, expected a number');  % If the input string is anything but a number, the str2double() function will return NaN.
                                                                        %Thankfully there is a function that can check for this and I use this with an if statement.
                                                                        %I then prompt the user if the input is wrong.
                                                                        
    elseif str2double(max) >10  %Since in this loop I am asking the user for how many graphs they want to plot, the value cannot be bigger than 10 so I put this line in to check it.
        disp('Input out of range, maximum input is 10'); %  I then prompt the user if the input is too large.
    else
        maxt = 1;       %If the input is valid, then I set the  condition variable to one, which will stop the while loop.
        max = str2double(max);  % I then assign the variable its double value from the string.
    end
end   
 % I do something similar with this function, however I make it so that
 % there is a default value if the user does not want to specify bounds for
 % the graph.
    lt = 0;
    while lt == 0
        L = input('Lower Bound of graph (default -5): ', 's');
        if isnan(str2double(L))
            L = -5; %Instead of prompting the user for another input due to an empty input field, I instead just assign the variable to its default value.
            % I use -5 and 5  as the default bounds since these were the default bounds
            % on previous graphing calculators I have used.
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
    %% User inputs
    for i = 1:max   % I initialise a for loop that is set to count up to the maximum value specified earlier in the program.
                            % The variable that I count with, i, is used to indicate the position in
                            % the matrixes for a,b, and c when assigning values, among other things.
        fprintf('Please enter a, b, and c for quadratic %.0f \n', i); %I use fprintf in place of disp() here so that I can insert a variable into the message.
         % With %.0f, I specify the format for the variable to be outputted and its
         % precision. I then specify the variable to be substituted after the
         % message.
        
        at = 0;
        while at == 0 % Similar for loop the the ones before.
            as = input('a: ', 's'); %This one differs in that I use a separate variable for the
                                            %string input, since a is defined as a matrix.
            if isnan(str2double(as)) % I use this separate variable in all of the checks.
                disp('Invalid input, expected a number');
            else
                at = 1;
                a(i) = str2double(as); %When I reassign the string, I assign it to a position within the matrix, defined by i.
            end
        end
        % I do the same for b and c
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
        
        fprintf('You have entered %d, %d, and %d \n', a(i),b(i),c(i)); % I then print out to the user what they have entered.
                                                                                                  % Originally this was a debug step, however I decided to leave it in since
                                                                                                  % it verifies to the user that their values have been saved correctly.
>>>>>>> 5e797811bbe837dd90b31f3852d2704d5c5948bb
    end
    %% Graph initialisation
    % Before we can plot the graphs, we need to set some things up first.
    hold on;    % If this was set to off, then every time we use the plot() function, it would delete the previous plot and redraw it,
                     % which is undesirable when we want to plot multiple graphs.
    grid on;    % I also initialise the grid, since it makes the scale of the graph clearer to the viewer.
    title(sprintf('Graph of %.0f equations',max)); % I set the title to tell the user how many equations are (or should be) on screen.
    plot([-500 500],[0 0],'k-'); % This and the command below draw the x and y axes for x = 0 and y = 0. I set the bounds so large
                                              % so that if the user wished to move the graph around then the axes won't
                                              % disappear.
    plot([0 0], [-500 500], 'k-'); % I specify 'k-' to tell the MATLAB to plot each line as a black straight line. For each new line on a graph MATLAB automatically assigns a new colour unless told otherwise.
    axis([L R 10*L 10*R]);            % MATLAB automatically scales the graph to fit everything. This isn't ideal if we're going to have potentially
                                                 % large graphs plotted so I tell MATLAB to only view within these bounds
                                                 % (the left and right bound previously specified and a ratio of these
                                                 % bounds for the top and bottom.    
    %% Graph plotting
    % I use this for loop to plot each individual graphs. I do this in a for
    % loop since I can use the counting variable to access the matrices at
    % desired indices.
    for i = 1:max
        x = L:0.01:R; %I set up a matrix to hold all of the values for x that will be used in the graph, between the bounds specified by the user.
  
            y = a(i)*x.^2 + b(i).*x + c(i); % I then calculate the y values, using the values at position i in each matrix a, b, and c.
            plot(x,y,'-.'); %With each matrix calculated, I can now plot the equation on to a graph.
            % The .- tells MATLAB to plot the graph as a dot-dash line. I do not
            % specify a colour, since, as mentioned before, matlab will assign a
            % different colour to each new graph, which is what I want to happen to
            % distinguish the graphs.
        
            %% Roots Calculations
            % Here is where we actually calculate the roots of each graph. I could have
            % just done x = (-b + sqrt(b(i)^2+4*a(i)*c(i))/2(a(i)), but I chose to do
            % some extra processing to make the results nicer.
            
            d = b(i)^2 - 4 * a(i) * c(i);        % I first calculate the discriminant of each formula.
                                                                             %Since the discriminant defines how many roots each graph will have, based
                                                                             %on 3 conditions, I can use it to form an if statement.
                                                                             
            if d > 0 %The first condition is where the discriminant is greater than zero, where it will have 2 real roots.
                
                [x1, x2] = twoRoots(a(i),b(i),d); %I pass the values for a and b at position i and the discriminant into a function that returns the roots.
                plot(x1,0, 'pg'); %I then plot each root on the star, using the parameter 'pg' to tell MATLAB to plot a green star.
                plot(x2,0, 'pg');
                fprintf('Equation %.0f: x1 = %.3f, x2 = %.3f \n',i , x1, x2);  % I then print the values of each roots in the console.
                                                                                                            % While I could plot the values of each root on the graph, if the values
                                                                                                            % are too close together, I feel that they will be difficult to read, so
                                                                                                            % printing them in the console makes it easier for the user to see the
                                                                                                            % roots of each equation.
                
            elseif d == 0 % The second condition is when the discriminant is zero, where it will only have 1 real root.
                
                x = oneRoot(a(i),b(i)); % I then pass the value of a and b at position i into a function I will explain below
                plot(x,0,'pr'); % I plot the root, defining it as a red star to differentiate it from other roots.
                fprintf('Equation %.0f: x = %.3f \n',i,x); % I then print the value of the root into the console.
                                                                              % Note that I specify .0f for the value of i and .3f for the value of the
                                                                              % root. This is because i will always be an integer, due to how the loop is
                                                                              % set up, and because x could be a fraction.
                                                                              
            else % The last condition (which I don't need to explicitly specify due to the nature of if statements) is when the discriminant is less than zero, where it will have no real roots.
                
                fprintf('Equation %.0f will have no x axis intersects \n', i); % I first let the user know that this is the case.
                [xi1, xi2] = twoRoots(a(i),b(i),d); % I still pass the values for a and b at position i and the discriminant into the same function, since the graph will still have roots,
                                                                   % except that they will not be real roots.
                fprintf('Equation %.0f: x1 = %.3f%+.3fi ,  x2 = %.3f-%.3fi \n',i,real(xi1),imag(xi1),real(xi2),imag(xi2));
                % I also print out the roots in the same way, although I have to substitute
                % them into fprintf slightly differenly. I know from the quadratic equation
                % that  one root will have a positive imaginary part and one will have a
                % negative imaginary part, so I can specify these in the fprintf function
                % manualle. fprintf doesn't have a mechanism for printing imaginary
                % numbers, so instead I have to substitute the real and imaginary parts
                % separately.
            end
    end
<<<<<<< HEAD
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
=======
    %% Calculations
    function [x1, x2] = twoRoots(a,b,d)
    x1 = (-b + sqrt(d))/(2*a); % The quadratic equation is defined with a ±in it. Since I can't code this directly, I split it into two equations and return the values as a matrix.
    x2 = (-b - sqrt(d))/(2*a); % I do not need to pass c into the function (even though the equation needs it) since I already used it in the calculation for the discriminant.
    end
    
    function x = oneRoot(a,b)
    x = (-b)/(2*a); % In this equation, I do not even pass the discriminant, since it will be 0 in either case. Because of this I can also just return a single value.
                            % (Note: One thing that I personally do not like about MATLAB is the
                            % implicit variable designations. MATLAB defaults to a number since it's a
                            % mathematical language, however I personally would like the ability to
                            % define what the output of a function would be. I guess that's the Java
                            % inside me speaking)
    end
    
    %% END
>>>>>>> 5e797811bbe837dd90b31f3852d2704d5c5948bb

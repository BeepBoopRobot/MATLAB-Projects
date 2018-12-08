clc;
clear;
a = input("pls fucking work", 's');
b = str2double(a);
disp(imag(b));
disp(abs(imag(b)));
disp(abs(imag(b)) > 0);
if(abs(imag(b)) > 0)
    disp("fuck");
end
clc;
clear;
a = struct;
b = struct;
for i = 1:10
    a(i).a = 2*i;
    a(i).b = sqrt(i);
end
for i = 1:10
    b.a(i) = 2*i;
    b.b(i) = sqrt(i);
end
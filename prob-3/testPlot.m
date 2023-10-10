clearvars
close all

x = linspace(0,pi,301);
y1 = sin(x);
y2 = cos(x);


plot(x,y1,'r-')
hold on
plot(x,-y2,'-b')
axis equal
grid on
aXis1 = gca;
hold off
f2 = figure
aXis2 = copy(aXis1, f2)
hold on
plot(x,y2,'g-')
hold off
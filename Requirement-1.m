%Requirement 1:
t = -pi:0.1:pi 
a.xh = 4*sin(t).^5+5;
a.yh = 3*cos(t)-1.7*cos(2*t)-cos(3*t)+1;
plot(a.xh,a.yh);
fill(a.xh, a.yh, 'r');
%Requirement 2:

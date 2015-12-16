%% 
% PERCEPTION V1 and Direction Selectivity Assignment
%
% 1)
% a)
deltaT = 1; % ms
duration = 1000; % ms
t = [0:deltaT:duration-deltaT];
x = zeros(size(t));
x(100)=1;
y1=0;
y1Save = zeros(size(t));
tau = 25; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT/tau) * (-y1 + x(tt));
    y1 = y1 + deltaY1;
    y1Save(tt) = y1;
end
plot(y1Save)
hold on;
% ? plot(exp(-t/tau))
tau = 50; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT/tau) * (-y1 + x(tt));
    y1 = y1 + deltaY1;
    y1Save(tt) = y1;
end
plot(y1Save)

b)
x(100:1000)=1;
y1=0;
y1Save = zeros(size(t));
tau = 25; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT/tau) * (-y1 + x(tt));
    y1 = y1 + deltaY1;
    y1Save(tt) = y1;
end
plot(y1Save)

for tt = 1:length(t)
    leftS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, 0, 1, 2*pi*tt/hz).*mkGaussian(length(x), (sigma*pixPerDeg)^2);
    upS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, pi/2, 1, 2*pi*tt/hz);
    rightS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, pi, 1, 2*pi*tt/hz);
    downS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, 3*pi/2, 1, 2*pi*tt/hz);
    showIm(leftS(:,:,tt));
    pause
end
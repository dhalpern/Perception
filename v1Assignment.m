addpath('../motionTutorial/')
addpath('../motionTutorial/pyrTools/')
addpath('../motionTutorial/pyrTools/MEX/')

%% 
% PERCEPTION V1 and Direction Selectivity Assignment
%
% 1)
% a)
deltaT = 1; % ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
x = zeros(size(t));
x(100)=1;
u = zeros(size(t));
u(100:1000)=1;
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
plot((deltaT/tau) * exp(-(t-99)/tau) .* u)
hold off;
tau = 50; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT/tau) * (-y1 + x(tt));
    y1 = y1 + deltaY1;
    y1Save(tt) = y1;
end
plot(y1Save)
hold on; 
plot((deltaT/tau) * exp(-(t-99)/tau) .* u)
%%
b)
y1=0;
y1Save = zeros(size(t));
tau = 25; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT/tau) * (-y1 + u(tt));
    y1 = y1 + deltaY1;
    y1Save(tt) = y1;
end
plot(y1Save)
hold on
plot((1 - exp(-(t-99)/tau)) .* u)
%%
% Since the form of the differential
% c)
x = sin(pi/32 * t);
y1=0;
y1Save = zeros(size(t));
tau = 25; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT/tau) * (-y1 + x(tt));
    y1 = y1 + deltaY1;
    y1Save(tt) = y1;
end
plot(y1Save)
hold on
plot(x)

%%
% 2
x = zeros(size(t));
x(100)=1;
y1Save = zeros(size(t));
y2Save = zeros(size(t));
y3Save = zeros(size(t));
y4Save = zeros(size(t));
y5Save = zeros(size(t));
y6Save = zeros(size(t));
y7Save = zeros(size(t));
y1=0;
y2=0;
y3=0;
y4=0;
y5=0;
y6=0;
y7=0;
tau = 25; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT / tau) * (-y1 + x(tt));
    y1 = y1 + deltaY1;
    deltaY2 = (deltaT / tau) * (-y2 + y1);
    y2 = y2 + deltaY2;
    deltaY3 = (deltaT / tau) * (-y3 + y2);
    y3 = y3 + deltaY3;
    deltaY4 = (deltaT / tau) * (-y4 + y3);
    y4 = y4 + deltaY4;
    deltaY5 = (deltaT / tau) * (-y5 + y4);
    y5 = y5 + deltaY5;
    deltaY6 = (deltaT / tau) * (-y6 + y5);
    y6 = y6 + deltaY6;
    deltaY7 = (deltaT / tau) * (-y7 + y6);
    y7 = y7 + deltaY7;
    y1Save(tt) = y1;
    y2Save(tt) = y2;
    y3Save(tt) = y3;
    y4Save(tt) = y4;
    y5Save(tt) = y5;
    y6Save(tt) = y6;
    y7Save(tt) = y7;
end
f1 = y3Save-y5Save;
f2 = y5Save-y7Save;
plot(f1)
plot(f2)

%%
% 3)
%deltaX = 1/120;
%deltaT = 1; % ms
deltaX = 1/20;
deltaT = 10; % ms
duration = 1000; % ms
Xrange = 2;
t = 0:deltaT:duration-deltaT;
x = -Xrange:deltaX:Xrange-deltaX;
y = -Xrange:deltaX:Xrange-deltaX;
sigma = .1;
cycPerDeg = 4;
pixPerDeg = 1/deltaX;
hz = 8 * deltaX; % 8Hz
odd = -sin(x*cycPerDeg*2*pi) .* normpdf(x,0,sigma);
even = cos(x*cycPerDeg*2*pi) .* normpdf(x,0,sigma);
vertOdd = conv2(odd, normpdf(x,0,sigma)');
vertEven = conv2(even, normpdf(x,0,sigma)');
horOdd = conv2(normpdf(x,0,sigma), odd');
horEven = conv2(normpdf(x,0,sigma), even');

sXT = [zeros(length(t)-length(x), length(x));eye(length(x))]; % x-t image
s = repmat(sXT, [1 1 length(x)]);
s = permute(s, [2 3 1]);
y1Save = zeros(length(x), length(y), length(t));
y2Save = zeros(length(x), length(y), length(t));
y3Save = zeros(length(x), length(y), length(t));
y4Save = zeros(length(x), length(y), length(t));
y5Save = zeros(length(x), length(y), length(t));
y6Save = zeros(length(x), length(y), length(t));
y7Save = zeros(length(x), length(y), length(t));
vertOddSlow = zeros(length(x), length(y), length(t));
vertEvenSlow = zeros(length(x), length(y), length(t));
vertOddFast = zeros(length(x), length(y), length(t));
vertEvenFast = zeros(length(x), length(y), length(t));
horOddSlow = zeros(length(x), length(y), length(t));
horEvenSlow = zeros(length(x), length(y), length(t));
horOddFast = zeros(length(x), length(y), length(t));
horEvenFast = zeros(length(x), length(y), length(t));
y1=0;
y2=0;
y3=0;
y4=0;
y5=0;
y6=0;
y7=0;
tau = 25; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT / tau) * (-y1 + s(:,:,tt));
    y1 = y1 + deltaY1;
    deltaY2 = (deltaT / tau) * (-y2 + y1);
    y2 = y2 + deltaY2;
    deltaY3 = (deltaT / tau) * (-y3 + y2);
    y3 = y3 + deltaY3;
    deltaY4 = (deltaT / tau) * (-y4 + y3);
    y4 = y4 + deltaY4;
    deltaY5 = (deltaT / tau) * (-y5 + y4);
    y5 = y5 + deltaY5;
    deltaY6 = (deltaT / tau) * (-y6 + y5);
    y6 = y6 + deltaY6;
    deltaY7 = (deltaT / tau) * (-y7 + y6);
    y7 = y7 + deltaY7;
    f1 = y3 - y5;
    f2 = y5 - y7;
    horOddFast(:,:,tt) = conv2(horOdd, f1, 'same');
    horEvenFast(:,:,tt) = conv2(horEven, f1, 'same');
    horOddSlow(:,:,tt) = conv2(horOdd, f2, 'same');
    horEvenSlow(:,:,tt) = conv2(horEven, f2, 'same');
    vertOddSlow(:,:,tt) = conv2(vertOdd, f2, 'same');
    vertEvenSlow(:,:,tt) = conv2(vertEven, f2, 'same');
    vertOddFast(:,:,tt) = conv2(vertOdd, f1, 'same');
    vertEvenFast(:,:,tt) = conv2(vertEven, f1, 'same');
    
    y1Save(:,:,tt) = y1;
    y2Save(:,:,tt) = y2;
    y3Save(:,:,tt) = y3;
    y4Save(:,:,tt) = y4;
    y5Save(:,:,tt) = y5;
    y6Save(:,:,tt) = y6;
    y7Save(:,:,tt) = y7;
end
imagesc(permute(vertOddFast(:,  40, :), [3 1 2]))
imagesc(permute(vertOddSlow(:,  40, :), [3 1 2]))
imagesc(permute(vertEvenFast(:,  40, :), [3 1 2]))
imagesc(permute(vertEvenSlow(:,  40, :), [3 1 2]))

% b) 
leftEven = horOddFast + horEvenSlow;
leftOdd = -horOddSlow + horEvenFast;
rightEven = -horOddFast + horEvenSlow;
rightOdd = horOddFast + horEvenFast;
upEven = vertOddFast + vertEvenSlow;
upOdd = -vertOddSlow + vertEvenFast;
downEven = -vertOddFast + vertEvenSlow;
downOdd = vertOddFast + vertEvenSlow;
imagesc(permute(rightEven(:,  40, :), [3 1 2]))
imagesc(permute(rightOdd(:,  40, :), [3 1 2]))
imagesc(permute(leftEven(:,  40, :), [3 1 2]))
imagesc(permute(leftOdd(:,  40, :), [3 1 2]))

% c)
leftEnergy = leftEven.^2  + leftOdd.^2;
rightEnergy = rightEven.^2 + rightOdd.^2;
upEnergy = upEven.^2  + upOdd.^2;
downEnergy = downEven.^2 + downOdd.^2;
showIm(permute(leftEnergy(:,  40, :), [3 1 2]))
showIm(permute(rightEnergy(:,  40, :), [3 1 2]))

% d)
deltaX = 1/20;
deltaT = 10; % ms
duration = 1000; % ms
Xrange = 2;
t = 0:deltaT:duration-deltaT;
x = -Xrange:deltaX:Xrange-deltaX;
y = -Xrange:deltaX:Xrange-deltaX;
sigma = .1;
cycPerDeg = 4;
pixPerDeg = 1/deltaX;
hz = 8 * deltaX; % 8Hz
leftS = zeros(length(x), length(y), length(t));
upS = zeros(length(x), length(y), length(t));
rightS = zeros(length(x), length(y), length(t));
downS = zeros(length(x), length(y), length(t));

odd = -sin(x*cycPerDeg*2*pi) .* normpdf(x,0,sigma);
even = cos(x*cycPerDeg*2*pi) .* normpdf(x,0,sigma);
vertOdd = conv2(odd, normpdf(x,0,sigma)');
vertEven = conv2(even, normpdf(x,0,sigma)');
horOdd = conv2(normpdf(x,0,sigma), odd');
horEven = conv2(normpdf(x,0,sigma), even');
for tt = 1:length(t)
    leftS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, 0, 1, 2*pi*tt/hz).*mkGaussian(length(x), (sigma*pixPerDeg)^2);
    upS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, pi/2, 1, 2*pi*tt/hz);
    rightS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, pi, 1, 2*pi*tt/hz);
    downS(:,:,tt) = mkSine(length(x), pixPerDeg/cycPerDeg, 3*pi/2, 1, 2*pi*tt/hz);
    showIm(leftS(:,:,tt));
    pause
end
s = leftS;
vertOddSlow = zeros(length(x), length(y), length(t));
vertEvenSlow = zeros(length(x), length(y), length(t));
vertOddFast = zeros(length(x), length(y), length(t));
vertEvenFast = zeros(length(x), length(y), length(t));
horOddSlow = zeros(length(x), length(y), length(t));
horEvenSlow = zeros(length(x), length(y), length(t));
horOddFast = zeros(length(x), length(y), length(t));
horEvenFast = zeros(length(x), length(y), length(t));
y1=0;
y2=0;
y3=0;
y4=0;
y5=0;
y6=0;
y7=0;
tau = 25; % ms
for tt = 1:length(t)
    deltaY1 = (deltaT / tau) * (-y1 + s(:,:,tt));
    y1 = y1 + deltaY1;
    deltaY2 = (deltaT / tau) * (-y2 + y1);
    y2 = y2 + deltaY2;
    deltaY3 = (deltaT / tau) * (-y3 + y2);
    y3 = y3 + deltaY3;
    deltaY4 = (deltaT / tau) * (-y4 + y3);
    y4 = y4 + deltaY4;
    deltaY5 = (deltaT / tau) * (-y5 + y4);
    y5 = y5 + deltaY5;
    deltaY6 = (deltaT / tau) * (-y6 + y5);
    y6 = y6 + deltaY6;
    deltaY7 = (deltaT / tau) * (-y7 + y6);
    y7 = y7 + deltaY7;
    f1 = y3 - y5;
    f2 = y5 - y7;
    horOddFast(:,:,tt) = conv2(horOdd, f1, 'same');
    horEvenFast(:,:,tt) = conv2(horEven, f1, 'same');
    horOddSlow(:,:,tt) = conv2(horOdd, f2, 'same');
    horEvenSlow(:,:,tt) = conv2(horEven, f2, 'same');
    vertOddSlow(:,:,tt) = conv2(vertOdd, f2, 'same');
    vertEvenSlow(:,:,tt) = conv2(vertEven, f2, 'same');
    vertOddFast(:,:,tt) = conv2(vertOdd, f1, 'same');
    vertEvenFast(:,:,tt) = conv2(vertEven, f1, 'same');
end
leftEven = horOddFast + horEvenSlow;
leftOdd = -horOddSlow + horEvenFast;
rightEven = -horOddFast + horEvenSlow;
rightOdd = horOddFast + horEvenFast;
upEven = vertOddFast + vertEvenSlow;
upOdd = -vertOddSlow + vertEvenFast;
downEven = -vertOddFast + vertEvenSlow;
downOdd = vertOddFast + vertEvenSlow;
leftEnergy = leftEven.^2  + leftOdd.^2;
rightEnergy = rightEven.^2 + rightOdd.^2;
upEnergy = upEven.^2  + upOdd.^2;
downEnergy = downEven.^2 + downOdd.^2;

impulse = zeros(length(t), length(x));
impulse(50, 20) = 1;
impulse = permute(repmat(impulse, [1 1 length(y)]), [2 3 1]);

s = impulse;
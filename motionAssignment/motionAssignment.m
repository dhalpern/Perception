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
odd = -sin(x*8*pi) .* normpdf(x,0,sigma);
even = cos(x*8*pi) .* normpdf(x,0,sigma);
%vertOdd = conv2(odd, normpdf(x,0,sigma)');
%vertEven = conv2(even, normpdf(x,0,sigma)');
%horOdd = conv2(normpdf(x,0,sigma), odd');
%horEven = conv2(normpdf(x,0,sigma), even');

vertOdd  = mkSine(length(x), pixPerDeg/cycPerDeg, 0, 1, 0).*mkGaussian(length(x), (sigma*pixPerDeg)^2);
horOdd  = mkSine(length(x), pixPerDeg/cycPerDeg, pi/2, 1, 0).*mkGaussian(length(x), (sigma*pixPerDeg)^2);
vertEven = mkSine(length(x), pixPerDeg/cycPerDeg, 0, 1, pi/2).*mkGaussian(length(x), (sigma*pixPerDeg)^2);
horEven = mkSine(length(x), pixPerDeg/cycPerDeg, pi/2,1, pi/2).*mkGaussian(length(x), (sigma*pixPerDeg)^2);


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

showIm(permute(vertOddFast(:,  1, :), [3 1 2]))
showIm(permute(vertOddSlow(:,  1, :), [3 1 2]))
showIm(permute(vertEvenFast(:,  1, :), [3 1 2]))
showIm(permute(vertEvenSlow(:,  1, :), [3 1 2]))
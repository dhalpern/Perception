addpath('../motionTutorial/')
addpath('../motionTutorial/pyrTools/')
addpath('../motionTutorial/pyrTools/MEX/')

close all;

%% 1: Recursive Temporal Filter
deltaT = 1; % ms
duration = 1000; % ms
t = [0:deltaT:duration-deltaT];

%% 1.a Impulse Response
x = zeros(size(t));
x(100)=1;


y1=0;
y1Save = zeros(size(t));
tau1 = 25; % ms
for tt = 1:length(t)
  deltaY1 = (deltaT/tau1) * (-y1 + x(tt));
  y1 = y1 + deltaY1;
  y1Save(tt) = y1;
end
yf1 = exp(-(t)./tau1);

%for tau2 = 1:5:800
y2=0;
y2Save = zeros(size(t));
tau2 = 50; % ms
for tt = 1:length(t)
  deltaY2 = (deltaT/tau2) * (-y2 + x(tt));
  y2 = y2 + deltaY2;
  y2Save(tt) = y2;
end
yf2 = exp(-(t)./tau2);

figure(1);hold on;
plot(t, y1Save, 'r--');
plot(t, yf1, 'r-');
plot(t, y2Save, 'b--');
plot(t, yf2, 'b-');
ylim([0 0.05]);
xlabel('Time (ms)'); ylabel('Response');title('Impulse Response')
legend(sprintf('tau1(red): %d', tau1), 'e^{t/tau1}', sprintf('tau2(blue): %d', tau2), 'e^{t/tau2}');
hold off
%pause
%close
%end

%% 1.b Step Response
x = zeros(size(t));
x(200:end)=1;


y1=0;
y1Save = zeros(size(t));
tau1 = 25; % ms
for tt = 1:length(t)
  deltaY1 = (deltaT/tau1) * (-y1 + x(tt));
  y1 = y1 + deltaY1;
  y1Save(tt) = y1;
end
yf1 = 1-exp(-t./tau1);


y2=0;
y2Save = zeros(size(t));
tau2 = 50; % ms
for tt = 1:length(t)
  deltaY2 = (deltaT/tau2) * (-y2 + x(tt));
  y2 = y2 + deltaY2;
  y2Save(tt) = y2;
end
yf2 = 1-exp(-(t)./tau2);

figure(2);hold on;
plot(t, y1Save, 'r--');
plot(t, yf1, 'r-');
plot(t, y2Save, 'b--');
plot(t, yf2, 'b-');
% ylim([0 0.05]);
xlabel('Time (ms)'); ylabel('Response');title('Step Response')
legend(sprintf('tau1(red): %d', tau1), 'e^{t/tau1}', sprintf('tau2(blue): %d', tau2), 'e^{t/tau2}');
hold off


%% 1.c Frequency Response

%for tau1=1:50:5000
w1 = 2^0/length(t);
w2 = 2^1/length(t);
w3 = 2^4/length(t);
w4 = 2^8/length(t);

x1 = sin(2*pi*w1*t);
x2 = sin(2*pi*w2*t);
x3 = sin(3*pi*w3*t);
x4 = sin(3*pi*w4*t);

y1=0;
y2=0;
y3=0;
y4=0;

y1Save = zeros(size(t));
y2Save = zeros(size(t));
y3Save = zeros(size(t));
y4Save = zeros(size(t));

tau1 = 25; % ms
tau = tau1;


for tt = 1:length(t)
  deltaY1 = (deltaT/tau) * (-y1 + x1(tt));
  deltaY2 = (deltaT/tau) * (-y2 + x2(tt));
  deltaY3 = (deltaT/tau) * (-y3 + x3(tt));
  deltaY4 = (deltaT/tau) * (-y4 + x4(tt));

  y1 = y1 + deltaY1;
  y2 = y2 + deltaY2;
  y3 = y3 + deltaY3;
  y4 = y4 + deltaY4;

  y1Save(tt) = y1;
  y2Save(tt) = y2;
  y3Save(tt) = y3;
  y4Save(tt) = y4;
end
Y = [y1Save; y2Save; y3Save; y4Save]';

x = zeros(size(t));
x(1) = 1;
y0=0;
y0Save = zeros(size(t));
tau1 = 25; % ms
for tt = 1:length(t)
  deltaY0 = (deltaT/tau1) * (-y0 + x(tt));
  y0 = y0 + deltaY0;
  y0Save(tt) = y0;
end
impulseResponse = y0Save;

disp('Freqeuncies (Hz):')
disp(1000*[w1 w2 w3 w4]);


disp('Simulated amplitudes:')
disp(max(Y(200:end,:)));

disp('Predicted amplitudes:')
% amplitude[y^(w)] = amplitude[d^d(w)] x amplitude[x^(w)]I
max(abs(fft(impulseResponse).*fft(y1Save)));

% this is wrong! produces exact same result as this one since
% fft(impulseReponse) does not change much of the peak, hence max()
%amp = 2*max(abs(repmat(fft(impulseResponse)',1,4).*[fft(y1Save)' fft(y2Save)' fft(y3Save)' fft(y4Save)']))./length(t)

amp = 2*max(abs([fft(y1Save)' fft(y2Save)' fft(y3Save)' fft(y4Save)']))./length(t);

disp('Phase:')
phase = pi/2 + 0


figure(3);
subplot(2,1,1);hold on;
plot(t, y1Save, 'r--');
plot(t, y2Save, 'b--');
plot(t, y3Save, 'm--');
plot(t, y4Save, 'k--');
hold off
xlabel('Time (ms)');  ylabel('Response');title(sprintf('Frequency Response (tau: %d)', tau1));
legend(cellstr([num2str(1e3*[w1 w2 w3 w4]') repmat(' Hz',4,1)]))


subplot(2,1,2); hold on
freq = 2*[w1 w2 w3 w4]/length(t);
semilogx(freq, max(Y(200:end,:)), 'ro');
semilogx(freq, amp, 'r--');
xlabel('Log Frequency (Hz)');  ylabel('Response');
legend('Simulated', 'Predicted')
hold off

%hold off
%pause
%close
%end

%% 2: Exponential Low Pass Filter
clear;

deltaT = 1; % ms
duration = 1000; % ms
t = [0:deltaT:duration-deltaT];
tau = 25; % ms

x = zeros(size(t));
x(1)=1;
y1 = 0;
y2 = 0;
y3 = 0;
y4 = 0;
y5 = 0;
y6 = 0;
y7 = 0;


y1Save = zeros(size(t));
y2Save = zeros(size(t));
y3Save = zeros(size(t));
y4Save = zeros(size(t));
y5Save = zeros(size(t));
y6Save = zeros(size(t));
y7Save = zeros(size(t));
f1Save = zeros(size(t));
f2Save = zeros(size(t));
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

  f1 = y3-y5;
  f2 = y5-y7;

  f1Save(tt) = f1;
  f2Save(tt) = f2;

end

figure(4);hold on;
plot(t, y1Save);
plot(t, y2Save);
plot(t, y3Save);
plot(t, y4Save);
plot(t, y5Save);
plot(t, y6Save);
plot(t, y7Save);
xlabel('Time (ms)');  ylabel('Response');title('A Cascade of Exponential Filters')
hold off

figure(5);hold on;
plot(t, f1Save);
plot(t, f2Save);
xlabel('Time (ms)');  ylabel('Response');title('A Cascade of Exponential Low-Pass Filters')
hold off


%% 3: Space-time filters and motion energy filters
clear;

% deltaT = 1; % ms
% deltaX = 1/200; %deg, 30 sec of arc, cone spacing in the retina

deltaT = 10; % ms
deltaX = 0.05; %deg, 30 sec of arc, cone spacing in the retina


duration = 1000; % ms
t = [0:deltaT:duration-deltaT];
tau = 25; % ms

% 3-D stimulus array (x,y,t)
rangeX = -2:deltaX:2;
XT = [zeros(length(t)-length(rangeX), length(rangeX));eye(length(rangeX))]; % x-t image
x = repmat(XT, [1 1 length(rangeX)]);
x = permute(x, [2 3 1]);
figure(101);
subplot(1,2,1);showIm(permute(x(:,1,:), [3 1 2]));xlabel('x-t slice of the input (y=1)');
subplot(1,2,2);showIm(permute(x(:,:,60), [2 1 3]));xlabel('x-y slice of the input (t=60)');
% print -dpng -r300 Q3InputStim.png

% make gabors
% [xx, yy] = meshgrid(rangeX, rangeX);
sigma = 0.1;
w = 4;
pixPerDeg = length(rangeX)/range(rangeX);

verOdd  = mkSine(length(rangeX), pixPerDeg/w, 0,    1, 0).*mkGaussian(length(rangeX),    (sigma*pixPerDeg)^2);
horOdd  = mkSine(length(rangeX), pixPerDeg/w, pi/2, 1, 0).*mkGaussian(length(rangeX),    (sigma*pixPerDeg)^2);
verEven = mkSine(length(rangeX), pixPerDeg/w, 0,    1, pi/2).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);
horEven = mkSine(length(rangeX), pixPerDeg/w, pi/2, 1, pi/2).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);

figure(6);
subplot(2,2,1);showIm(horOdd);ylabel('Odd Phase');set(get(gca,'Ylabel'),'Visible','on')
subplot(2,2,2);showIm(verOdd);
subplot(2,2,3);showIm(horEven);ylabel('Even Phase');xlabel('Horizontal');set(get(gca,'Ylabel'),'Visible','on')
subplot(2,2,4);showIm(verEven);xlabel('Vertical');
% print -dpng -r300 Q3GaborFilters.png



% init values
y1 = zeros(length(rangeX));
y2 = zeros(length(rangeX));
y3 = zeros(length(rangeX));
y4 = zeros(length(rangeX));
y5 = zeros(length(rangeX));
y6 = zeros(length(rangeX));
y7 = zeros(length(rangeX));

% containers
horOddFast  = zeros(size(x));
horEvenFast = zeros(size(x));
verOddFast  = zeros(size(x));
verEvenFast = zeros(size(x));
horOddSlow  = zeros(size(x));
horEvenSlow = zeros(size(x));
verOddSlow  = zeros(size(x));
verEvenSlow = zeros(size(x));

%f1Save = zeros(size(x));
%f2Save = zeros(size(x));

for tt = 1:length(t)
  deltaY1 = (deltaT / tau) * (-y1 + x(:,:,tt));
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

  %y1Save(tt) = y1;
  %y2Save(tt) = y2;
  %y3Save(tt) = y3;
  %y4Save(tt) = y4;
  %y5Save(tt) = y5;
  %y6Save(tt) = y6;
  %y7Save(tt) = y7;

  f1 = y3-y5;
  f2 = y5-y7;

  %f1Save(:,:,tt) = f1;
  %f2Save(:,:,tt) = f2;

  horOddFast(:,  :, tt) = conv2(horOdd,  f1, 'same');
  horEvenFast(:, :, tt) = conv2(horEven, f1, 'same');
  verOddFast(:,  :, tt) = conv2(verOdd,  f1, 'same');
  verEvenFast(:, :, tt) = conv2(verEven, f1, 'same');
  horOddSlow(:,  :, tt) = conv2(horOdd,  f2, 'same');
  horEvenSlow(:, :, tt) = conv2(horEven, f2, 'same');
  verOddSlow(:,  :, tt) = conv2(verOdd,  f2, 'same');
  verEvenSlow(:, :, tt) = conv2(verEven, f2, 'same');

end



%% 3.a x-t slices of the impulse response for 4 filters that prefer vertical
figure(7);
subplot(2, 2, 1);showIm(permute(verOddFast(:,  1, :), [3 1 2]));ylabel('Odd Phase');set(get(gca,'Ylabel'), 'Visible', 'on')
subplot(2, 2, 2);showIm(permute(verOddSlow(:,  1, :), [3 1 2]));
subplot(2, 2, 3);showIm(permute(verEvenFast(:, 1, :), [3 1 2]));ylabel('Even Phase');xlabel('Fast');set(get(gca, 'Ylabel'), 'Visible', 'on')
subplot(2, 2, 4);showIm(permute(verEvenSlow(:, 1, :), [3 1 2]));xlabel('Slow');
% print -dpng -r300 Q3avertical.png

% figure;
% subplot(2, 2, 1);showIm(permute(horOddFast(:,  1, :), [3 1 2]));ylabel('Odd Phase');set(get(gca,'Ylabel'), 'Visible', 'on')
% subplot(2, 2, 2);showIm(permute(horOddSlow(:,  1, :), [3 1 2]));
% subplot(2, 2, 3);showIm(permute(horEvenFast(:, 1, :), [3 1 2]));ylabel('Even Phase');xlabel('Fast');set(get(gca, 'Ylabel'), 'Visible', 'on')
% subplot(2, 2, 4);showIm(permute(horEvenSlow(:, 1, :), [3 1 2]));xlabel('Slow');



%% 3.b x-t slices of the impulse response for 2 left and 2 right
leftEven  = horOddFast  + horEvenSlow;
leftOdd   = -horOddSlow + horEvenFast;
rightEven = -horOddFast + horEvenSlow;
rightOdd  = horOddSlow  + horEvenFast;
upEven    = verOddFast  + verEvenSlow;
upOdd     = -verOddSlow + verEvenFast;
downEven  = -verOddFast + verEvenSlow;
downOdd   = verOddSlow  + verEvenFast;

figure(8);
subplot(2, 2, 1);showIm(permute(leftOdd(:,   1, :), [3 1 2]));ylabel('Odd Phase');set(get(gca,'Ylabel'), 'Visible', 'on')
subplot(2, 2, 2);showIm(permute(rightOdd(:,  1, :), [3 1 2]));
subplot(2, 2, 3);showIm(permute(leftEven(:, 1, :), [3 1 2]));ylabel('Even Phase');xlabel('Left');set(get(gca, 'Ylabel'), 'Visible', 'on')
subplot(2, 2, 4);showIm(permute(rightEven(:, 1, :), [3 1 2]));xlabel('Right');
% print -dpng -r300 Q3bLeftRight.png


%% 3.c x-t slices of the impulse response of the energy response
leftEnergy  = leftEven.^2  + leftOdd.^2;
rightEnergy = rightEven.^2 + rightOdd.^2;
upEnergy    = upEven.^2    + upOdd.^2;
downEnergy  = downEven.^2  + downOdd.^2;

figure(9);
subplot(1, 4, 1);showIm(permute(leftEnergy(:,  1, :), [3 1 2]));xlabel('Left');
subplot(1, 4, 2);showIm(permute(rightEnergy(:, 1, :), [3 1 2]));xlabel('Right');
subplot(1, 4, 3);showIm(permute(upEnergy(:,    1, :), [3 1 2]));xlabel('Up');
subplot(1, 4, 4);showIm(permute(downEnergy(:,  1, :), [3 1 2]));xlabel('Down');
% print -dpng -r300 Q3cEnergy.png


%% 3.d drifting sinosoid stimuli

% create drifting gabors
v = 80; % 8Hz
leftX = zeros(size(x));
for tt = 1:length(t)
  leftX(:,:,tt) = mkSine(length(rangeX), pixPerDeg/w, 0, 1, -2*pi*tt/v).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);
%   showIm(leftX(:,:,tt));
%   pause
end
downX = rot90(leftX);
rightX = rot90(downX);
upX = rot90(rightX);

% run the Q3.a-c simulation above for each gabor
stim = {leftX, rightX, upX, downX};
for iStim = 1:numel(stim)
  x = stim{iStim};


  % init values
  y1 = zeros(length(rangeX));
  y2 = zeros(length(rangeX));
  y3 = zeros(length(rangeX));
  y4 = zeros(length(rangeX));
  y5 = zeros(length(rangeX));
  y6 = zeros(length(rangeX));
  y7 = zeros(length(rangeX));

  % containers
  horOddFast  = zeros(size(x));
  horEvenFast = zeros(size(x));
  verOddFast  = zeros(size(x));
  verEvenFast = zeros(size(x));
  horOddSlow  = zeros(size(x));
  horEvenSlow = zeros(size(x));
  verOddSlow  = zeros(size(x));
  verEvenSlow = zeros(size(x));

  for tt = 1:length(t)
    deltaY1 = (deltaT / tau) * (-y1 + x(:,:,tt));
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


    f1 = y3-y5;
    f2 = y5-y7;

    horOddFast(:,  :, tt) = conv2(horOdd,  f1, 'valid');
    horEvenFast(:, :, tt) = conv2(horEven, f1, 'valid');
    verOddFast(:,  :, tt) = conv2(verOdd,  f1, 'valid');
    verEvenFast(:, :, tt) = conv2(verEven, f1, 'valid');
    horOddSlow(:,  :, tt) = conv2(horOdd,  f2, 'valid');
    horEvenSlow(:, :, tt) = conv2(horEven, f2, 'valid');
    verOddSlow(:,  :, tt) = conv2(verOdd,  f2, 'valid');
    verEvenSlow(:, :, tt) = conv2(verEven, f2, 'valid');

  end

  leftEven  = horOddFast  + horEvenSlow;
  leftOdd   = -horOddSlow + horEvenFast;
  rightEven = -horOddFast + horEvenSlow;
  rightOdd  = horOddSlow  + horEvenFast;
  upEven    = verOddFast  + verEvenSlow;
  upOdd     = -verOddSlow + verEvenFast;
  downEven  = -verOddFast + verEvenSlow;
  downOdd   = verOddSlow  + verEvenFast;

  leftEnergy  = leftEven.^2  + leftOdd.^2;
  rightEnergy = rightEven.^2 + rightOdd.^2;
  upEnergy    = upEven.^2    + upOdd.^2;
  downEnergy  = downEven.^2  + downOdd.^2;

  figure(9 + iStim);
  subplot(4, 1, 1)
  hold on
  plot(squeeze(leftOdd(   floor(end/2), floor(end/2), :)));
  plot(squeeze(leftEven(  floor(end/2), floor(end/2), :)));
  plot(squeeze(leftEnergy(floor(end/2), floor(end/2), :)));
  hold off
  xlabel('Left');
  legend('Odd', 'Even', 'Energy');

  subplot(4, 1, 2)
  hold on
  plot(squeeze(rightOdd(   floor(end/2), floor(end/2), :)));
  plot(squeeze(rightEven(  floor(end/2), floor(end/2), :)));
  plot(squeeze(rightEnergy(floor(end/2), floor(end/2), :)));
  hold off
  xlabel('Right');
  legend('Odd', 'Even', 'Energy');

  subplot(4, 1, 3)
  hold on
  plot(squeeze(upOdd(   floor(end/2), floor(end/2), :)));
  plot(squeeze(upEven(  floor(end/2), floor(end/2), :)));
  plot(squeeze(upEnergy(floor(end/2), floor(end/2), :)));
  hold off
  xlabel('Up');
  legend('Odd', 'Even', 'Energy');

  subplot(4, 1, 4)
  hold on
  plot(squeeze(downOdd(   floor(end/2), floor(end/2), :)));
  plot(squeeze(downEven(  floor(end/2), floor(end/2), :)));
  plot(squeeze(downEnergy(floor(end/2), floor(end/2), :)));
  hold off
  xlabel('Down');
  legend('Odd', 'Even', 'Energy');
  % print -dpng -r300 Q3cEnergy.png


end

%% 4.a Normalization, contrast
clear

deltaContrast = 0.1;
contrasts = 10.^([-1:deltaContrast:0]);

deltaT = 10; % ms
deltaX = 0.05; %deg, 30 sec of arc, cone spacing in the retina
% deltaT = 1; % ms
% deltaX = 1/200; %deg, 30 sec of arc, cone spacing in the retina

duration = 1000; % ms
t = [0:deltaT:duration-deltaT];
tau = 25; % ms

% ignore first 300 ms
ignoreTimeBlock = 300;
validT = [floor((ignoreTimeBlock/duration) * length(t)):length(t)];

% 3-D stimulus array (x,y,t)
rangeX = -2:deltaX:2;
x = zeros(length(rangeX), length(rangeX), length(t));

% make gabors
sigma = 0.1;
w = 4;
pixPerDeg = length(rangeX)/range(rangeX);

verOdd  = mkSine(length(rangeX), pixPerDeg/w, 0,    1, 0).*mkGaussian(length(rangeX),    (sigma*pixPerDeg)^2);
horOdd  = mkSine(length(rangeX), pixPerDeg/w, pi/2, 1, 0).*mkGaussian(length(rangeX),    (sigma*pixPerDeg)^2);
verEven = mkSine(length(rangeX), pixPerDeg/w, 0,    1, pi/2).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);
horEven = mkSine(length(rangeX), pixPerDeg/w, pi/2, 1, pi/2).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);


% create drifting gabors at different contrasts
v = 80; % 8Hz
stim = cell(1, length(contrasts));
for iContrast = 1:length(contrasts)
rightX = zeros(size(x));
for tt = 1:length(t)
  rightX(:,:,tt) = mkSine(length(rangeX), pixPerDeg/w, 0, contrasts(iContrast), -2*pi*tt/v).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);
  %showIm(rightX(:,:,tt));
  %pause
end
stim{iContrast} = rightX;
end


% run the Q3.a-c simulation above for each gabor
centralLRUDEnergyNorm = NaN(length(stim), 4);
for iStim = 1:numel(stim)
  x = stim{iStim};


  % init values
  y1 = zeros(length(rangeX));
  y2 = zeros(length(rangeX));
  y3 = zeros(length(rangeX));
  y4 = zeros(length(rangeX));
  y5 = zeros(length(rangeX));
  y6 = zeros(length(rangeX));
  y7 = zeros(length(rangeX));

  % containers
  horOddFast  = zeros(size(x));
  horEvenFast = zeros(size(x));
  verOddFast  = zeros(size(x));
  verEvenFast = zeros(size(x));
  horOddSlow  = zeros(size(x));
  horEvenSlow = zeros(size(x));
  verOddSlow  = zeros(size(x));
  verEvenSlow = zeros(size(x));

  for tt = 1:length(t)
    deltaY1 = (deltaT / tau) * (-y1 + x(:,:,tt));
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


    f1 = y3-y5;
    f2 = y5-y7;

    horOddFast(:,  :, tt) = conv2(horOdd,  f1, 'valid');
    horEvenFast(:, :, tt) = conv2(horEven, f1, 'valid');
    verOddFast(:,  :, tt) = conv2(verOdd,  f1, 'valid');
    verEvenFast(:, :, tt) = conv2(verEven, f1, 'valid');
    horOddSlow(:,  :, tt) = conv2(horOdd,  f2, 'valid');
    horEvenSlow(:, :, tt) = conv2(horEven, f2, 'valid');
    verOddSlow(:,  :, tt) = conv2(verOdd,  f2, 'valid');
    verEvenSlow(:, :, tt) = conv2(verEven, f2, 'valid');

  end

  leftEven  = horOddFast  + horEvenSlow;
  leftOdd   = -horOddSlow + horEvenFast;
  rightEven = -horOddFast + horEvenSlow;
  rightOdd  = horOddSlow  + horEvenFast;
  upEven    = verOddFast  + verEvenSlow;
  upOdd     = -verOddSlow + verEvenFast;
  downEven  = -verOddFast + verEvenSlow;
  downOdd   = verOddSlow  + verEvenFast;

  leftEnergy  = leftEven.^2  + leftOdd.^2;
  rightEnergy = rightEven.^2 + rightOdd.^2;
  upEnergy    = upEven.^2    + upOdd.^2;
  downEnergy  = downEven.^2  + downOdd.^2;

  % Normalization
  sumEnergy = leftEnergy + rightEnergy + upEnergy + downEnergy;
  sigma = std([leftEnergy(:); rightEnergy(:); upEnergy(:); downEnergy(:)]);

  leftEnergyNorm  = leftEnergy ./ (sumEnergy + sigma^2);
  rightEnergyNorm = rightEnergy./ (sumEnergy + sigma^2);
  upEnergyNorm    = upEnergy   ./ (sumEnergy + sigma^2);
  downEnergyNorm  = downEnergy ./ (sumEnergy + sigma^2);


  leftEnergyNorm  = mean(leftEnergyNorm,  3);
  rightEnergyNorm = mean(rightEnergyNorm, 3);
  upEnergyNorm    = mean(upEnergyNorm,    3);
  downEnergyNorm  = mean(downEnergyNorm,  3);

  centralLeftEnergyNorm  = leftEnergyNorm( floor(end/2), floor(end/2));
  centralRightEnergyNorm = rightEnergyNorm(floor(end/2), floor(end/2));
  centralUpEnergyNorm    = upEnergyNorm(   floor(end/2), floor(end/2));
  centralDownEnergyNorm  = downEnergyNorm( floor(end/2), floor(end/2));

  centralLRUDEnergyNorm(iStim, :) = [ ...
    centralLeftEnergyNorm, ...
    centralRightEnergyNorm, ...
    centralUpEnergyNorm, ...
    centralDownEnergyNorm];

end

figure(14);
hold on
semilogx(contrasts, centralLRUDEnergyNorm);
hold off
xlabel('Grating contrast (log)');
ylabel('Normalized energy of the central neuron');
legend('Left', 'Right', 'Up', 'Down');
% print -dpng -r300 Q4aContrast.png



%% 4.b Normalization, cross-orientation

% run the same simulation as above, adding new orientation
clear

deltaContrast = 0.1;
contrast50 = 1-log10(2);
contrasts = 10.^([-1:deltaContrast:contrast50]);

deltaT = 10; % ms
deltaX = 0.05; %deg, 30 sec of arc, cone spacing in the retina
% deltaT = 1; % ms
% deltaX = 1/200; %deg, 30 sec of arc, cone spacing in the retina

duration = 1000; % ms
t = [0:deltaT:duration-deltaT];
tau = 25; % ms

% ignore first 300 ms
ignoreTimeBlock = 300;
validT = [floor((ignoreTimeBlock/duration) * length(t)):length(t)];

% 3-D stimulus array (x,y,t)
rangeX = -2:deltaX:2;
x = zeros(length(rangeX), length(rangeX), length(t));

% make gabors
sigma = 0.1;
w = 4;
pixPerDeg = length(rangeX)/range(rangeX);

verOdd  = mkSine(length(rangeX), pixPerDeg/w, 0,    1, 0).*mkGaussian(length(rangeX),    (sigma*pixPerDeg)^2);
horOdd  = mkSine(length(rangeX), pixPerDeg/w, pi/2, 1, 0).*mkGaussian(length(rangeX),    (sigma*pixPerDeg)^2);
verEven = mkSine(length(rangeX), pixPerDeg/w, 0,    1, pi/2).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);
horEven = mkSine(length(rangeX), pixPerDeg/w, pi/2, 1, pi/2).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);


% create drifting gabors at different contrasts
v = 80; % 8Hz
stim = cell(1, 2*length(contrasts)); % twice as long, half for cross-oriented
for iContrast = 1:length(contrasts) % only generate for one direction
rightX = zeros(size(x));
for tt = 1:length(t)
  rightX(:,:,tt) = mkSine(length(rangeX), pixPerDeg/w, 0, contrasts(iContrast), -2*pi*tt/v).*mkGaussian(length(rangeX), (sigma*pixPerDeg)^2);
  %showIm(rightX(:,:,tt));
  %pause
end
stim{iContrast} = rightX;
end

% last gabor is at 50% contrast
upX = rot90(rightX);
% fill in the rest
for iContrast = 1:length(contrasts) % cross orientation
  crossX = stim{iContrast} + upX; % add upward motion to the previous ones
  %for tt=1:length(t)
  %showIm(crossX(:,:,tt));
  %pause
  %end
  stim{iContrast+length(contrasts)} = crossX;
end


% run the Q3.a-c simulation above for each gabor
centralLRUDEnergyNorm = NaN(length(stim), 4);
for iStim = 1:numel(stim)
  x = stim{iStim};


  % init values
  y1 = zeros(length(rangeX));
  y2 = zeros(length(rangeX));
  y3 = zeros(length(rangeX));
  y4 = zeros(length(rangeX));
  y5 = zeros(length(rangeX));
  y6 = zeros(length(rangeX));
  y7 = zeros(length(rangeX));

  % containers
  horOddFast  = zeros(size(x));
  horEvenFast = zeros(size(x));
  verOddFast  = zeros(size(x));
  verEvenFast = zeros(size(x));
  horOddSlow  = zeros(size(x));
  horEvenSlow = zeros(size(x));
  verOddSlow  = zeros(size(x));
  verEvenSlow = zeros(size(x));

  for tt = 1:length(t)
    deltaY1 = (deltaT / tau) * (-y1 + x(:,:,tt));
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


    f1 = y3-y5;
    f2 = y5-y7;

    horOddFast(:,  :, tt) = conv2(horOdd,  f1, 'valid');
    horEvenFast(:, :, tt) = conv2(horEven, f1, 'valid');
    verOddFast(:,  :, tt) = conv2(verOdd,  f1, 'valid');
    verEvenFast(:, :, tt) = conv2(verEven, f1, 'valid');
    horOddSlow(:,  :, tt) = conv2(horOdd,  f2, 'valid');
    horEvenSlow(:, :, tt) = conv2(horEven, f2, 'valid');
    verOddSlow(:,  :, tt) = conv2(verOdd,  f2, 'valid');
    verEvenSlow(:, :, tt) = conv2(verEven, f2, 'valid');

  end

  leftEven  = horOddFast  + horEvenSlow;
  leftOdd   = -horOddSlow + horEvenFast;
  rightEven = -horOddFast + horEvenSlow;
  rightOdd  = horOddSlow  + horEvenFast;
  upEven    = verOddFast  + verEvenSlow;
  upOdd     = -verOddSlow + verEvenFast;
  downEven  = -verOddFast + verEvenSlow;
  downOdd   = verOddSlow  + verEvenFast;

  leftEnergy  = leftEven.^2  + leftOdd.^2;
  rightEnergy = rightEven.^2 + rightOdd.^2;
  upEnergy    = upEven.^2    + upOdd.^2;
  downEnergy  = downEven.^2  + downOdd.^2;

  % Normalization
  sumEnergy = leftEnergy + rightEnergy + upEnergy + downEnergy;
  sigma = std([leftEnergy(:); rightEnergy(:); upEnergy(:); downEnergy(:)]);

  leftEnergyNorm  = leftEnergy ./ (sumEnergy + sigma^2);
  rightEnergyNorm = rightEnergy./ (sumEnergy + sigma^2);
  upEnergyNorm    = upEnergy   ./ (sumEnergy + sigma^2);
  downEnergyNorm  = downEnergy ./ (sumEnergy + sigma^2);


  leftEnergyNorm  = mean(leftEnergyNorm,  3);
  rightEnergyNorm = mean(rightEnergyNorm, 3);
  upEnergyNorm    = mean(upEnergyNorm,    3);
  downEnergyNorm  = mean(downEnergyNorm,  3);

  centralLeftEnergyNorm  = leftEnergyNorm( floor(end/2), floor(end/2));
  centralRightEnergyNorm = rightEnergyNorm(floor(end/2), floor(end/2));
  centralUpEnergyNorm    = upEnergyNorm(   floor(end/2), floor(end/2));
  centralDownEnergyNorm  = downEnergyNorm( floor(end/2), floor(end/2));

  centralLRUDEnergyNorm(iStim, :) = [ ...
    centralLeftEnergyNorm, ...
    centralRightEnergyNorm, ...
    centralUpEnergyNorm, ...
    centralDownEnergyNorm];

end

% [single orientation, cross orientation]
centralLRUDEnergyNorm = [centralLRUDEnergyNorm(1:end/2,:) centralLRUDEnergyNorm(1+end/2:end,:)];

figure(15);
subplot(4, 1, 1)
hold on
semilogx(contrasts, centralLRUDEnergyNorm(:,[1 4+1]));
hold off
title('Normalized energy of the central neuron');
%xlabel('Grating contrast (log)');
ylabel('Left');
legend('rightward grating', 'rightward+upward grating');

subplot(4, 1, 2)
hold on
semilogx(contrasts, centralLRUDEnergyNorm(:,[2 4+2]));
hold off
%title('Normalized energy of the central neuron');
%xlabel('Grating contrast (log)');
ylabel('Right');
%legend('rightward grating', 'rightward+upward grating');

subplot(4, 1, 3)
hold on
semilogx(contrasts, centralLRUDEnergyNorm(:,[3 4+3]));
hold off
%title('Normalized energy of the central neuron');
%xlabel('Grating contrast (log)');
ylabel('Up');
%legend('rightward grating', 'rightward+upward grating');

subplot(4, 1, 4)
hold on
semilogx(contrasts, centralLRUDEnergyNorm(:,[4 4+4]));
hold off
%title('Normalized energy of the central neuron');
xlabel('Grating contrast (log)');
ylabel('Down');
%legend('rightward grating', 'rightward+upward grating');

% print -dpng -r300 Q4bCross.png
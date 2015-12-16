function [leftresp, rightresp, upresp, downresp] = temporalresponses(inputstimulus, dT, dX)

tau = 25;
spatialspan = 2;
timelength = 1000;
space = -spatialspan:dX:spatialspan;
nX = length(space);
nT = timelength/dT;
nY = 7;

% ====== making gabors ======
[X,Y]=meshgrid(space,space);
freq = 4; % cycles per degree. 

% gaussian window
sigma = 0.1; % SD of gaussian window
gauss1d = normpdf(space,0,sigma);
gauss2d = gauss1d'*gauss1d;
gauss2d = gauss2d/sum(gauss2d(:));

% four gabors
evenhorz = cos(2*pi*freq*Y).*gauss2d;
evenvert = cos(2*pi*freq*X).*gauss2d;
oddhorz = sin(2*pi*freq*Y).*gauss2d;
oddvert = sin(2*pi*freq*X).*gauss2d;

% imagesc(oddvert); colormap('gray')

% compute f1, f2, and images
y = zeros(nX,nX,nY);
for it = 1:nT;
    dyn  = (dT/tau)* (-y(:,:,1) + squeeze(inputstimulus(:,:,it)));
    y(:,:,1) = y(:,:,1) + dyn;
    for iy = 2:nY;
        y(:,:,iy) = y(:,:,iy) + (dT/tau) * (-y(:,:,iy) + y(:,:,iy-1)); 
    end
    f1 = y(:,:,3) - y(:,:,5); 
    f2 = y(:,:,5) - y(:,:,7); 
    
    % convolving gabor with impulse response
    evenhorzfast(:,:,it) = conv2(f1,evenhorz);
    evenvertfast(:,:,it) = conv2(f1,evenvert);
    oddhorzfast(:,:,it) = conv2(f1,oddhorz);
    oddvertfast(:,:,it) = conv2(f1,oddvert);
    
    evenhorzslow(:,:,it) = conv2(f2,evenhorz);
    evenvertslow(:,:,it) = conv2(f2,evenvert);
    oddhorzslow(:,:,it) = conv2(f2,oddhorz);
    oddvertslow(:,:,it) = conv2(f2,oddvert);
end

leftresp.even = oddhorzfast + evenhorzslow;
leftresp.odd = -oddhorzslow + evenhorzfast;
rightresp.even = -oddhorzfast + evenhorzslow;
rightresp.odd = oddhorzslow + evenhorzfast;

upresp.even = oddvertfast + evenvertslow;
upresp.odd = -oddvertslow + evenvertfast;
downresp.even = -oddvertfast + evenvertslow;
downresp.odd = oddvertslow + evenvertfast;

leftresp.energy = leftresp.even.^2 + leftresp.odd.^2;
rightresp.energy = rightresp.even.^2 + rightresp.odd.^2;
upresp.energy = upresp.even.^2 + upresp.odd.^2;
downresp.energy = downresp.even.^2 + downresp.odd.^2;
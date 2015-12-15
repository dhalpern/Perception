%% 
% PERCEPTION COLOR ASSIGNMENT
%
% 1)
% Making a graph of the 18th Macbeth surface
%
spectrum = linspace(400,700,31);
load surfaces  
load illuminants
spect18fl = macbeth(18,:)' .* flourescent';
spect18a = macbeth(18,:)' .* cie_a';
plot(spectrum,[spect18fl, spect18a]);
xlabel('Wavelength (nm)');
ylabel('Reflected Energy');
title('Reflected energy of surfaces under flourescent light and illuminant A')
legend('cie a','flourescent');
 
%%
% 
% Cone responses under flourescent light
% 

load cones
coneSignals18fl = cones * spect18fl
%%
% 
% Cone responses under illuminant A
% 

coneSignals18a = cones * spect18a

%%
% 
% Flourescent looks more blueish because S cells are much more active with
% little change in L and M


%% 
% 2)
% Computing the monitor intensities for color matching under illuminant A
% 

load phosphors
monitor_to_cones = cones * phosphors';
cones_to_monitor = inv(monitor_to_cones);
monitorSignals = cones_to_monitor * coneSignals18a

%%
% 
% There are intensities greater than 255 so cannot be displayed on this
% monitor
% 


%% 
% 3)
% Cone responses to baseline stimulus
% 
baseline = phosphors' * [1 1 1]';
baseline_cones =  cones * baseline

%%
% 
% Monitor signal required to increase S cone response by .5
% 

coneSignalsDeltaS = baseline_cones + [0 0 0.5]'; %increment S
monitorSignalsDeltaS = cones_to_monitor * coneSignalsDeltaS % 

%%
% 
% If we increment by .8, we get a negative monitor intensity which is
% impossible
% 

coneSignalsDeltaS = baseline_cones + [0 0 0.8]'; %increment S
monitorSignalsDeltaS = cones_to_monitor * coneSignalsDeltaS
%%
% 
% .7 seems to be highest possible increment (rounded to the nearest tenth)
% 

coneSignalsDeltaS = baseline_cones + [0 0 0.7]'; %increment S
monitorSignalsDeltaS = cones_to_monitor * coneSignalsDeltaS

%% 
% 4)
% Computing color matching function
%
% Given a test light %t%, we can figure out the intensities $e_{cie}$ of the cie
% lights needed to match as $Ct = CM_{cie}e_{cie}$
%
% Rewritten, this becomes $(CM_{cie})^{-1}Ct = e_{cie}$ so the color
% matching function is $(CM_{cie})^{-1}C$
% 
cie = zeros(3, 31);
cie(sub2ind(size(cie), [1 2], [4 16])) = 1;
cie(3, 31) = 80;
color_matching_function = inv(cones * cie') * cones
plot(spectrum,color_matching_function(1,:),'b');
hold on
plot(spectrum,color_matching_function(2,:),'g');
plot(spectrum,color_matching_function(3,:),'r');
xlabel('Wavelength (nm)');
ylabel('Relative Intensity');
title('Color Matching Function')
legend('B','G','R');

%% 
% 5)
% Given a vector of cie intensities, we can compute the requisite phosphor
% intensities withouth reference to a test SPD.
% This is because we can write the SPD of a set of cie intensities as
% $t = M_{cie}e_{cie}$
% So the above equation for the phosphors becomes $(CM_{phos})^{-1}CM_{cie}e_{cie} = e_{phos}$
% which gives us $(CM_{phos})^{-1}CM_{cie}$ as the function to go from cie
% intensities to phosphor intensities
cie_to_cones = cones * cie';
cones_to_phosphors = inv(cones * phosphors');
cie_to_phosphors = cones_to_phosphors * cie_to_cones
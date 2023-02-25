% Run this function to generate plots for AE502_AssignmentOne
set(0,'defaultAxesFontSize',15);

%% Build Mex function
command = '/Applications/MATLAB_R2021b.app/bin/mex AssignmentOneMex.cpp UniversalPropagator.cpp';
[status,cmdout] = system(command);
disp(cmdout)

%% Question 3
A = AssignmentOneMex(int32(1),true); dims = size(A); xDep = 1:dims(1); yArr = 1:dims(2);
DV_1I_Rendezvous_Clipped = A; 
for ii = 1:numel(A)
    if DV_1I_Rendezvous_Clipped(ii) > 50
        DV_1I_Rendezvous_Clipped(ii) = NaN;
    end 
end 

hFig1 = figure; surf(xDep,yArr,DV_1I_Rendezvous_Clipped','EdgeColor','interp','FaceColor','interp');
title('Rendezvous with 1I/’Oumouamoua')
xlabel('Departure Date, days from Jan 2017')
ylabel('Arrival Date, days from Aug 2017')
c = colorbar;
set(get(c,'label'),'string','\DeltaV (km/s)');
view(2)

A = AssignmentOneMex(int32(1),false); dims = size(A); xDep = 1:dims(1); yArr = 1:dims(2);
DV_1I_FlybyClipped = A;
for ii = 1:numel(A)
    if DV_1I_FlybyClipped(ii) > 20
        DV_1I_FlybyClipped(ii) = NaN;
    end 
end 
hFig2 = figure; surf(xDep,yArr,DV_1I_FlybyClipped','EdgeColor','interp','FaceColor','interp');
title('Fly by of 1I/’Oumouamoua')
xlabel('Departure Date, days from Jan 2017')
ylabel('Arrival Date, days from Aug 2017')
c = colorbar;
set(get(c,'label'),'string','\DeltaV (km/s)');
view(2)

%% Question 4 2I/Borisov
A = AssignmentOneMex(int32(2),true);dims = size(A); xDep = 1:dims(1); yArr = 1:dims(2);
DV_2I_Rendezvous_Clipped = A; 
for ii = 1:numel(A)
    if DV_2I_Rendezvous_Clipped(ii) > 60
        DV_2I_Rendezvous_Clipped(ii) = NaN;
    end 
end 

hFig3 = figure; 
surf(xDep,yArr,DV_2I_Rendezvous_Clipped','EdgeColor','interp','FaceColor','interp');
title('Rendezvous with 2I/Borisov')
xlabel('Departure Date, days from Jan 2017')
ylabel('Arrival Date, days from June 2019')
c = colorbar;
set(get(c,'label'),'string','\DeltaV (km/s)');
view(2)

A = AssignmentOneMex(int32(2),false);
DV_2I_FlybyClipped = A;
for ii = 1:numel(A)
    if DV_2I_FlybyClipped(ii) > 20
        DV_2I_FlybyClipped(ii) = NaN;
    end 
end 
hFig4 = figure; 
surf(xDep,yArr,DV_2I_FlybyClipped','EdgeColor','interp','FaceColor','interp');
title('Fly by of 2I/Borisov')
xlabel('Departure Date, days from Jan 2017')
ylabel('Arrival Date, days from June 2019')
c = colorbar;
set(get(c,'label'),'string','\DeltaV (km/s)');
view(2)

% Get directory two levels up
saveDir = '/Users/himmatpanag/Documents/UIUC/AE502 Advanced Orbital Mechanics';
originalDir = pwd; 
hFigs = {hFig1,hFig2,hFig3,hFig4};
cd(saveDir);
for ii = 1:numel(hFigs)
    saveas(hFigs{ii},sprintf('Fig%d.png',ii))
end
cd(originalDir);

%% Question 5
r1I = [3.515868886595499e-2, -3.162046390773074, 4.493983111703389];
v1I = [-2.317577766980901e-3, 9.843360903693031e-3, -1.541856855538041e-2];
r2I = [7.249472033259724, 14.61063037906177, 14.24274452216359];
v2I = [-8.241709369476881e-3, -1.156219024581502e-2, -1.317135977481448e-2];

AU = 149597870.7; % km
day = 86400; % s
mu = 1.32712440018e11; % km^3/s^2
elements1I = GetOrbitElements(r1I*AU,v1I*AU/day,mu);
elements2I = GetOrbitElements(r2I*AU,v2I*AU/day,mu);

%% Question 6 
vEsc = sqrt(2*mu/AU)
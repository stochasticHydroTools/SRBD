clear
clc
close all

% Reset defaults for plotting
set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)
set(0, 'DefaultFigureRenderer', 'opengl');


%% Step 0: Initialize Parameters
blobs = 1500;   % Number of blobs (actin) in the simulation
dimers = 50;    % Number of dimers in the simulation
numParticles = blobs + 2*dimers;
skip = 1;

% Cross linker radii:
a1 = 2*0.01;   %CL Radii
a2 = 1.5*0.008;
%Doi box size (I assume that this is the same across all dimensions)
boxSize = 1;    

%% Step 1: Generate Data
% Read in the file output from fortran code
A = importdata('DimerSimulationIsotropic.txt');
t = linspace(0,length(A)/numParticles,10001);       

% Read in blobs from the file. We only read in a specific amount though.
% as it is too costly to read in all of them (also the simulation in
% fortran wasn't done with all loaded)
actin = importdata('FibersDiameterApart2.txt');
actin = actin(1:blobs,:);

% Create vector of indixies for each new time step.
% Note that this occurs when the first
% column has value 1. Each time step is organized so that it lists the
% blobs in increasing order, so each time step would list (for 100 particles)
% 1-100, then 1-100 again and so on (look through
% DimerSimulationIsotropic.txt to see
timeStep = find(A(:,1) == 1);



%% Step 2: Draw/Render Scenario
% Render spheres
[sx, sy, sz] = sphere(20);
fig1 = figure(1)
k = 0;   % k is for my own debugging purposes, but is not important for the program at present
Lp=40;
[X, Y] = meshgrid([-Lp:0.5:Lp],[-Lp:0.5:Lp]);

counter = 0;
set(fig1,'units','normalized','outerposition',[0 0 1 1])

%iterate across all time steps measured
for i = 1:skip:length(timeStep)   %Not correct in long term simulations
    %Wipe the slate clean so we are plotting with a blank figure
    clf
    counter = counter + 1;
    
    % p - particle, s - species, x,y,z - euclidean coordinates of particle
    p = A(timeStep(i):timeStep(i+1)-1,1); 
    s = A(timeStep(i):timeStep(i+1)-1,2);
    x = A(timeStep(i):timeStep(i+1)-1,3);
    y = A(timeStep(i):timeStep(i+1)-1,4);    
    z = A(timeStep(i):timeStep(i+1)-1,5);
     
    % Set camera light to make the light reflect
    camlight('right');
    
    % Plot the unbound actin seperately from the particles themselves. 
    for k = 1:length(actin)
        % Below is just what I thought looked nice, including a
        % transparency gradient. This could be modified. 
        actinTransparency = 0.5*((norm(actin(k,:)) / sqrt(3)));
        h = surface(actin(k,1) + a2*sx, actin(k,2) + a2*sy, actin(k,3) + a2*sz, 'facecolor',[1.0000, 0.874, 0], 'edgecolor','none','FaceAlpha',actinTransparency,'EdgeAlpha',0);
            set(h,'FaceLighting','gouraud',...%'facealpha',0.2, ...
            'AmbientStrength',0.5, ...
             'DiffuseStrength',0.3, ... 
             'Clipping','off',...
             'BackFaceLighting','lit', ...
             'SpecularStrength',0.3, ...
             'SpecularColorReflectance',0.7, ...
             'SpecularExponent',1)
             hold on
       
       % Commented below were my attempts to draw periodic extensions to
       % the actin. I dropped this for the movie as I decided it wasn't
       % that important, but it is something to improve on.
             
       % Periodic extensions on SOME actin.
       %if (actin(k,1) < 0.1) actinCopyX = actin(k,1) + boxSize; end
       %if (actin(k,2) < 0.1) actinCopyY = actin(k,2) + boxSize; end
       %if (actin(k,3) < 0.1) actinCopyZ = actin(k,3) + boxSize; end
       %if (actin(k,1) > 0.9) actinCopyX = actin(k,1) - boxSize; end
       %if (actin(k,2) > 0.9) actinCopyY = actin(k,2) - boxSize; end
       %if (actin(k,3) > 0.9) actinCopyZ = actin(k,3) - boxSize; end 
       
       %if (actinCopyX ~= actin(k,1) || actinCopyY ~= actin(k,2) || actinCopyZ ~= actin(k,3)) 
       %     h = surface(actinCopyX + a2*sx, actinCopyY + a2*sy, actinCopyZ + a2*sz, 'facecolor',[1.0000, 0.874, 0], 'edgecolor','none','FaceAlpha',actinTransparency,'EdgeAlpha',0);
       %     set(h,'FaceLighting','gouraud',...%'facealpha',0.2, ...
       %     'AmbientStrength',0.5, ...
       %     'DiffuseStrength',0.3, ... 
       %      'Clipping','off',...
       %      'BackFaceLighting','lit', ...
       %      'SpecularStrength',0.3, ...
       %      'SpecularColorReflectance',0.7, ...
       %      'SpecularExponent',1)
       %end
    end
    
        
    % Iterate over all (non-blob) particles
    for j = 1:length(x)
        %Plot the current location of the particle
        % Booleans for whether particle has exited the box. 
        % Calculate Distance Components (FOR PBC - periodic boundary condition)
        if (mod(j,2) == 1 && (s(j) == 1 || s(j) == 3))
            % s(j) = 1, or 3 are cross linkers that are unbound or bound
            dx = x(j) - x(j+1);
            dy = y(j) - y(j+1);
            dz = z(j) - z(j+1);
            
            % Apply periodic wrapping
            
            % particle to the left
            if (dx > boxSize/2)
                dx = dx - boxSize; 
            end
            if (dy > boxSize/2)
                dy = dy - boxSize; 
            end
            if (dz > boxSize/2)
                dz = dz - boxSize; 
            end

            % particle to the right
            if (dx < -boxSize/2)
                dx = dx + boxSize; 
            end
            if (dy < -boxSize/2)
                dy = dy + boxSize; 
            end
            if (dz < -boxSize/2)
                dz = dz + boxSize; 
            end
            
            % If bound CL, then get rid of periodic wrapping (this is so
            % you do not have two ends of a bound CL on opposite sides of
            % the box)
            if (s(j) == 3)
                x(j+1) = x(j) - dx;
                y(j+1) = y(j) - dy;
                z(j+1) = z(j) - dz;
            else 
                x(j) = dx + x(j+1); 
                y(j) = dy + y(j+1); 
                z(j) = dz + z(j+1); 
            end
            
        end
        
        k = k + 1;
        % Set color schemes
        pink = [1.000000000000000   0.074509803921569   0.650980392156863];
        cyan = [0.039215686274510   0.796078431372549   0.811764705882353];
        
        % Iterate over all cases for the species (not case == 2 should
        % never happen as we handle blobs seperately, but I have it
        % included from a previous version)
        switch s(j)
            % Unbound CL
            case 1
                % Distinguish between cases where one end of the CL is bound:
                % 1. ) If one end of the CL is bound, then both ends of the CL
                % (bound and unbound) must be more opaque
                % 2. ) If both ends are unbound, then must be more
                % translucent
                if ((j < (length(x)) && ((s(j+1) == 3 && mod(j,2) == 1) || (mod(j,2) == 0 && s(j-1) == 3))))
                    h = surface(x(j) + a1*sx, y(j) + a1*sy, z(j) + a1*sz, 'facecolor',cyan, 'edgecolor','none','FaceAlpha',1.0,'EdgeAlpha',1.0);
                    set(h,'FaceLighting','gouraud',...%'facealpha',0.2, ...
                    'AmbientStrength',0.5, ...
                    'DiffuseStrength',0.3, ... 
                    'Clipping','off',...
                    'BackFaceLighting','lit', ...
                    'SpecularStrength',0.3, ...
                    'SpecularColorReflectance',0.7, ...
                    'SpecularExponent',1)
                else
                    h = surface(x(j) + a1*sx, y(j) + a1*sy, z(j) + a1*sz, 'facecolor',cyan, 'edgecolor','none','FaceAlpha',0.2,'EdgeAlpha',0.2);
                    set(h,'FaceLighting','gouraud',...%'facealpha',0.2, ...
                    'AmbientStrength',0.5, ...
                    'DiffuseStrength',0.3, ... 
                    'Clipping','off',...
                    'BackFaceLighting','lit', ...
                    'SpecularStrength',0.3, ...
                    'SpecularColorReflectance',0.7, ...
                    'SpecularExponent',1)
                end
            % Unbound actin
            case 2
                h = surface(x(j) + a2*sx, y(j) + a2*sy, z(j) + a2*sz, 'facecolor',[0.878431372549020   0.737254901960784   0.317647058823529], 'edgecolor','none','FaceAlpha',0.2,'EdgeAlpha',0.2);
                set(h,'FaceLighting','gouraud',...%'facealpha',0.2, ...
                'AmbientStrength',0.5, ...
                'DiffuseStrength',0.3, ... 
                'Clipping','off',...
                'BackFaceLighting','lit', ...
                'SpecularStrength',0.3, ...
                'SpecularColorReflectance',0.7, ...
                'SpecularExponent',1)
            % Bound CL
            case 3
                h = surface(x(j) + a1*sx, y(j) + a1*sy, z(j) + a1*sz, 'facecolor',pink, 'edgecolor','none','FaceAlpha',1.0,'EdgeAlpha',1.0);
                set(h,'FaceLighting','gouraud',...%'facealpha',0.2, ...
                'AmbientStrength',0.5, ...
                'DiffuseStrength',0.3, ... 
                'Clipping','off',...
                'BackFaceLighting','lit', ...
                'SpecularStrength',0.3, ...
                'SpecularColorReflectance',0.7, ...
                'SpecularExponent',1)
            % Bound actin
            case 4
                h = surface(x(j) + a2*sx, y(j) + a2*sy, z(j) + a2*sz, 'facecolor',[0.149019607843137   0.149019607843137   0.149019607843137], 'edgecolor','none','FaceAlpha',1.0,'EdgeAlpha',1.0);
                set(h,'FaceLighting','gouraud',...%'facealpha',0.2, ...
                'AmbientStrength',0.5, ...
                'DiffuseStrength',0.3, ... 
                'Clipping','off',...
                'BackFaceLighting','lit', ...
                'SpecularStrength',0.3, ...
                'SpecularColorReflectance',0.7, ...
                'SpecularExponent',1) 
        end
        
        % All particles at one time step should be plotted
        hold on 

        %Plot lines between cross-linkers, but light gray if not bound, and
        %black if bound
        if (mod(j,2) == 1 && s(j) == 1 && s(j+1) == 1)
            pts = [x(j), y(j), z(j);x(j+1),y(j+1),z(j+1)];
            plot3(pts(:,1),pts(:,2),pts(:,3),'color',[0.8 0.8 0.8],'Linewidth',3)
        elseif(mod(j,2) == 1 && (s(j) == 1 || s(j) == 3))
            pts = [x(j), y(j), z(j);x(j+1),y(j+1),z(j+1)];
            plot3(pts(:,1),pts(:,2),pts(:,3),'k','Linewidth',3)
        end
        
        
        % Decorate the plot
        set(gca,'visible','off');
        daspect([1 1 1])
        view([-20 35]) %view([0 90]) %view([-140 10])% 
        axis([-0.1 1.1 -0.1 1.1 -0.1 1.1])

        %% Step 4: Advance time
        % Accomplished by the for loop
    end
    % get current axis and turn it off
    set(gca,'visible','off')
    drawnow
    
    % Print out figures to the directory below (may need to modify this)
    print('-dpng',['./Isotropic/crossLinker_' num2str(counter) '.png'],'-r200')

end





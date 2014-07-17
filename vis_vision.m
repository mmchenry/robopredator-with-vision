function vis_vision
% Visualizes the visual stimulus presented to larvae in our experiments


%% Set paths

% This should be specific to each computer that executes this code
if isdir('/Users/mmchenry/')
    root = '/Users/mmchenry/Documents/Projects/Robopredator with vision';
else
    root = '/Users/arjunnair0513/Dropbox/Shared with Arjun';
end


%% Code execution


% Save frames for a movie
save_images = 1;

% Clean look for presentations 
vis_clean = 1;


%% Parameter values

% Color of 3D rendering of predator
p_clr = .5.*[1 1 1];

% Color of 3D rendering of prey
py_clr = .5.*[1 1 1];

% Colors for the 3 speeds
sp_clr(1,:) = [0 0 1];
sp_clr(2,:) = [0 1 0];
sp_clr(3,:) = [1 0 0];

% Scale factor for spherical heads of larvae
scl_fctr = 0.0005;

% Receptive field of visual system (m, rad). Easter (1996) give the functional 
% retinal field of 163 deg for 4 dpf. Though I give an additional 20 deg
% for saccades (based on vergenece angles of Patterson et al, 2013)
vis_depth  = 0.15;
vis_az     = (163+20)/180*pi;
vis_el     = (163)/180*pi;
off_el     = 0;
vis_num    = 100;

% Deg per photorecptor pair (Easter, 1996) for 4 dpf
ren_den = 2.9235/180*pi;

% Number of points to define predator's shape
p_num = 200;

% Number of position points for the predator's motion
num_pos = 100;

% Vergence angle: mean = 32 deg, max = 70 deg (during predation) (Patterson et al, 2013)
verg_ang = 32/180*pi;

% Radial coordinates 
theta = linspace(0,2*pi,vis_num)';

% Which approach speed from which to choose a sequence to visualize
spd_num = 2;

% Indivual to visualize
indiv_num = 10;

% Grid color
g_clr = .5.*[1 1 1];


%% Load and organize data

% Load Kelsey's pooled data ('L')
load([root filesep 'kelsey data' filesep 'pooled_transdata.mat']) 

% Load predator data 
load([root filesep 'bills data' filesep 'Pred3Dbodyshape.mat'])
M.x = pred3DshapeX./100;
M.y = pred3DshapeY./100;
M.z = pred3DshapeZ./100;

clear pred3DshapeX pred3DshapeY pred3DshapeZ

% Load prey morphology data, define dimensions ('m')
load([root filesep 'bills data' filesep '6_02_L19_metrics.mat']);

% Load visual system modeling data ('F')
load([root filesep 'kelsey data' filesep 'Visual data.mat'])   

% Extract data for visualizing
pred_pos   = F(spd_num).pred_pos(:,indiv_num);
imR        = reshape(F(spd_num).imR(:,:,indiv_num,:),...
                     size(F(spd_num).imR(:,:,indiv_num,:),1),...
                     size(F(spd_num).imR(:,:,indiv_num,:),2),...
                     size(F(spd_num).imR(:,:,indiv_num,:),4));
imL        = reshape(F(spd_num).imL(:,:,indiv_num,:),...
                     size(F(spd_num).imL(:,:,indiv_num,:),1),...
                     size(F(spd_num).imL(:,:,indiv_num,:),2),...
                     size(F(spd_num).imL(:,:,indiv_num,:),4));
time       = F(spd_num).time(:,indiv_num);
head       = L(spd_num).head(indiv_num,:);
tail       = L(spd_num).tail(indiv_num,:);
com        = L(spd_num).com(indiv_num,:);
Raz_f      = F(spd_num).R_azF./pi.*180;
Rel_f      = F(spd_num).R_elF./pi.*180;
Raz_field  = F(spd_num).R_azField./pi.*180;
Rel_field  = F(spd_num).R_elField./pi.*180;

Laz_f      = F(spd_num).L_azF./pi.*180;
Lel_f      = F(spd_num).L_elF./pi.*180;
Laz_field  = F(spd_num).L_azField./pi.*180;
Lel_field  = F(spd_num).L_elField./pi.*180;

% Coordinates for rendered prey
[Xprey,Yprey,Zprey] = prey_surf(m.s./100,m.h./100,m.w./100,m.c./100,...
                                p_num,head,com,tail);

%% Render simulations                            
                            
% Step through time
for i = 1:length(pred_pos)
    
    % Initialize figure
    if (i==1) 
        hF = figure;
        set(hF,'DoubleBuffer','on')
        light_on = 1;
        
        if vis_clean
            set(hF,'Color','k')
        end
    else
        light_on = 0;
    end
    
    % Render the predator in 3D from side
    subplot(4,4,1:4) %-----------------------------------------------------
    hP1 = render_surf(M.x + pred_pos(i),M.y,M.z,p_clr,light_on);
    hold on
    
    % Render the prey in 3D
    if i==1
        hPrey = render_surf(Xprey,Yprey,Zprey,py_clr,light_on);
    end
    
    % Set title (time)
    hT = title([num2str(round(time(i).*1000)) ' ms']);
    
    % Set axis appearence
    if ~vis_clean
        xlabel('X'); ylabel('Y');zlabel('Z')
    else
        grid off
        set(gca,'Color',0.*[1 1 1]) 
        set(gca,'XTickLabel',' ')
        set(gca,'YTickLabel',' ')
        set(gca,'ZTickLabel',' ')
        set(hT,'Color','w')
    end    
    
    view([0 0])
    if i == 1
        axis tight
        xL = xlim;
        pos1 = get(gca,'Position');
    else
        axis equal
        axis tight
        xlim(xL)
    end
    
    % Render the predator in 3D from above
    subplot(4,4,5:8) %-----------------------------------------------------
    hP2 = render_surf(M.x + pred_pos(i),M.y,M.z,p_clr,light_on);
    if vis_clean
       set(gca,'Color','k') 
    end
    hold on
    
    % Render the prey in 3D
    if i==1
        hPrey = render_surf(Xprey,Yprey,Zprey,py_clr,light_on);
    end
    
    % Set axis appearence
    if ~vis_clean
        xlabel('X'); ylabel('Y');zlabel('Z')
    else
        grid off
        set(gca,'Color',0.*[1 1 1]) 
        set(gca,'XTickLabel',' ')
        set(gca,'YTickLabel',' ')
        set(gca,'ZTickLabel',' ')

    end
    view([0 90])
    
    if i == 1
        axis equal
        axis tight
        xlim(xL)
        pos2 = get(gca,'Position');
        set(gca,'Position',[pos2(1) pos2(2) pos1(3) pos1(4)]);
    else
        axis equal
        axis tight
        xlim(xL)
    end
    
    
    % Plot view from left eye
    subplot(4,4,[9 10 13 14]) %------------------------------------------------------
    
    ax = set_polar_plot(Lel_field,Laz_field,vis_clean);    
    h = geoshow(Lel_field,Laz_field,double(~imL(:,:,i).*255),...
        'DisplayType','texturemap');
    %set(gca,'XDir','reverse')
    
    axis equal
    colormap([.3.*[1 1 1]; [1 1 1]])
    %set(gca,'YDir','normal')
    hT = title('Left');
    if vis_clean
        set(hT,'Color',0*[1 1 1])
    end
    axis tight
    hold off
    
    % Plot view from right eye
    subplot(4,4,[11 12 15 16]) %------------------------------------------------------
    ax = set_polar_plot(Rel_field,Raz_field,vis_clean);    
    h = geoshow(Rel_field,Raz_field,double(~imR(:,:,i).*255),...
        'DisplayType','texturemap');
    %set(gca,'XDir','reverse')
    
    axis equal
    colormap([.3.*[1 1 1]; [1 1 1]])
    %set(gca,'YDir','normal')
    hold off
    hT = title('Right');
    if vis_clean
        set(hT,'Color',0*[1 1 1])
    end
    axis tight
    
    pause(.001)
    
    % Save image (Frame)
    if save_images
        
        % Capture frame image
        im = getframe(gcf);
        
        
        % Define filename
        F_txt = ['00' num2str(i)];
        F_txt = F_txt(end-2:end);
        fname = ['spd ' num2str(L(spd_num).spd) ' larva ' ...
                 num2str(indiv_num) ' frame ' F_txt '.jpg'];
        
        % Write file
        imwrite(im.cdata,[root filesep 'animations' filesep 'trial' fname],'jpg');
    end
    
    % Set up for next iteration
    delete(hP1), delete(hP2)
    
    clear s_txt F_txt light_on im h hT hP1 hP2
end

beep;beep;



function ax = set_polar_plot(el,az,vis_clean)

ax = axesm('ortho',...
    'Grid','on',...
    'GLineStyle','-',...
    'GColor',.8.*[1 1 1],...
    'LabelFormat','signed');

%'MapLatLimit',[min(az(:)) max(az(:))],...
%    'MapLonLimit',[min(el(:)) max(el(:))],...

if ~vis_clean
    setm(ax,'meridianlabel','on')
    setm(ax,'parallellabel','on')
    setm(ax,'plabellocation',90)
    setm(ax,'mlabellocation',0)
    setm(ax,'mlinelocation',30)
    setm(ax,'plinelocation',30)
    setm(ax,'mlabelparallel',-90)
    setm(ax,'plabelmeridian',0)
end
axis off
framem; 
%gridm;


function plot_larvae(head,tail,clr,msize)

if size(head,2) == 2
    for i = 1:size(head,1)
        
        h = plot(head(i,1),head(i,2),'o',[head(i,1) tail(i,1)],...
                 [head(i,2) tail(i,2)],'-');
        set(h(1),'MarkerFaceColor',clr,'MarkerEdgeColor',clr,...
                 'MarkerSize',msize)
        set(h(2),'Color',clr)      
        
    end
end


function h = render_surf(x,y,z,p_clr,light_on)

% Plot surfaces
h = patch(x,y,z,z*0);

% Adjust parameters
set(h,'FaceLighting','gouraud',...
    'LineStyle','none',...
    'BackFaceLighting','reverselit',...
    'FaceColor',p_clr,...
    'AmbientStrength',.5);

% Lighting, & view control --------------------------------
if light_on
    hL(1) = light('position',[0 0 20]);
    hL(2) = light('position',[0 0 -20]);
end

axis equal

xlabel('X'); ylabel('Y');zlabel('Z')
view([-180 0]);
grid on



function [X,Y,Z]= prey_surf(s,h,w,c,numPts,rost,com,tail)
% Provides 3d coordinates of the body

% Define radial positions along vector
theta = linspace(0,2*pi,numPts)';

% Define empty vectors for coordinates
x=[];y=[];z=[];

% Make mouth cap  _______________________________________
n = numPts/10;
phi = linspace(0,.75*pi/2,n)';
ds = .02.*range(s); %2*s(2)-s(1);
%sC = linspace(s(1)-ds,s(1),n);
hC = h(1) .* sin(phi)./max(sin(phi));
wC = w(1) .* sin(phi)./max(sin(phi));
sC = -(ds.*cos(phi)-ds.*cos(phi(end)));

% Loop down the body length
for i=1:length(sC)-1  
    
  % Draw first ellipse   
    yTemp1 = sC(i)*ones(size(theta));
    xTemp1 = (wC(i)/2) .* cos(theta);
    zTemp1 = (hC(i)/2) .* sin(theta) + c(1);
    
  % Draw second ellipse  
    yTemp2 = sC(i+1)*ones(size(theta));
    xTemp2 = (wC(i+1)/2) .* cos(theta);
    zTemp2 = (hC(i+1)/2) .* sin(theta) + c(1);
    
  % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
end 

clear xTemp1 yTemp1 zTemp1 xTemp2 yTemp2 zTemp2
clear n phi ds sC hC wC


% Make body coordinates _______________________________________

% Loop down the body length
for i=1:length(s)-1  
    
  % Draw first ellipse  
    yTemp1      = s(i)*ones(size(theta));
    xTemp1      = (w(i)/2) .* cos(theta);
    zTemp1      = (h(i)/2) .* sin(theta) + c(i);
    
  % Draw second ellipse    
    yTemp2      = s(i+1)*ones(size(theta));
    xTemp2      = (w(i+1)/2) .* cos(theta);
    zTemp2      = (h(i+1)/2) .* sin(theta) + c(i+1);
    
  % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
end  

% Make tail cap  _______________________________________
n = numPts/10;
phi = linspace(0,0.75*pi/2,n)';
ds = .02.*range(s); %2*s(2)-s(1);
%sC = linspace(s(1)-ds,s(1),n);
hC = h(end) .* sin(phi)./max(sin(phi));
wC = w(end) .* sin(phi)./max(sin(phi));
sC = s(end) + ds.*cos(phi)-+ ds.*cos(phi(end));

% Loop down the body length
for i=1:length(sC)-1  
    
  % Draw first ellipse   
    yTemp1 = sC(i)*ones(size(theta));
    xTemp1 = (wC(i)/2) .* cos(theta);
    zTemp1 = (hC(i)/2) .* sin(theta) + c(end);
    
  % Draw second ellipse  
    yTemp2 = sC(i+1)*ones(size(theta));
    xTemp2 = (wC(i+1)/2) .* cos(theta);
    zTemp2 = (hC(i+1)/2) .* sin(theta) + c(end);
    
  % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
end 

clear xTemp1 yTemp1 zTemp1 xTemp2 yTemp2 zTemp2
clear n phi ds sC hC wC

% Transform for units and wrt rostrum
X  = (y-min(y(:))).*100;
Y  = x.*100;
Z  = -(z-z(1)+max(h(:))/2).*100;

% Transform in the global FOR
[X,Y,Z] = prey_global(rost,com,tail,X,Y,Z);


function [pXg,pYg,pZg] = prey_global(rost,com,tail,pX,pY,pZ)
% Transforms prey points in global FOR

tmp1 = local_to_global(rost,com,tail,...
                        [pX(1,:)' pY(1,:)' pZ(1,:)']);
tmp2 = local_to_global(rost,com,tail,...
                        [pX(2,:)' pY(2,:)' pZ(2,:)']);
tmp3 = local_to_global(rost,com,tail,...
                        [pX(3,:)' pY(3,:)' pZ(3,:)']);
tmp4 = local_to_global(rost,com,tail,...
                        [pX(4,:)' pY(4,:)' pZ(4,:)']);

pXg = [tmp1(:,1)'; tmp2(:,1)'; tmp3(:,1)'; tmp4(:,1)'];
pYg = [tmp1(:,2)'; tmp2(:,2)'; tmp3(:,2)'; tmp4(:,2)'];
pZg = [tmp1(:,3)'; tmp2(:,3)'; tmp3(:,3)'; tmp4(:,3)'];


function ptsT = local_to_global(rost,com,tail,pts)
% Transforms coordinates (pts, in n x 3) from global to local coordinate
% system, assuming the y-axis of larvae runs perpendicular to gravity

% Check dimensions of landmark coordinates
if size(rost,1)~=1 || size(rost,2)~=3 || size(com,1)~=1 | size(com,2)~=3 ...
   size(tail,1)~=1 || size(tail,2)~=3
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = tail(1) - rost(1);
xaxis(1,2) = tail(2) - rost(2);
xaxis(1,3) = tail(3) - rost(3);

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis' yaxis' zaxis'];

% If points are not meshgridded
if size(pts,2)==3
    % Rotate points
    ptsT = [R * pts']';
    
    % Translate global coordinates wrt rostrum
    ptsT(:,1) = ptsT(:,1) + rost(1);
    ptsT(:,2) = ptsT(:,2) + rost(2);
    ptsT(:,3) = ptsT(:,3) + rost(3);
    
else
    error('points need to be arranged in 3 columns')
    
end


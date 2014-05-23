function model_vision
% Models the visual stimulus for each experiment.


%% Set paths

% This should be specific to each computer that executes this code
root = '/Users/mmchenry/Dropbox/Projects/Robopredator with vision';

%% Code execution

% Model visual stimuli (takes a while)
do_model_visual = 1;

% Visualize results of visual cue tests
vis_cue_tests = 1;

% Save frames for a movie
save_images = 0;

% Analyze visual cues from visual calculation
do_ana_cue = 1;


%% Parameter values

% Color of 3D rendering of predator
p_clr = .5.*[1 1 1];

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

% Initial position of the predator (larva responded at x = 0)
x0 = -12e-2;

% Number of points to define predator's shape
p_num = 500;

% Number of position points for the predator's motion
num_pos = 100;

% Vergence angle: mean = 32 deg, max = 70 deg (during predation) (Patterson et al, 2013)
verg_ang = 32/180*pi;

% Radial coordinates 
theta = linspace(0,2*pi,vis_num)';

% Latency btwn stim & response (s)
latency = 200e-3;


%% Organize data

% Load Kelsey's pooled data ('L')
load([root filesep 'kelsey data' filesep 'pooled_transdata.mat']) 

% Load predator data ('M')
load([root filesep 'bills data' filesep 'Pred3Dbodyshape.mat'])
M.x = pred3DshapeX./100;
M.y = pred3DshapeY./100;
M.z = pred3DshapeZ./100;

clear pred3DshapeX pred3DshapeY pred3DshapeZ


%% Model the sight of the predator in the prey FOR (do_model_visual)

if do_model_visual %isempty(dir([root filesep 'Visual data.mat']))
    
tic

clear tmpR tmpL ptsL ptsR a b

% Calculate time vector
%time = linspace(-7e-2,0,num_pos)';

% Coordinates for the ellipse of the FOV
%xFOV = (vis_az/2).*cos(theta);
%yFOV = (vis_el/2).*sin(theta);

% Max potential retinal field
el_f = (vis_el/2):-ren_den:-(vis_el/2);
az_f = (vis_az/2):-ren_den:-(vis_az/2);
[azField,elField] = meshgrid(az_f,el_f);

% Loop thru speeds
for i = 1:3
    
    disp(' ');disp(['Starting spd = ' num2str(L(i).spd) ' cm/s']) 
    
    % Position of predator
    %pred_pos = P(i).spd .* time;
    
    % Store stats
    F(i).spd       = L(i).spd;
    F(i).L_elF     = el_f;
    F(i).L_azF     = az_f;
    F(i).L_azField = azField;
    F(i).L_elField = elField;
    F(i).R_elF     = el_f;
    F(i).R_azF     = -az_f;
    F(i).R_azField = -azField;
    F(i).R_elField = elField;
    
    %clear el_f az_f azField elField
    
    
    % Loop thru individual larvae
    for j = 1:size(L(i).head,1)
        
        % Define head and tail coordinates for current larva
        head = L(i).head(j,:);
        tail = L(i).tail(j,:);
        
        % Predator position to end sequence
        x_end = max([0 head(1)]);
        
        % Response time
        t_resp = -head(1)/L(i).spd;
        
        % Predator displacement
        dx = (x_end - x0)/num_pos;
        
        % Predator position
        pred_pos = (x0+dx):dx:x_end;
        
        % Change in time
        dt = dx/(L(i).spd./100);
        
        % Time (s)
        time = [(-(x_end - x0)/(L(i).spd./100)+dt):dt:0] + max([0 t_resp]);
        
        % Store data for current individual
        F(i).time(:,j)      = time;
        F(i).t_resp(1,j)    = t_resp;
        F(i).pred_pos(:,j)  = pred_pos;
        F(i).head(j,:)      = head;
        F(i).tail(j,:)      = tail;
        
        clear x_end t_resp dx pred_pos dt time
        
        disp(' ')
        disp([' Starting ' num2str(j) ' of ' ...
                  num2str(size(L(i).head,1)) ' larvae ('...
                  num2str(L(i).spd.*100) ' cm/s)']); disp(' ')
        
        % Loop thru predator position values
        for k = 1:length(F(i).pred_pos(:,j))
            
            % Test coordinates
            %head = [20e-3 0 10e-3]
            %tail = [(20e-3+3e-3*sin(pi/6)) 3e-3*cos(pi/6) 0]
            %tail = [(20e-3+3e-3) 0 10e-3]
            %vis_az = pi
            %vis_el = pi/20
            %verg_ang = pi/4+0*10/180*pi;

            % Define coordinate systems for the two eyes
            [rS,lS] = eye_systems(head,tail,verg_ang,vis_az);
            
            % Rotation matrix 
            R = prey_system(head,tail);
            
            % Define predator downsampled inertial coords for current frame
            pdX = M.x + F(i).pred_pos(k,j);
            pdY = M.y;
            pdZ = M.z;
            
            % Transform predator in prey FOR (Right Eye)
            [pdXR,pdYR,pdZR] = prey_local(head,rS,pdX,pdY,pdZ);
            
            % Transform predator in prey FOR (Left Eye)
            [pdXL,pdYL,pdZL] = prey_local(head,lS,pdX,pdY,pdZ);   
            
            % Accept points along positive y-axis of eye FOR (right eye)
            idx = pdYR>=0;
            pdXR = pdXR(idx);
            pdYR = pdYR(idx);
            pdZR = pdZR(idx);
            
            % Accept points along positive y-axis of eye FOR (left eye)
            idx = pdYL>=0;
            pdXL = pdXL(idx);
            pdYL = pdYL(idx);
            pdZL = pdZL(idx);

            % Convert coordinates into angular dimensions
            [R_th,R_phi,R_r] = cart2sph(pdYR,pdXR,pdZR);
            [L_th,L_phi,L_r] = cart2sph(pdYL,pdXL,pdZL);
                  
            clear idx R_r L_r

            % Find delaunay triangulations
            warning off
            R_dt = delaunayTriangulation(R_th,R_phi);
            L_dt = delaunayTriangulation(L_th,L_phi);
            warning on
            
            % Find coordinates, if any present
            if ~isempty(R_th)
                % Identify points in the triangulation
                R_id = pointLocation(R_dt,F(i).R_azField(:),F(i).R_elField(:));
            else
                R_id = nan(size(F(i).R_azField(:)));
            end
            
            % Find coordinates, if any present
            if ~isempty(L_th)
                % Identify points in the triangulation
                L_id = pointLocation(L_dt,F(i).L_azField(:),F(i).L_elField(:));

            else
                L_id = nan(size(F(i).L_azField(:)));
            end
            
            % Create blank image data to show identified images
            imR = false(size(F(i).R_azField)); 
            imL = false(size(F(i).L_azField)); 
            
            % Fill in predator image
            imR(~isnan(R_id)) = true;
            imL(~isnan(L_id)) = true;

            % Save images       
            F(i).imR(:,:,j,k) = imR;
            F(i).imL(:,:,j,k) = imL;
          
            % Visually confirm image creation
            if 0
                figure
                
                subplot(2,2,1)
                h = plot(L_th,L_phi,'k.');
                set(h,'MarkerEdgeColor',.5.*[1 1 1])
                axis equal
                hold on
                h = plot(F(i).azField,F(i).elField,'ok');
                set(h,'MarkerSize',5)
                h = plot(F(i).azField(~isnan(L_id)),F(i).elField(~isnan(L_id)),'ro');
                set(h,'MarkerFaceColor','r')
                set(h,'MarkerSize',5)
                title('Left eye')
                
                subplot(2,2,2)
                h = plot(R_th,R_phi,'k.');
                set(h,'MarkerEdgeColor',.5.*[1 1 1])
                axis equal
                hold on
                h = plot(F(i).azField,F(i).elField,'ok');
                set(h,'MarkerSize',5)
                h = plot(F(i).azField(~isnan(R_id)),F(i).elField(~isnan(R_id)),'ro');
                set(h,'MarkerFaceColor','r')
                set(h,'MarkerSize',5)
                title('Right eye')    
                
                subplot(2,2,3)
                h = image(F(i).az_f,F(i).el_f,imL);
                axis equal
                colormap([0 0 0; 1 1 1])
                set(gca,'YDir','normal')
                   
                subplot(2,2,4)
                h = image(F(i).az_f,F(i).el_f,imR);
                axis equal
                colormap([0 0 0; 1 1 1])
                set(gca,'YDir','normal')
            end
            
            
            disp(['     Done ' num2str(k) ' of ' ...
                  num2str(length(F(i).pred_pos(:,j))) ' positions'])
              
            clear imR imL 
        end
    end   
end

clear az_f el_f azField elField i j k 


% Save data
save([root filesep 'kelsey data' filesep 'Visual data.mat'],'F')

% Report processing time
tlapse = toc;
disp(' ');
disp(['Total processing time = ' num2str(tlapse/60) ' min'])
beep;pause(.3);beep;pause(.3);beep;


end


%% Organize stimulus data

if do_ana_cue
            
% Load data in 'F' structure
load([root filesep 'kelsey data' filesep 'Visual data.mat'])    
    
% Loop thru approach speeds
for i = 1:length(F)
   
    % Loop thru individuals
    for j = 1:size(F(i).head,1)
        
        % Store sequence data
        D.spd(j,1)      = F(i).spd;
        D.time(j,:)     = F(i).time(:,j)';
        D.pred_pos(j,:) = F(i).pred_pos(:,j)';
        D.t_resp(j,1)   = F(i).t_resp(j);
        D.t_stim(j,1)   = F(i).t_resp(j) - latency;
        D.behav(j,1)    = L(i).behav(j);
        
        % Loop thru images for each instant of time
        for k = 1:size(F(i).imR,4)
            
            % Stimulated area (deg^2)
            D.Rarea(j,k) = sum(sum(F(i).imR(:,:,j,k))).*(ren_den/pi*180)^2;
            D.Larea(j,k) = sum(sum(F(i).imL(:,:,j,k))).*(ren_den/pi*180)^2;
            
            % Angular position of the area's centroid
            
            
            % Angular range (deg) of stimulus in az and el
            D.RazRange(j,k) = max([0 range(F(i).R_azField(F(i).imR(:,:,j,k)))/pi*180]);
            D.RelRange(j,k) = max([0 range(F(i).R_elField(F(i).imR(:,:,j,k)))/pi*180]);
            D.LazRange(j,k) = max([0 range(F(i).L_azField(F(i).imL(:,:,j,k)))/pi*180]);
            D.LelRange(j,k) = max([0 range(F(i).L_elField(F(i).imL(:,:,j,k)))/pi*180]);
        end  
        
        % Logical for whether right eye got the bigger stim
        D.R_facing(j,1) = max(D.Rarea(j,:)) > max(D.Larea(j,:));
        
        % Logical for binocular stimulus
        D.binoc(j,1)    = sum(D.Rarea(j,:))~=0 & sum(D.Larea(j,:))~=0;
        
        % Find first derivative of each parameter
        D.time_rate(j,:)  = D.time(j,1:(end-1)) + mean(diff(D.time(j,:)))/2;
        D.Rarea_rate(j,:) = diff(D.Rarea(j,:)) ./ diff(D.time(j,:));
        D.Larea_rate(j,:) = diff(D.Larea(j,:)) ./ diff(D.time(j,:));
        D.RazRange_rate(j,:) = diff(D.RazRange(j,:))./ diff(D.time(j,:));
        D.RelRange_rate(j,:) = diff(D.RelRange(j,:))./ diff(D.time(j,:));
        D.LazRange_rate(j,:) = diff(D.LazRange(j,:))./ diff(D.time(j,:));
        D.LelRange_rate(j,:) = diff(D.LelRange(j,:))./ diff(D.time(j,:));
    end
end

% Save data
save([root filesep 'kelsey data' filesep 'Visual data (organized).mat'],'D')

% Clear variables
clear F
   
else
    % Display results
    disp(' '); disp('Loading organized visual data . . .'); disp(' '); 
    
    % Load data in 'D' structure
    load([root filesep 'kelsey data' filesep 'Visual data (organized).mat'])
end


%% Visualize stimulus data

if vis_cue_tests
    
    % TODO: display fast start responses separately from swimming responses
    
    % Plot simple variables
    vars{1} = 'area';
    vars{2} = 'azRange';
    vars{3} = 'elRange';
    
    figure
    plot_series(D,vars)
    
    figure
    plot_hists(D,vars)
    
    clear vars
    
    % Plot rate of change in simple variables
    vars{1} = 'area_rate';
    vars{2} = 'azRange_rate';
    vars{3} = 'elRange_rate';
    
    figure
    plot_series(D,vars)
    
    figure
    plot_series(D,vars)
    
    clear vars
    
end


%% Visualize comparison boxplots

if vis_cue_tests

vars{1} = 'area';
vars{2} = 'azRange';
vars{3} = 'elRange';
vars{4} = 't_resp';

vars{5} = 'area_rate';
vars{6} = 'azRange_rate';
vars{7} = 'elRange_rate';


figure
subplot(2,1,1)
plot_box(D,vars)

clear vars

end


function plot_box(D,vars)

% Step thru each variable
for i = 1:length(vars)
    
    % Response time
    if strcmp(vars{i},'t_resp')
        
        thresh_val(:,i) = D.t_resp;
    
    % Variables for both R and L eyes
    else
        eval(['Rvals = D.R' vars{i} ';'])
        eval(['Lvals = D.L' vars{i} ';'])
        
        if strcmp(vars{i}(end-3:end),'rate')
            time = D.time_rate;
        else
            time = D.time;
        end
        
        % Interpolate threshold value for each individual
        for j = 1:size(Rvals,1)
            
            if D.R_facing(j)
                thresh_val(j,i) = interp1(time(j,:),Rvals(j,:),D.t_stim(j));
            else
                thresh_val(j,i) = interp1(time(j,:),Lvals(j,:),D.t_stim(j));
            end
        end
        
    end
    
    thresh_val(:,i) = abs(thresh_val(:,i));
    
    outliers = ( (thresh_val(:,i) - mean(thresh_val(:,i))) ...
                 > 3.*std(thresh_val(:,i)) ) & ... 
               (thresh_val(:,i)==0);
    
    % Normalize by median
    thresh_val(:,i) = thresh_val(:,i)./mean(thresh_val(~outliers,i));
end

% Remove zeros
thresh_val(thresh_val==0) = nan;
boxplot(thresh_val,vars)
ylabel('Normalized Thresh. value')


function plot_hists(D,vars)

for i = 1:length(vars)
    
    eval(['Rvals = D.R' vars{i} ';']) 
    eval(['Lvals = D.L' vars{i} ';']) 
    
    if strcmp(vars{i}(end-3:end),'rate')
        time = D.time_rate;
    else
        time = D.time;
    end
        

    for j = 1:size(Rvals,1)
        
        if D.R_facing(j)           
            thresh_val(j,1) = interp1(time(j,:),Rvals(j,:),D.t_stim(j));
        else         
            thresh_val(j,1) = interp1(time(j,:),Lvals(j,:),D.t_stim(j));          
        end
    end
    
    subplot(length(vars),1,i)
    hist(thresh_val ./ mean(thresh_val))
    
    hold off
    ylabel('Frequency')
    xlabel(vars{i})
end


function plot_series(D,vars)

for i = 1:length(vars)
    
    eval(['Rvals = D.R' vars{i} ';']) 
    eval(['Lvals = D.L' vars{i} ';']) 
    
    if strcmp(vars{i}(end-3:end),'rate')
        time = D.time_rate;
    else
        time = D.time;
    end
        
    
    subplot(length(vars),1,i)
    
    for j = 1:size(Rvals,1)
        
        if D.R_facing(j)
            
            tmp = interp1(time(j,:),Rvals(j,:),D.t_stim(j));
            h = plot(time(j,:),Rvals(j,:),'k-',D.t_stim(j),tmp,'ro');
            set(h,'MarkerFaceColor','r')
            hold on
        else
            
            tmp = interp1(time(j,:),Lvals(j,:),D.t_stim(j));
            h = plot(time(j,:),Lvals(j,:),'k-',D.t_stim(j),tmp,'ro');
            set(h,'MarkerFaceColor','r')
            hold on
            
        end
    end
    
    hold off
    xlabel('Time (s)')
    ylabel(vars{i})
end
        

        
function [xRT,yRT,zRT] = convert_eye(xR,yR,zR,verg_ang,fov)

% psi - forward tilt angle of an eye relative to the body
psi = pi/2 + verg_ang - fov/2;

% unit vector axes
yaxis = [-sin(psi) cos(psi) 0]./norm([-sin(psi) cos(psi) 0]);
xaxis = [cos(psi) sin(psi) 0]./norm([cos(psi) sin(psi) 0]);
zaxis = [0 0 1];

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis' yaxis' zaxis'];

% Package points in a matrix
pts = [xR yR zR]';

% Rotate points
ptsT = [R * pts]';

% Extract points
xRT = ptsT(:,1);
yRT = ptsT(:,2);
zRT = ptsT(:,3);

if 1
   figure
   plot3(xRT,yRT,zRT,'.'); %axis equal
   xlabel('X'); ylabel('Y');zlabel('Z')
   view(2)  
end

        
% function [Xo,Yo,Zo] = find_bound(X,Y,Z)
% 
% for i = 1:length(Y)
%     
% end


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



function h = surf_pred(x,y,z,p_clr,light_on)

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





function [sR,sL] = eye_systems(head,tail,verg_ang,fov)

% Check dimensions
if size(head,1)~=1 || size(head,2)~=3 || size(tail,1)~=1 || size(tail,2)~=3 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = tail(1) - head(1);
xaxis(1,2) = tail(2) - head(2);
xaxis(1,3) = tail(3) - head(3);

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

% Create rotation matrix (from inertial axes to local axes)
sH = [xaxis; yaxis; zaxis];

% psi - forward tilt angle of an eye relative to the body
psi = pi/2 + verg_ang - fov/2;

% Unit vector axes for R eye (in body system)
yaxisR = [-sin(psi) cos(psi) 0]./norm([-sin(psi) cos(psi) 0]);
xaxisR = [cos(psi) sin(psi) 0]./norm([cos(psi) sin(psi) 0]);
zaxisR = [0 0 1];

% Create rotation matrix (from inertial axes to local axes)
sR = [xaxisR; yaxisR; zaxisR];

% Rotate the local system according to the eye position
sR = [sR * sH];

% Unit vector axes for L eye (in body system)
yaxisL = [-sin(psi) -cos(psi) 0]./norm([-sin(psi) cos(psi) 0]);
xaxisL = [-cos(psi) sin(psi) 0]./norm([-cos(psi) sin(psi) 0]);
zaxisL = [0 0 1];

% Create rotation matrix (from inertial axes to local axes)
sL = [xaxisL; yaxisL; zaxisL];

% Rotate the local system according to the eye position
sL = [sL * sH];




function R = prey_system(head,tail)

% Check dimensions
if size(head,1)~=1 || size(head,2)~=3 || size(tail,1)~=1 || size(tail,2)~=3 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = tail(1) - head(1);
xaxis(1,2) = tail(2) - head(2);
xaxis(1,3) = tail(3) - head(3);

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
R = [xaxis; yaxis; zaxis];

% Take inverse
%R = inv(R);


function [pXg,pYg,pZg] = prey_local(head,R,pX,pY,pZ)
% Transforms polygon points from global to local FOR

tmp1 = global_to_local(head,R,[pX(1,:)' pY(1,:)' pZ(1,:)']);
tmp2 = global_to_local(head,R,[pX(2,:)' pY(2,:)' pZ(2,:)']);
tmp3 = global_to_local(head,R,[pX(3,:)' pY(3,:)' pZ(3,:)']);
tmp4 = global_to_local(head,R,[pX(4,:)' pY(4,:)' pZ(4,:)']);

pXg = [tmp1(:,1)'; tmp2(:,1)'; tmp3(:,1)'; tmp4(:,1)'];
pYg = [tmp1(:,2)'; tmp2(:,2)'; tmp3(:,2)'; tmp4(:,2)'];
pZg = [tmp1(:,3)'; tmp2(:,3)'; tmp3(:,3)'; tmp4(:,3)'];

function [pXg,pYg,pZg] = prey_global(head,R,pX,pY,pZ)
% Transforms polygon points from local to global FOR

tmp1 = local_to_global(head,R,[pX(1,:)' pY(1,:)' pZ(1,:)']);
tmp2 = local_to_global(head,R,[pX(2,:)' pY(2,:)' pZ(2,:)']);
tmp3 = local_to_global(head,R,[pX(3,:)' pY(3,:)' pZ(3,:)']);
tmp4 = local_to_global(head,R,[pX(4,:)' pY(4,:)' pZ(4,:)']);

pXg = [tmp1(:,1)'; tmp2(:,1)'; tmp3(:,1)'; tmp4(:,1)'];
pYg = [tmp1(:,2)'; tmp2(:,2)'; tmp3(:,2)'; tmp4(:,2)'];
pZg = [tmp1(:,3)'; tmp2(:,3)'; tmp3(:,3)'; tmp4(:,3)'];


function ptsT = global_to_local(head,R,pts)

% Translate
pts(:,1) = pts(:,1) - head(1);
pts(:,2) = pts(:,2) - head(2);
pts(:,3) = pts(:,3) - head(3);

% Rotate points
ptsT = [R * pts']';



function ptsT = local_to_global(head,R,pts)

% Rotate points
ptsT = [inv(R) * pts']';

% Translate global coordinates wrt head
ptsT(:,1) = ptsT(:,1) + head(1);
ptsT(:,2) = ptsT(:,2) + head(2);
ptsT(:,3) = ptsT(:,3) + head(3);


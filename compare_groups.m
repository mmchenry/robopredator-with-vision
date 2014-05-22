function compare_groups
% Analyzes the responses of different treatment groups of larvae


%% Set paths

% Set root path
if isdir('/Users/mmchenry/Dropbox/Projects/Kelsey')
    root = '/Users/mmchenry/Dropbox/Projects/Kelsey';
    root_bill = '/Users/mmchenry/Dropbox/Projects/Robopredator';
else
    root = '/Volumes/Flow HD/Dropbox/Projects/Kelsey';
    root_bill = '/Volumes/Flow HD/Dropbox/Projects/Robopredator';
end

%% Code execution

% Plots distribution of different behavioral responses
do_behave_comp = 1;

% Plot distribution of treatment groups
do_treat_groups = 1;

% Plot distutions of LLS treated and not (with lights on)
do_LLS_treat = 1;

% 3D rendering of the distribution of larvae and predator
vis_3D_dist = 0;


%% Parameter values

% Color of 3D rendering of predator
p_clr = .5.*[1 1 1];

% Colors for the 3 speeds
sp_clr(1,:) = [0 0 1];
sp_clr(2,:) = [0 1 0];
sp_clr(3,:) = [1 0 0];

% Scale factor for spherical heads of larvae
scl_fctr = 0.0005;

% Period between LL detection and fast start initiation (s)
latency = 5e-3;


%% Organize data

% Load compiled data from Kelsey's experiments (speed2 speed11 speed20)
load([root filesep 'kelseysData.mat']);

% Load swim bladder points (swim_blad2 swim_blad22 swim_blad20)
%load('swim_bladder_pts.mat')

% Store coordinates in structure (P)
P(1).spd  = 2e-2;
P(2).spd  = 11e-2;
P(3).spd  = 20e-2;

P(1).head = speed2{1};
P(2).head = speed11{1};
P(3).head = speed20{1};

P(1).tail = speed2{2};
P(2).tail = speed11{2};
P(3).tail = speed20{2};

P(1).behav = speed2{3};
P(2).behav = speed11{3};
P(3).behav = speed20{3};

clear head_pts2 head_pts11 head_pts20 tail_pts2 tail_pts11 tail_pts20 

% Load behavior data from bill's experiments ('b')
load([root_bill filesep 'behavior' filesep 'Transformed_Prey_Coords.mat'])

% Convert units
L.preyx   = b.preyx./100;
L.preyy   = b.preyy./100;
L.preyz   = b.preyz./100;
L.preyx2  = b.preyx2./100;
L.preyy2  = b.preyy2./100;
L.preyz2  = b.preyz2./100;
L.speed   = b.speed./100;
L.lit     = b.lit;
L.LL      = b.LL;

clear b

% Number of sequences
num_seq = length(L.preyx(:,1));

% Indices for each speed of sequences in the dark, with lateral line intact
index{1} = (L.speed(1:num_seq)==2e-2) & (L.LL(1:num_seq)==1) ...
     & ~isnan(L.preyx(:,1)) & (L.lit(1:num_seq)==0);

index{2} = (L.speed(1:num_seq)==11e-2) & (L.LL(1:num_seq)==1) ...
     & ~isnan(L.preyx(:,1)) & (L.lit(1:num_seq)==0);
 
index{3} = (L.speed(1:num_seq)==20e-2) & (L.LL(1:num_seq)==1) ...
     & ~isnan(L.preyx(:,1)) & (L.lit(1:num_seq)==0);

% Load predator morphology data 
load([root filesep 'Pred3Dbodyshape.mat'])
M.x = pred3DshapeX./100;
M.y = pred3DshapeY./100;
M.z = pred3DshapeZ./100;

clear pred3DshapeX pred3DshapeY pred3DshapeZ num_seq


%% Plot distribution of different behaviors

if do_behave_comp
    
    % Index for the convex hull of the predator (frontal plane)
    iP = convhull(M.x,M.y);
    
    % New figure
    figure;
    
    % Marker size
    m_size = 2;
    
    for i = 1:length(P)
        
        subplot(2,1,1)
        
        iB = P(i).behav== 'f';
         
        hP = fill(M.x(iP),M.y(iP),0.*M.y(iP));
        set(hP,'FaceColor',p_clr,'LineStyle','none')
        hold on
        plot_larvae([P(i).head(iB,1) P(i).head(iB,2)],...
                    [P(i).tail(iB,1) P(i).tail(iB,2)],sp_clr(i,:),m_size);      
        axis equal
        xlim([-.02 .06])
        ylim([-.04 .04])
        title('Escape')
        
        subplot(2,1,2)
        
        iB = (P(i).behav== 's') | (P(i).behav== 'c');
         
        hP = fill(M.x(iP),M.y(iP),0.*M.y(iP));
        set(hP,'FaceColor',p_clr,'LineStyle','none')
        hold on
        plot_larvae([P(i).head(iB,1) P(i).head(iB,2)],...
                    [P(i).tail(iB,1) P(i).tail(iB,2)],sp_clr(i,:),m_size);      
        axis equal
        xlim([-.02 .06])
        ylim([-.04 .04])
        title('Avoidance')
 
    end
end


%% Compare distributions of escapes with and without light

if do_treat_groups


% Index for the convex hull of the predator (frontal plane)
iP = convhull(M.x,M.y);

% Index for the convex hull of the predator (sagittal plane)
iPz = convhull(M.x,M.z);

% New figure
figure;

% Marker size
m_size = 2;

for i = 1:length(P)
    
    % Choose index for sequences that have stage 2 data
    idx = index{i} & ~isnan(L.preyx(:,2)) & L.preyx(:,2)>0  & L.preyx(:,3)>0;
    
    % Bill's: Get body points for all seqnences in current speed
    [rost0,com0,tail0,rost1,com1,tail1,rost2,com2,tail2,i_vent,...
            i_wrong,prey_dir] = give_points(L,idx,latency,P(i).spd);
    
        
    
    subplot(2,2,1)
    
    iB = P(i).behav== 'f';
    
    hP = fill(M.x(iP),M.y(iP),0.*M.y(iP));
    set(hP,'FaceColor',p_clr,'LineStyle','none')
    hold on
    plot_larvae([P(i).head(iB,1) (P(i).head(iB,2))],...
                [P(i).tail(iB,1) (P(i).tail(iB,2))],sp_clr(i,:),m_size);
    axis equal
    xlim([-.02 .06])
    ylim([-.04 .04])
    title('Escapes, Lights on')
    
    
    subplot(2,2,2)
    
    iB = (P(i).behav== 's') | (P(i).behav== 'c');
    
    hP = fill(M.x(iP),M.y(iP),0.*M.y(iP));
    set(hP,'FaceColor',p_clr,'LineStyle','none')
    hold on
    plot_larvae([rost1(:,1) (rost1(:,2))],...
                [tail1(:,1) (tail1(:,2))],sp_clr(i,:),m_size);
    axis equal
    xlim([-.02 .06])
    ylim([-.04 .04])
    title('Escapes, Lights off')
    
    
    subplot(2,2,3)
    
    iB = P(i).behav== 'f';
    
    hP = fill(M.x(iPz),M.z(iPz),0.*M.z(iPz));
    set(hP,'FaceColor',p_clr,'LineStyle','none')
    hold on
    plot_larvae([P(i).head(iB,1) (P(i).head(iB,3))],...
                [P(i).tail(iB,1) (P(i).tail(iB,3))],sp_clr(i,:),m_size);
    axis equal
    xlim([-.02 .06])
    ylim([-.04 .04])
    title('Escapes, Lights on')
    
    
    subplot(2,2,4)
    
    iB = (P(i).behav== 's') | (P(i).behav== 'c');
    
    hP = fill(M.x(iPz),M.z(iPz),0.*M.z(iPz));
    set(hP,'FaceColor',p_clr,'LineStyle','none')
    hold on
    plot_larvae([rost1(:,1) (rost1(:,3))],...
                [tail1(:,1) (tail1(:,3))],sp_clr(i,:),m_size);
    axis equal
    xlim([-.02 .06])
    ylim([-.04 .04])
    title('Escapes, Lights off')
    
    
    
end

end



%% Render the data in 3D

if vis_3D_dist

figure

if 1
% Render the predator -----------------------------------------------------
for i = 1:2
    subplot(2,1,i)
    h(i) = patch(M.x,M.y,M.z,M.z*0);
    
    set(h(i),'FaceLighting','gouraud',...
        'LineStyle','none',...
        'BackFaceLighting','reverselit',...
        'FaceColor',p_clr,...
        'AmbientStrength',.5);
    
    % Lighting, & view control --------------------------------------------
    hL(1) = light('position',[0 0 20]);
    hL(2) = light('position',[0 0 -20]);
    axis tight
    axis equal
   
end
end

% Lateral view
subplot(2,1,1)
view([0 0])
xlabel('X'); ylabel('Z')
hold on
xlim([-2e-2 5e-2])
zlim([-4e-2 4e-2])
ylim([-6e-2 6e-2])

% Dorsal view
subplot(2,1,2)
view(2)
xlabel('X'); ylabel('Y')
hold on
xlim([-2e-2 5e-2])
ylim([-6e-2 6e-2])
zlim([-4e-2 4e-2])

[X,Y,Z] = sphere(30);

X = X .* scl_fctr;
Y = Y .* scl_fctr;
Z = Z .* scl_fctr;

for i = 1:length(P)
    for j = 1:length(P(i).head)
        
        if ~strcmp(P(i).behav,'n')
            
            subplot(2,1,1)
            hH = surf(X+P(i).head(j,1),Y+P(i).head(j,2),Z+P(i).head(j,3),Z.*0);
            set(hH,'FaceColor',sp_clr(i,:),'LineStyle','none')
            hT = plot3([P(i).head(j,1) P(i).tail(j,1)],...
                [P(i).head(j,2) P(i).tail(j,2)],...
                [P(i).head(j,3) P(i).tail(j,3)],'k-');
            set(hT,'Color',sp_clr(i,:));
            
            subplot(2,1,2)
            hH = surf(X+P(i).head(j,1),Y+P(i).head(j,2),Z+P(i).head(j,3),Z.*0);
            set(hH,'FaceColor',sp_clr(i,:),'LineStyle','none')
            hT = plot3([P(i).head(j,1) P(i).tail(j,1)],...
                [P(i).head(j,2) P(i).tail(j,2)],...
                [P(i).head(j,3) P(i).tail(j,3)],'k-');
            set(hT,'Color',sp_clr(i,:));
            
            bLength = sqrt((P(i).head(j,1)-P(i).tail(j,1))^2 + ...
                (P(i).head(j,2)-P(i).tail(j,2))^2 + ...
                (P(i).head(j,3)-P(i).tail(j,3))^2);
        end
    end
end

end



function  [rost0,com0,tail0,rost1,com1,tail1,rost2,com2,tail2,i_vent,...
           i_wrong,prey_dir] = give_points(b,idx,latency,spd)
% Returns the points of the prey body for all sequences denoted by 'idx'
       
% Offset in x, due to latency
lat_offset = latency * spd;

% Position (wrt predator) of body points when flow sensed
rost0 = [b.preyx(idx,1)+lat_offset b.preyy(idx,1) b.preyz(idx,1)];
com0  = [b.preyx(idx,2)+lat_offset b.preyy(idx,2) b.preyz(idx,2)];
tail0 = [b.preyx(idx,3)+lat_offset b.preyy(idx,3) b.preyz(idx,3)];

% Position (wrt predator) of body points when larva first moves
rost1 = [b.preyx(idx,1) b.preyy(idx,1) b.preyz(idx,1)];
com1  = [b.preyx(idx,2) b.preyy(idx,2) b.preyz(idx,2)];
tail1 = [b.preyx(idx,3) b.preyy(idx,3) b.preyz(idx,3)];

% Position (wrt predator) of body points at end of stage 2
rost2 = [b.preyx2(idx,1) b.preyy2(idx,1) b.preyz2(idx,1)];
com2  = [b.preyx2(idx,2) b.preyy2(idx,2) b.preyz2(idx,2)];
tail2 = [b.preyx2(idx,3) b.preyy2(idx,3) b.preyz2(idx,3)];

% Find direction of response
prey_dir(:,1) = com2(:,1) - com1(:,1);
prey_dir(:,2) = com2(:,2) - com1(:,2);
prey_dir(:,3) = com2(:,3) - com1(:,3);

% Calculate dist from midline at t1 and t2
dist1 = sqrt(com1(:,2).^2 + com1(:,3).^2);
dist2 = sqrt(com2(:,2).^2 + com2(:,3).^2);

% Index of individuals positioned ventral to predator
i_vent = com0(:,3)<=0;

% Index of responses in the 'wrong' direction
i_wrong = dist2 < dist1;


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


function h = surf_pred(x,y,z,p_clr)

% Plot surfaces
h = patch(x,y,z,z*0);

% Adjust parameters
set(h,'FaceLighting','gouraud',...
    'LineStyle','none',...
    'BackFaceLighting','reverselit',...
    'FaceColor',p_clr,...
    'AmbientStrength',.5);

% Lighting, & view control --------------------------------
hL(1) = light('position',[0 0 20]);
hL(2) = light('position',[0 0 -20]);

axis equal

xlabel('X'); ylabel('Y');zlabel('Z')
view([-180 0]);
grid on


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
R = [xaxis' yaxis' zaxis'];

% Take inverse
R = inv(R);


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


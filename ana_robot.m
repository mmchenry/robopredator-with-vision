function ana_robot
% Performs anlayses that compare the results of kelsey's experiments (with
% vision) and Bill's experiments (just lateral line).



%% Code execution

% Rose plot for elevation in global FOR 
vis_el_global = 0;

% Rose plot for azimuth in local FOR (Old way of sorting orientation)
vis_az_local = 0;

% Rose plot for azimuth in global FOR 
vis_az_global = 1;

% Group responses by binocular vision
vis_az_binoc = 1;

%Scatterplot of stimulus angle vs. azimuth escape angle
scatter_stim_az = 0;

% Visualize the position of responses (visual system)
vis_visResp = 0;

% Plots that compare the response positions with and without light
do_treat_groups = 0;

% Distance boxplots (with and without light)
do_distbox = 0;

% Exract visual statistics 
vis_stats = 0;

% Execute rose plots in global FOR
do_global = 0;

% Pool L/R responses in azimuth of responses
do_LR_pool = 1;


%% Parameters

% Number of bins used in rose plots
num_bin = 20;

% Marker Size 
mSize = 3;

% Period between LL detection and fast start initiation for the lateral line (s)
latency = 0.*5e-3;

% Color of 3D rendering of predator
p_clr = .5.*[1 1 1];

% Colors for the 3 speeds
sp_clr(1,:) = [0 0 1];
sp_clr(2,:) = [0 1 0];
sp_clr(3,:) = [1 0 0];

% Scale factor for spherical heads of larvae
scl_fctr = 0.0005;

% Deg per photorecptor pair (Easter, 1996) for 4 dpf
ren_den = 2.9235/180*pi;


%% Load data

% This should be specific to the computer that executes this code
root = '/Users/arjunnair0513/Dropbox/Shared with Arjun';

% Load Kelsey's pooled data ('L')
load([root filesep 'kelsey data' filesep 'pooled_transdata.mat'])

% Load behavior data from Bill's experiments ('b')
load([root filesep 'bills data' filesep 'Transformed_Prey_Coords.mat'])

% Convert units of Bill's data
B.preyx   = b.preyx./100;
B.preyy   = b.preyy./100;
B.preyz   = b.preyz./100;
B.preyx2  = b.preyx2./100;
B.preyy2  = b.preyy2./100;
B.preyz2  = b.preyz2./100;
B.speed   = b.speed./100;
B.lit     = b.lit;
B.LL      = b.LL;

% Load Bill's flow cue data ('r') 
load([root filesep 'bills data' filesep 'bodycue data.mat'])

% Bills: Number of sequences
num_seq = length(b.preyx(:,1));

% Bills: Values of speeds
spds = [2 11 20];

% Load predator morphology data ('M')
load([root filesep 'bills data' filesep 'Pred3Dbodyshape.mat'])
M.x = pred3DshapeX./100;
M.y = pred3DshapeY./100;
M.z = pred3DshapeZ./100;

% Bills: Indices for each speed of sequences in the dark, with lateral line intact
index{1} = (b.speed(1:num_seq)==2) & (b.LL(1:num_seq)==1) ...
     & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0);

index{2} = (b.speed(1:num_seq)==11) & (b.LL(1:num_seq)==1) ...
     & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0);
 
index{3} = (b.speed(1:num_seq)==20) & (b.LL(1:num_seq)==1) ...
     & ~isnan(b.preyx(:,1)) & (b.lit(1:num_seq)==0);

clear num_seqs

% Load visual system modeling data in 'F' structure
load([root filesep 'kelsey data' filesep 'Visual data.mat'])    


%% Visualize position of responses (vision) 

if vis_visResp

figure

clrs = {'k','k','k'};
clrs2 = {'r','r','r'};

for i = 1:3
    
    for j = 1:size(L(i).head,1)
        
        subplot(2,1,1) %-------------------------------
        h = plot([L(i).head(j,1) L(i).tail(j,1)],...
                 [L(i).head(j,2) L(i).tail(j,2)],'-');
        set(h,'Color',clrs{i})
        hold on  
        h = plot(L(i).head(j,1),L(i).head(j,2),'o');
        set(h,'Color',clrs{i})
        set(h,'MarkerFaceColor',clrs{i})
        set(h,'MarkerSize',mSize)
        
        %h = plot([L(i).head2(j,1) L(i).com2(j,1)],...
        %         [L(i).head2(j,2) L(i).com2(j,2)],'-');
        %set(h,'Color',clrs2{i}) 
        %h = plot([L(i).head2(j,1)],...
        %         [L(i).head2(j,2)],'o');
        %set(h,'Color',clrs2{i})        
        h = plot([L(i).com(j,1) L(i).com2(j,1)],...
                 [L(i).com(j,2) L(i).com2(j,2)],'-');
        set(h,'Color',clrs2{i})
        
        axis equal
        
        
        subplot(2,1,2) %-------------------------------
        
        
        h = plot([L(i).head(j,1) L(i).tail(j,1)],...
                 [L(i).head(j,3) L(i).tail(j,3)],'-');
        set(h,'Color',clrs{i})
        hold on
        
        h = plot(L(i).head(j,1),L(i).head(j,3),'o');
        set(h,'Color',clrs{i})
        set(h,'MarkerFaceColor',clrs{i})
        set(h,'MarkerSize',mSize)
        
        h = plot([L(i).com(j,1) L(i).com2(j,1)],...
                 [L(i).com(j,3) L(i).com2(j,3)],'-');
        set(h,'Color',clrs2{i})
        
        axis equal
        xlabel('X'); ylabel('Z')
    end
end

subplot(2,1,1)
xlabel('X'); ylabel('Y')

subplot(2,1,2)
xlabel('X'); ylabel('Z')

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

for i = 1:length(L)
    
    % Choose index for sequences that have stage 2 data
    idx = index{i} & ~isnan(B.preyx(:,2)) & B.preyx(:,2)>0  & B.preyx(:,3)>0;
    
    % Bill's: Get body points for all seqnences in current speed
    [rost0,com0,tail0,rost1,com1,tail1,rost2,com2,tail2,i_vent,...
            i_wrong,prey_dir] = give_points(B,idx,latency,spds(i)); 
        
    
    % Plot Kelsey's data (XY)
    subplot(2,2,1)
    
    iB = L(i).behav== 'f';
    
    hP = fill(M.x(iP),M.y(iP),0.*M.y(iP));
    set(hP,'FaceColor',p_clr,'LineStyle','none')
    hold on
    plot_larvae([L(i).head(iB,1) (L(i).head(iB,2))],...
                [L(i).tail(iB,1) (L(i).tail(iB,2))],sp_clr(i,:),m_size);
    axis equal
    xlim([-.02 .06])
    ylim([-.04 .04])
    title('Escapes, Lights on')
    
    % Plot Bill's data (XY)
    subplot(2,2,2)
    
    iB = (L(i).behav== 's') | (L(i).behav== 'c');
    
    hP = fill(M.x(iP),M.y(iP),0.*M.y(iP));
    set(hP,'FaceColor',p_clr,'LineStyle','none')
    hold on
    plot_larvae([rost1(:,1) (rost1(:,2))],...
                [tail1(:,1) (tail1(:,2))],sp_clr(i,:),m_size);
    axis equal
    xlim([-.02 .06])
    ylim([-.04 .04])
    title('Escapes, Lights off')
    
    % Plot Kelsey's data (XZ)
    subplot(2,2,3)
    
    iB = L(i).behav== 'f';
    
    hP = fill(M.x(iPz),M.z(iPz),0.*M.z(iPz));
    set(hP,'FaceColor',p_clr,'LineStyle','none')
    hold on
    plot_larvae([L(i).head(iB,1) (L(i).head(iB,3))],...
                [L(i).tail(iB,1) (L(i).tail(iB,3))],sp_clr(i,:),m_size);
    axis equal
    xlim([-.02 .06])
    ylim([-.04 .04])
    title('Escapes, Lights on')
    
    % Plot Bill's data (XZ)
    subplot(2,2,4)
    
    iB = (L(i).behav== 's') | (L(i).behav== 'c');
    
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


%% Box plots of distance for light & dark

if do_distbox
    
    % Loop thru speeds
    for i = 1:3
        
       % Extract distance data from Kelsey's 
       dHead_kel = sqrt(L(i).head(:,1).^2 + ...
                        L(i).head(:,2).^2 + ...
                        L(i).head(:,3).^2);
       dTail_kel = sqrt(L(i).head(:,1).^2 + ...
                        L(i).head(:,2).^2 + ...
                        L(i).head(:,3).^2); 
       d_kel{i} = [min([dHead_kel dTail_kel]')]';
       
       % Extract distance from Bill's
       idx = B.speed.*100 == spds(i);
       dHead_bill = sqrt(B.preyx(idx,1).^2 + ...
                         B.preyy(idx,1).^2 + ...
                         B.preyz(idx,1).^2);
       dTail_bill = sqrt(B.preyx(idx,3).^2 + ...
                         B.preyy(idx,3).^2 + ...
                         B.preyz(idx,3).^2);
       d_bill{i} = [min([dHead_bill dTail_bill]')]';
       
       clear dHead_kel dTail_kel dHead_bill dTail_bill
    end
    
    
    % Make figure
    figure
    
    % Loop thru speeds
    for i = 1:3
        
       % Plot left panel (lights on)
       h(1) = subplot(3,2,(i-1).*2+1);
       boxplot(d_kel{i}.*1000)
       title([num2str(spds(i)) ' cm/s, Lights on'])
       ylabel('Distance (mm)')
       ylims(1,:) = ylim;
       
       % Plot right panel (lights off)
       h(2) = subplot(3,2,(i-1).*2+2);
       boxplot(d_bill{i}.*1000)
       title([num2str(spds(i)) ' cm/s, Lights off'])
       ylabel('Distance (mm)')
       ylims(2,:) = ylim;
       
       % Set same limits
       axes(h(1))
       ylim([min(ylims(:,1)) max(ylims(:,2))])
       
       axes(h(2))
       ylim([min(ylims(:,1)) max(ylims(:,2))])
    end
    
    clear d_bill d_kel
end


%% Calc stats on visual stimulus, pool all data, store in 'D'
   
if vis_stats
n = 1;

% Sum of all larvae
all_larvae = length(F(1).head) + length(F(2).head) + length(F(3).head);

% Loop thru approach speeds
for i = 1:length(F)
   
    % Loop thru individuals
    for j = 1:size(F(i).head,1)
         
        % Loop thru images for each instant of time
        for k = 1:size(F(i).imR,4)
            
            % Time 
            time(k,1) = F(i).time(k,j);
            
            % Stimulated area (deg^2)
            Rarea(k,1) = sum(sum(F(i).imR(:,:,j,k))).*(ren_den/pi*180)^2;
            Larea(k,1) = sum(sum(F(i).imL(:,:,j,k))).*(ren_den/pi*180)^2;
            
            % Angular position of the area's centroid (Right)
            [iX,iY] = find(F(i).imR(:,:,j,k));
            RazCent(k,1)  = mean(mean(F(i).R_azField(iX,iY)));
            RelCent(k,1)  = mean(mean(F(i).R_elField(iX,iY)));
            RangDist(k,1) = sqrt(RazCent(k,1)^2 + RelCent(k,1)^2);
            
            % Angular position of the area's centroid (Left)
            [iX,iY] = find(F(i).imL(:,:,j,k));
            LazCent(k,1)  = mean(mean(F(i).L_azField(iX,iY)));
            LelCent(k,1)  = mean(mean(F(i).L_elField(iX,iY)));
            LangDist(k,1) = sqrt(LazCent(k,1)^2 + LelCent(k,1)^2);
            
            % Angular range (deg) of stimulus in az and el
            RazRange(k,1) = max([0 range(F(i).R_azField(F(i).imR(:,:,j,k)))/pi*180]);
            RelRange(k,1) = max([0 range(F(i).R_elField(F(i).imR(:,:,j,k)))/pi*180]);
            LazRange(k,1) = max([0 range(F(i).L_azField(F(i).imL(:,:,j,k)))/pi*180]);
            LelRange(k,1) = max([0 range(F(i).L_elField(F(i).imL(:,:,j,k)))/pi*180]);
            
            clear iX iY
        end  
        
        % Store sequence data
        D.spd(n,1)      = F(i).spd;
        D.pred_pos(n,:) = F(i).pred_pos(:,j)';
        D.t_resp(n,1)   = F(i).t_resp(j);
        D.behav(n,1)    = L(i).behav(j);
        D.az(n,1)       = L(i).az(j);
        D.el(n,1)       = L(i).el(j);
        D.azL(n,1)      = L(i).azL(j);
        D.elL(n,1)      = L(i).elL(j);
        D.head(n,:)     = L(i).head(j,:);
        D.tail(n,:)     = L(i).tail(j,:);
        D.com(n,:)      = L(i).com(j,:);
        D.com2(n,:)     = L(i).com2(j,:);
        D.wrong(n,1)    = L(i).wrong(j);
        
        % Index for time interval to be considered
        idx = (time < -.1) & (time > -.3);
        
        % Logical for whether right eye got the bigger area stim
        D.R_facing(n,1) = max(Rarea(idx)) > max(Larea(idx));
        
        % Logical for binocular stimulus
        D.binoc(n,1)    = sum(Rarea(idx))~=0 & sum(Larea(idx))~=0;
        
        % Retinal position of predator
        D.RazCent(n,1) = mean(RazCent(idx));
        D.RelCent(n,1) = mean(RelCent(idx));
        D.LazCent(n,1) = mean(LazCent(idx));
        D.LelCent(n,1) = mean(LelCent(idx));
        
        % Update status
        disp(['Done ' num2str(n) ' of ' num2str(all_larvae)])
        
        % Clear for next iteration
        clear time Rarea Larea RazCent RelCent LazCent LelCent RazRange RelRange LazRange
        clear LelRange
        
%         % Find first derivative of each parameter
%         D(i).time_rate(j,:)     = D(i).time(j,1:(end-1)) + mean(diff(D(i).time(j,:)))/2;
%         D(i).Rarea_rate(j,:)    = diff(D(i).Rarea(j,:)) ./ diff(D(i).time(j,:));
%         D(i).Larea_rate(j,:)    = diff(D(i).Larea(j,:)) ./ diff(D(i).time(j,:));
%         D(i).RazRange_rate(j,:) = diff(D(i).RazRange(j,:))./ diff(D(i).time(j,:));
%         D(i).RelRange_rate(j,:) = diff(D(i).RelRange(j,:))./ diff(D(i).time(j,:));
%         D(i).LazRange_rate(j,:) = diff(D(i).LazRange(j,:))./ diff(D(i).time(j,:));
%         D(i).LelRange_rate(j,:) = diff(D(i).LelRange(j,:))./ diff(D(i).time(j,:));
       
        % Advance index
        n = n + 1;
    end
end

% Save processed pooled data ('D')
save([root filesep 'pooled_processed.mat'],'D')

disp(' ');disp(' ');disp(' ');disp(' ')

else
    % Load processed pooled data ('D')
    load([root filesep 'pooled_processed.mat'])
end


%% Visulize rose plots of elevation wrt speed (vis_el_global)

if vis_el_global
    
    clrs1 = {'k','k','k'};
    clrs2 = {'b','b','b'};
    clrs3 = {'r','r','r'};

    figure
    
    % Loop thru speeds
    for i = 1:3
        
        % Index for current speed
        idx = (D.behav=='f') & (D.spd==spds(i));
        
        % Index for downward elevation
        idx1 = (D.RelCent < 0) | (D.LelCent < 0) ;
           
        % Loop thru individuals
        for j = 1:length(idx)
            if idx(j)
                
                % Plot position at time of response & resp direction
                subplot(3,3,i)
                h = plot([D.head(j,1) D.tail(j,1)],...
                         [D.head(j,3) D.tail(j,3)],'-');
                if idx1(j)
                    set(h,'Color',clrs2{i})
                else
                    set(h,'Color',clrs3{i})
                end
                
                hold on
                h = plot(D.head(j,1),D.head(j,3),'o');
                
                if idx1(j)
                    set(h,'Color',clrs2{i})
                    set(h,'MarkerFaceColor',clrs2{i})
                else
                    set(h,'Color',clrs3{i})
                    set(h,'MarkerFaceColor',clrs1{i})
                end
                
                set(h,'MarkerSize',mSize)
                
                h = plot([D.com(j,1) D.com2(j,1)],...
                         [D.com(j,3) D.com2(j,3)],'-');
                set(h,'Color',clrs1{i})
            end
        end
        
        % Axis labels for position plots
        axis equal; xlabel('X'); ylabel('Z');
        
        % Title
        if do_global
            title(['Global FOR, ' num2str(spds(i)) 'cm/s'])
        else
            title(['Local FOR, ' num2str(spds(i)) 'cm/s'])
        end
        
        % Rose plot for when predator is below
        subplot(3,3,i+3)
        if do_global
            rose_plot(D.el(idx & idx1),D.wrong(idx & idx1),num_bin)  
        else
            rose_plot(D.elL(idx & idx1),D.wrong(idx & idx1),num_bin)
        end
        
        h = title('Predator below (el)');
        set(h,'Color',clrs2{i})
        
        % Rose plot for when predator is above
        subplot(3,3,i+6)
        if do_global
            rose_plot(D.el(idx & ~idx1),D.wrong(idx & ~idx1),num_bin)
        else
            rose_plot(D.elL(idx & ~idx1),D.wrong(idx & ~idx1),num_bin)
        end
        
        h = title('Predator above (el)');
        set(h,'Color',clrs3{i})
 
    end
end % vis_el_global


%% Visulize rose plots of azimuth wrt speed (vis_az_local)

if vis_az_local
    
    clrs1 = {'k','k','k'};
    clrs2 = {'b','b','b'};
    clrs3 = {'r','r','r'};

    figure
    
    % Loop thru speeds
    for i = 1:3       
        
        % Index for current speed
        idx = (D.behav=='f') & (D.spd==spds(i));
        
        % Index for azimuth more on center for the right
        idx1 = (~isnan(D.RazCent) & isnan(D.LazCent)) | ...
               (~isnan(D.RazCent) & (abs(D.RazCent) < abs(D.LazCent)));
        
        % Index for azimuth more on center for the left
        idx2 = (~isnan(D.LazCent) & isnan(D.RazCent)) | ...
               (~isnan(D.RazCent) & (abs(D.LazCent) < abs(D.RazCent)));
           
        % Loop thru individuals
        for j = 1:length(idx)
            if idx(j)
                
                % Plot position at time of response & resp direction 
                if do_LR_pool
                    subplot(2,3,i)
                else
                    subplot(3,3,i)
                end
                
                h = plot([D.head(j,1) D.tail(j,1)],...
                         [D.head(j,2) D.tail(j,2)],'-');
                     
                if idx1(j)
                    set(h,'Color',clrs2{i})
                elseif idx2(j)
                    set(h,'Color',clrs3{i})
                else
                    set(h,'Color',.8.*[1 1 1])
                end
                
                hold on
                h = plot(D.head(j,1),D.head(j,2),'o');
                
                if idx1(j)
                    set(h,'Color',clrs2{i})
                    set(h,'MarkerFaceColor',clrs2{i})
                elseif idx2(j)
                    set(h,'Color',clrs3{i})
                    set(h,'MarkerFaceColor',clrs3{i})
                else
                    set(h,'Color',.8.*[1 1 1])
                end
                
                set(h,'MarkerSize',mSize)
                
                h = plot([D.com(j,1) D.com2(j,1)],...
                         [D.com(j,2) D.com2(j,2)],'-');
                set(h,'Color',clrs1{i})
            end
        end
        
        axis equal; xlabel('X'); ylabel('Y');

        % Title
        title(['Local FOR, ' num2str(spds(i)) 'cm/s'])
        
        % Rose plots with pooled responses
        if do_LR_pool
            
            wrong_pool = [D.wrong(idx & idx1); D.wrong(idx & idx2)];
            
            % Pool responses (flip response of right eye stim)
            az_pool    = [-D.azL(idx & idx1); D.azL(idx & idx2)];
            
            subplot(2,3,i+3)
            rose_plot(az_pool,wrong_pool,num_bin)
            title('Response away from stim')
            
        % Rose plots with L and R separated
        else
            % Rose plot for when predator is on right
            subplot(3,3,i+3)
            rose_plot(D.azL(idx & idx1),D.wrong(idx & idx1),num_bin)
            
            h = title('Predator on right (az)');
            set(h,'Color',clrs2{i})
            
            % Rose plot for when predator on left
            subplot(3,3,i+6)
            rose_plot(D.azL(idx & idx2),D.wrong(idx & idx2),num_bin)
            
            h = title('Predator on left (az)');
            set(h,'Color',clrs3{i})
        end
    end
    
    clear az_pool
end % vis_az_local

%% Visulize rose plots of azimuth wrt speed (vis_az_global)

if vis_az_global
    
    clrs1 = {'k','k','k'};
    clrs2 = {'b','b','b'};
    clrs3 = {'r','r','r'};

    figure
    
    % Loop thru speeds
    for i = 1:3       
        
        % Index for current speed
        idx = ((D.behav=='f') | (D.behav=='s')) & (D.spd==spds(i));
        
        % Index for larvae that are left of the predator
        idx1 = D.com(:,2) >= 0;
        
        % Index for larvae that are right of the predator
        idx2 = ~idx1;
           
        % Loop thru individuals
        for j = 1:length(idx)
            if idx(j)
                
                % Plot position at time of response & resp direction 
                if do_LR_pool
                    subplot(2,3,i)
                else
                    subplot(3,3,i)
                end
                
                h = plot([D.head(j,1) D.tail(j,1)],...
                         [D.head(j,2) D.tail(j,2)],'-');
                     
                if idx1(j)
                    set(h,'Color',clrs2{i})
                elseif idx2(j)
                    set(h,'Color',clrs3{i})
                else
                    set(h,'Color',.8.*[1 1 1])
                end
                
                hold on
                h = plot(D.head(j,1),D.head(j,2),'o');
                
                if idx1(j)
                    set(h,'Color',clrs2{i})
                    set(h,'MarkerFaceColor',clrs2{i})
                elseif idx2(j)
                    set(h,'Color',clrs3{i})
                    set(h,'MarkerFaceColor',clrs3{i})
                else
                    set(h,'Color',.8.*[1 1 1])
                end
                
                set(h,'MarkerSize',mSize)
                
                h = plot([D.com(j,1) D.com2(j,1)],...
                         [D.com(j,2) D.com2(j,2)],'-');
                set(h,'Color',clrs1{i})
            end
        end
        
        axis equal; xlabel('X'); ylabel('Y');

        % Title
        title(['Global FOR, ' num2str(spds(i)) 'cm/s'])
        
        % Rose plots with pooled responses
        if do_LR_pool
            
            wrong_pool = [D.wrong(idx & idx1); D.wrong(idx & idx2)];
            
            % Pool responses (flip response of larvae on the right side)
            az_pool    = [D.az(idx & idx1); -D.az(idx & idx2)];
            
            subplot(2,3,i+3)
            rose_plot(az_pool,wrong_pool,num_bin)
            title('Response away from stim')
            
        % Rose plots with L and R separated
        else
            % Rose plot for when predator is on right
            subplot(3,3,i+3)
            rose_plot(D.az(idx & idx1),D.wrong(idx & idx1),num_bin)
            
            h = title('Predator on right (az)');
            set(h,'Color',clrs2{i})
            
            % Rose plot for when predator on left
            subplot(3,3,i+6)
            rose_plot(D.az(idx & idx2),D.wrong(idx & idx2),num_bin)
            
            h = title('Predator on left (az)');
            set(h,'Color',clrs3{i})
        end
    end
    
    clear az_pool
end % vis_az_global


%% Visulize rose plots of azimuth wrt binocular stim (vis_az_binoc)

if vis_az_binoc
    
    % Create figure
    figure
    
    % Index for binocular stimulus
    idxB = D.binoc & (D.behav=='s');
    
    % Index for left eye
    idxL = ~D.binoc & ~isnan(D.LazCent) & (D.behav=='s');
    
    % Index for right eye
    idxR = ~D.binoc & ~isnan(D.RazCent) & (D.behav=='s');
    
    % Loop thru individuals
    for j = 1:length(D.spd)
        if idxB(j)
            subplot(2,3,1)
            
        elseif idxL(j)
            subplot(2,3,2)
        elseif idxR(j)
            subplot(2,3,3)
        end
        
        h = plot([D.head(j,1) D.tail(j,1)],...
            [D.head(j,2) D.tail(j,2)],'-k');
        hold on
        h = plot(D.head(j,1),D.head(j,2),'ok');
        set(h,'Color','k')
        set(h,'MarkerFaceColor','k')
        set(h,'MarkerSize',mSize)
        h = plot([D.com(j,1) D.com2(j,1)],...
            [D.com(j,2) D.com2(j,2)],'-');
        set(h,'Color','r')
    end
    
    axis equal; xlabel('X'); ylabel('Y');
    
    % Title
    if do_global
        title(['Global FOR'])
    else
        title(['Local FOR'])
    end
    
    
    % Rose plots (global FOR)
    if do_global
        subplot(2,3,4)
        rose_plot(D.az(idxB),D.wrong(idxB),num_bin)
        title('Binocular stimulus')
        
        subplot(2,3,5)
        rose_plot(D.az(idxL),D.wrong(idxL),num_bin)
        title('Left eye stimulus')
        
        subplot(2,3,6)
        rose_plot(D.az(idxR),D.wrong(idxR),num_bin)
        title('Right eye stimulus')
        
        % Rose plots (local FOR)
    else
        subplot(2,3,4)
        rose_plot(D.azL(idxB),D.wrong(idxB),num_bin)
        title('Binocular stimulus')
        
        subplot(2,3,5)
        rose_plot(D.azL(idxL),D.wrong(idxL),num_bin)
        title('Left eye stimulus')
        
        subplot(2,3,6)
        rose_plot(D.azL(idxR),D.wrong(idxR),num_bin)
        title('Right eye stimulus')
        
    end
    
    % Scatterplots of response wrt stimulus -------------------------------
    if scatter_stim_az
        % Make figure window
        figure
        
        subplot(2,2,1)
        plot(D.LazCent(idxB),D.az(idxB),'or', D.LazCent(idxL),D.az(idxL),'ob')
        xlabel('angle of stimulus')
        ylabel('Azimuth of resposne (global)')
        title('Left Eye: GLOBAL')
        legend('Binoc','One eye')
        grid on
        axis equal
        
        subplot(2,2,3)
        plot(D.LazCent(idxB),D.azL(idxB),'or', D.LazCent(idxL),D.azL(idxL),'ob')
        xlabel('angle of stimulus')
        ylabel('Azimuth of response (local)')
        title('Left Eye: LOCAL')
        grid on
        axis equal
        
        subplot(2,2,2)
        plot(D.RazCent(idxB),D.az(idxB),'or',D.RazCent(idxR),D.az(idxR),'ob')
        xlabel('angle of stimulus')
        ylabel('Azimuth of resposne (global)')
        title('Right Eye: GLOBAL')
        grid on
        axis equal
        
        subplot(2,2,4)
        plot(D.RazCent(idxB),D.azL(idxB),'or',D.RazCent(idxR),D.azL(idxR),'ob')
        xlabel('angle of stimulus')
        ylabel('Azimuth of response (local)')
        title('Right Eye: LOCAL')
        grid on
        axis equal
    end
    
    clear idxB idxR idxL
end % vis_az_global



function rose_plot(ang,wrong,num_bin)
% Creates customized rose plots (includes 95% confidence intervals)  

% Define 'wrong' index
wrong = logical(wrong);

% Calculate and plot mean and CIs
warning off
[mu,l1,l2] = circ_mean(ang);
warning on

% Rose plot in correct direction
h = rose(ang(~wrong),num_bin);
set(h,'Color','g')
hold on

% Rose plot of wrong direction
if sum(wrong)>0
    h = rose(ang(wrong),num_bin);
    set(h,'Color','r')
end

% Find axis limits
ymax = abs(max(ylim));

% Create series of values between CIs
r_val = linspace(l1,l2,50);

% Create polar plots: Mean
%h = polar([mu mu],[0 ymax]);
%set(h,'Color','k')

% Create polar plots: CIs
h = polar(r_val,[0 ymax.*ones(1,length(r_val)-2) 0]);
set(h,'Color',.5.*[1 1 1])
set(h,'LineStyle','-')
hold off



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
% Plots larvae

if size(head,2) == 2
    for i = 1:size(head,1)      
        h = plot(head(i,1),head(i,2),'o',[head(i,1) tail(i,1)],...
                 [head(i,2) tail(i,2)],'-');
        set(h(1),'MarkerFaceColor',clr,'MarkerEdgeColor',clr,...
                 'MarkerSize',msize)
        set(h(2),'Color',clr)           
    end
end


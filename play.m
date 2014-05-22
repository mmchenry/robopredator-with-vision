function play


qrypts=[xgridcoords, ygridccords]; %grid x and y coordinates
triids = pointLocation(dt, qrypts); %associate each grid point to its Delaunay triangle
Lia = ismember(triids,dtsubset); %keep subset of points contained in the desired triangles (dtsubset contains the indices of desired triangles)
IM=false(size(grid)); 
IM(Lia)=1;
function [collision,errPose]=postAnalyze(obj,g_sum)
    %Test if the whole manipulator collide with an obstacle.
    %Also, output the position and orientation of the topmost plate

    collision=false;    
    rDiff=3; % connecting point of the vertical link to the plate is not located exactly at the vertex of the plate.
    protrusionPlate=2; % height of protrusion on the plates. 
    thicknessPlate=1.5;  %1.5mm, thickness of the plates.
    g_protrusion=[0,0,0,0,0,protrusionPlate]'; %transformation vector for the protrusion itself (considering the height of the protrusion).
    g_plate=[0,0,0,0,0,thicknessPlate]'; %transformation vector for the plate itself (considering the thickness of the plate).
    if nargin<2
        g_sum=obj.G_vector(:,end);
    end
    
    base_plate = obj.r_vec(obj.r+rDiff);    %3X3 matrix, each column represents the xyz coordinates of a vertex.
    base_plate_vlink=obj.r_vec();
    
    %-----------------------------------for the bottom most plate (including its upward protrusions)
    top_plate=base_plate;
    test_plate_collision(top_plate); 
    g_top_matrix=Modules.multi_T(eye(4),g_plate);    
    top_plate_vlink_up_protrusion_bot=Modules.transform(g_top_matrix,base_plate_vlink);
    g_top_matrix=Modules.multi_T(g_top_matrix,g_protrusion);
    top_plate_vlink_up_protrusion_top=Modules.transform(g_top_matrix,base_plate_vlink);
    for i=1:obj.N
        two_vertices=[top_plate_vlink_up_protrusion_bot(:,i),top_plate_vlink_up_protrusion_top(:,i)]; %3X2 matrix
        test_link_collision(two_vertices);
    end
    %-----------------------------------for the bottom most plate (including its upward protrusions)
    
    for jj=1:obj.num
       g=g_sum(((jj-1)*6+1):6*jj);
       %g_matrix is the corresponding transformation matrix of vector g
       g_top_matrix=Modules.multi_T(g_top_matrix,g);
    
       %testing the collision of the vertical links with the obstacle
       top_plate_vlink_down_protrusion_bot=Modules.transform(g_top_matrix,base_plate_vlink);
       for ij=1:obj.N
          test_link_collision([top_plate_vlink_up_protrusion_top(:,ij),top_plate_vlink_down_protrusion_bot(:,ij)]); 
       end  
       %testing the collision of the vertical links with the obstacle
    
       %testing the collision of the protrusions and plates with the obstacle
       g_top_matrix=Modules.multi_T(g_top_matrix,g_protrusion);
       top_plate_vlink_down_protrusion_top=Modules.transform(g_top_matrix,base_plate_vlink);
       top_plate=Modules.transform(g_top_matrix,base_plate);
       for ij=1:obj.N
          test_link_collision([top_plate_vlink_down_protrusion_bot(:,ij),top_plate_vlink_down_protrusion_top(:,ij)]); 
       end
       test_plate_collision(top_plate); 
    
       if ~ (jj==obj.num)
           g_top_matrix=Modules.multi_T(g_top_matrix,g_plate);   
           top_plate_vlink_up_protrusion_bot=Modules.transform(g_top_matrix,base_plate_vlink);
           g_top_matrix=Modules.multi_T(g_top_matrix,g_protrusion);
           top_plate_vlink_up_protrusion_top=Modules.transform(g_top_matrix,base_plate_vlink);
           for ij=1:obj.N
              test_link_collision([top_plate_vlink_up_protrusion_bot(:,ij),top_plate_vlink_up_protrusion_top(:,ij)]);            
           end
       else
           g_top_matrix=Modules.multi_T(g_top_matrix,g_plate);  %this is the output for the top most plate.
       end
    end
    
%     rotationMatrix=Modules.R(g_top_matrix(1:3));
%     errOrientation=100*(rotationMatrix(1,3)^2+rotationMatrix(2,3)^2);
    %we expect the Z axis of the coordinate system attached to the top plate is perpendicular to the fixed bottom most plate.
    %In this case, the last column of rotationMatrix should be [0,0,-1];
    errPosition=g_top_matrix(1:3,4)-[-35,0,0]';
    %errPose=norm(errPosition);  %as LI Chen suggested, we do not calculate the distance in the matlab code.
    errPose=errPosition;

    function test_plate_collision(planar_plate)
        % This function is used to determine if the planar plate collide with obstacle
        isCollision=Modules.Face2FaceCollision(obj.obstacle,planar_plate);
        if isCollision
           collision=1;
        end
    end

    function test_link_collision(two_vertices)
        % This function is used to determine if the axis (cylinder) collide with obstacle
        isCollision=Modules.Ray2FaceCollision(obj.obstacle,two_vertices);
        if isCollision
           collision=1;
        end
    end

end      

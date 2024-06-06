function modulesPlot(obj,g_sum)
    %plot the module in 3 dimension
    clf;
    rDiff=3; % connecting point of the vertical link to the plate is not located exactly at the vertex of the plate.
    protrusionPlate=2; % height of protrusion on the plates. 
    thicknessPlate=1.5;  %1.5mm, thickness of the plates.
    g_protrusion=[0,0,0,0,0,protrusionPlate]'; %transformation vector for the protrusion itself (considering the height of the protrusion).
    g_plate=[0,0,0,0,0,thicknessPlate]'; %transformation vector for the plate itself (considering the thickness of the plate).
    radiusLinks=1; %0.8mm, radius of the protrusion and the links.
    do_plot=true;
    if nargin<2
        g_sum=obj.G_vector(:,end);
    end
    view([30,45]);
    hold on
    get_plate(obj.obstacle,thicknessPlate,do_plot,'m',1); 
    
    base_plate = obj.r_vec(obj.r+rDiff);    %3X3 matrix, each column represents the xyz coordinates of a vertex.
    base_plate_vlink=obj.r_vec();
    
    %-----------------------------------for the bottom most plate (including its upward protrusions)
    top_plate=base_plate;
    get_plate(top_plate,thicknessPlate,do_plot,[0.2,0.2,0.2],obj.faceAlpha_plates); %plot the base plate (bottom plate of the bottom module)
    g_top_matrix=Modules.multi_T(eye(4),g_plate);    
    top_plate_vlink_up_protrusion_bot=Modules.transform(g_top_matrix,base_plate_vlink);
    g_top_matrix=Modules.multi_T(g_top_matrix,g_protrusion);
    top_plate_vlink_up_protrusion_top=Modules.transform(g_top_matrix,base_plate_vlink);
    for i=1:obj.N
        two_vertices=[top_plate_vlink_up_protrusion_bot(:,i),top_plate_vlink_up_protrusion_top(:,i)]; %3X2 matrix
        get_link(two_vertices,radiusLinks,do_plot,[0.2,0.2,0.2],obj.faceAlpha_vlinks);
    end
    %-----------------------------------for the bottom most plate (including its upward protrusions)
    
    for jj=1:obj.num
       g=g_sum(((jj-1)*6+1):6*jj);
       %g_matrix is the corresponding transformation matrix of vector g
       g_top_matrix=Modules.multi_T(g_top_matrix,g);
    
       %plotting the vertical links  
       top_plate_vlink_down_protrusion_bot=Modules.transform(g_top_matrix,base_plate_vlink);
       for ij=1:obj.N
          get_link([top_plate_vlink_up_protrusion_top(:,ij),top_plate_vlink_down_protrusion_bot(:,ij)],radiusLinks,do_plot,[0.753,0.753,0.753],obj.faceAlpha_vlinks); 
          %[0.753,0.753,0.753] is the rgb set for silver
       end  
       %plotting the vertical links
    
       %plotting the bottom protrusion of the top plate and the top plate itself
       g_top_matrix=Modules.multi_T(g_top_matrix,g_protrusion);
       top_plate_vlink_down_protrusion_top=Modules.transform(g_top_matrix,base_plate_vlink);
       top_plate=Modules.transform(g_top_matrix,base_plate);
       for ij=1:obj.N
          get_link([top_plate_vlink_down_protrusion_bot(:,ij),top_plate_vlink_down_protrusion_top(:,ij)],radiusLinks,do_plot,[0.2,0.2,0.2],obj.faceAlpha_vlinks); 
       end
       get_plate(top_plate,thicknessPlate,do_plot,[0.2,0.2,0.2],obj.faceAlpha_plates); 
       %plotting the bottom protrusion of the top plate and the top plate itself       

       %plotting the top protrusion of the top plate
       if ~ (jj==obj.num)
           g_top_matrix=Modules.multi_T(g_top_matrix,g_plate);   
           top_plate_vlink_up_protrusion_bot=Modules.transform(g_top_matrix,base_plate_vlink);
           g_top_matrix=Modules.multi_T(g_top_matrix,g_protrusion);
           top_plate_vlink_up_protrusion_top=Modules.transform(g_top_matrix,base_plate_vlink);
           for ij=1:obj.N
              get_link([top_plate_vlink_up_protrusion_bot(:,ij),top_plate_vlink_up_protrusion_top(:,ij)],radiusLinks,do_plot,[0.2,0.2,0.2],obj.faceAlpha_vlinks);            
           end
       else
           g_top_matrix=Modules.multi_T(g_top_matrix,g_plate);  %this is the output for the top most plate.
       end
       %plotting the top protrusion of the top plate
    end
    
    axis equal
    axis off;
    hold off;


    function thick_plate=get_plate(planar_plate,thickness,do_plot,color,face_alpha)
        % This function is used to generate the data for a 3-d plate (a block), based on the given data for a planar plate.
        % Basically, the 2-d plate extends by thickness perpendicularly to the upward direction.
        %planar-plate-------3X3 matrix. Each column represents a vertex.
        %thickness-----------scalar
        %do_plot-----------boolean, indicating whether or not plot the plate using patch
        %thick_plate---------3X6 matrix.
       
        %if the plate collide with the obstacle (surface), turn the color to red.
        if size(planar_plate,2)==3  %only applies on the triangular plates, not on the obstacle itself
           isCollision=Modules.Face2FaceCollision(obj.obstacle,planar_plate);
           if isCollision
               color='red';
           end
        end

        normal=obj.get_normal(planar_plate);
        thick_plate=[planar_plate,planar_plate+thickness.*normal];
        if do_plot
            if size(planar_plate,2)==3
               faces_bot=[1,2,3,NaN;4,5,6,NaN;1,4,5,2;2,5,6,3;3,6,4,1];
               patch('Faces',faces_bot,'Vertices',thick_plate','EdgeColor','none','FaceColor',color,'FaceAlpha',face_alpha); 
            elseif size(planar_plate,2)==4
               faces_bot=[1,2,3,4;5,6,7,8;1,5,6,2;2,6,7,3;3,7,8,4;4,8,5,1];
               patch('Faces',faces_bot,'Vertices',thick_plate','EdgeColor','b','FaceColor',color,'FaceAlpha',face_alpha); 
            end
        end
    end

    function thick_link=get_link(two_vertices,radius,do_plot,color,face_alpha)
        % This function is used to generate the data for a cylindrical vertical link, 
        %             based on the given data for an axis (determined by two vertices).
        %two_vertices-------3X2 matrix. Each column represents a vertex (the first column is the starting point of the vector).
        %radius-----------scalar
        %thick_link---------3X2n matrix.

        %if the link collide with the obstacle (surface), turn the color to red.
        isCollision=Modules.Ray2FaceCollision(obj.obstacle,two_vertices);
        if isCollision
            color='red';
        end

        n=40; %40 points along the circumference
        axis=two_vertices(:,2)-two_vertices(:,1);
        opts = optimset('Display','off');
        x=fsolve(@(x)f_circle(x,axis,radius),[0,1]',opts); 
        initial_dot=[x(:);0]; % a column vector. the dot is located on the bottom plate, whose center of circle is located at the origin
        circular_dots=zeros(3,2*n);
        for i=1:n
            axang=[axis(:)',(i-1)*2*pi/n]; %a row vector
            rotm=axang2rotm(axang);
            circular_dots(:,i)=rotm*initial_dot+two_vertices(:,1);
        end
        circular_dots(:,n+1:2*n)=circular_dots(:,1:n)+axis;
        thick_link=circular_dots;

        if do_plot
            faces=[1:n;n+1:2*n];
            patch('Faces',faces,'Vertices',thick_link','EdgeColor','none','FaceColor',color,'FaceAlpha',face_alpha);    
            %plot the bottom and top surfaces of the cylinder (vertical link)

            faces2=zeros(n,4);
            for j=1:n-1
               faces2(j,:)=[j,j+1,j+n+1,j+n];
            end
            faces2(n,:)=[n,1,n+1,2*n];
            patch('Faces',faces2,'Vertices',thick_link','EdgeColor','none','FaceColor',color,'FaceAlpha',face_alpha);
            %plot the cylindrical body of the cylinder (vertical link)
        end
    end

    function F=f_circle(x,axis,radius)
        %a special case is that the circle is parallel to the xy plane (or axis is perpendiculat to the xy plane).
        F(1)=x(1)^2+x(2)^2-radius^2;
        F(2)=x(1)*axis(1)+x(2)*axis(2);
    end
end      
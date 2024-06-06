classdef Modules < handle
     %in this version of the modules, we consider the thickness of the plates, height of the protrusions. 
     %Also, both clockwise and counterclockwise modules are considered. We note that bs always sticks to the numbering of the face.
     %as sticks with the vertical link. So for a counterclockwise module, as(1) and bs(1) are in face(1). 
     %for a clockwise module, as(2) and bs(1) are in face(1). as(1) pairs with bs(3) instead, in face(3).
     properties
         Stiffnesses (:,1) double {mustBeNonnegative} = 1
         Cable_displacements (3,1) double {mustBeNonnegative} = 1
         r=17.8979  %distance from the center to the vertex of the triangular plate
         h (1,1)    %height of the module (distance between the top plate and the bottom plate)
         k (1,1)    %ratio of height to dimension of triangular plate
         phi=10/180*pi     %pre-rotation of the modules
         l (1,1)         %refer to page 2 of the slides.
         N (1,1)=3               %number of sides of plolygon. In the case of triangular plates, N=3;
         num (1,1) =4     %number of modules
         stiffness_plateJoints=0.5;  % the stiffness of the 30A silicone tube acting as the soft joints on bottom and top plates
         b_ini (1,1)      %refer to page 2 of the slides.
         a_ini (1,1)      %refer to page 2 of the slides.
         as (3,1)         %as, a 3X1 column vector, denotes the lengths of vertical links for a single module
         As (:,1)         %As, a num*3 column vector, denotes the lengths of vertical links for all module stacked in series
         Bs (:,1)         %Bs ,a num*3 column vector, is varying throught the deformation. 
                          % It denotes the lengths of diagonal links for all module stacked in series
         base_vertices (3,3) %coordinates of the three vertices on the bottom plates, for a single module          
         Bs_vector        %matrix used to store the values of Bs for every frame/iteration of the simulation.
         G_vector         %store G (for multiple modules stacked in series) for every frame/iteration of the simulation.
         g_straight (6,1)  %value of g (for a single module in its straight/initial shape).
         Fval_vector  (1,:)  %store fval (for multiple modules stacked in series) for every frame/iteration of the simulation.
         Thetas_vector    %store Thetas (for multiple modules stacked in series) for every frame/iteration of the simulation. 
         count (1,1)  =1
         n_iter (1,1) =3 %number of iterations of the forward calculation, default to 1. 
         b_del (3,1)    %step size of each cable displacement. if n_iter=1, then b_del is equal to cable_displacements.  
         display_error =false  %in the step method, if display_error is true, then display error if any.
         isSuccessful =false %if the forward model runs smoothly without any errors

         faceAlpha_plates=0.8
         faceAlpha_vlinks=0.4
         obstacle=[-20,-30,0; -20,30,0; -20,30,50; -20,-30,50]'; %the quadrilateral obstacle, without thickness.
     end

     methods
         function obj=Modules(Stiffnesses,Cable_displacements)
             obj.h=32;
             obj.k=obj.h/obj.r;
             Stiffnesses(Stiffnesses<0.1)=0.1;
             Stiffnesses(Stiffnesses>5)=5;
             Cable_displacements(Cable_displacements<0)=0;
             Cable_displacements(Cable_displacements>170)=170;
             obj.Stiffnesses=Stiffnesses;
             obj.Cable_displacements=Cable_displacements;
             obj.num=numel(obj.Stiffnesses)/(obj.N);
             obj.a_ini=obj.r*sqrt(2*(1-cos(obj.phi))+obj.k^2); 
             obj.b_ini=obj.r*sqrt(2*(1-cos(-obj.phi+2*pi/obj.N))+obj.k^2); 
             obj.l=2*obj.r*sin(pi/obj.N);         
             obj.as=obj.a_ini*ones(obj.N,1);
             obj.As=obj.a_ini*ones(obj.N*obj.num,1);
             obj.Bs=obj.b_ini*ones(obj.N*obj.num,1);
             obj.g_straight=[0;0;obj.phi;0;0;obj.k*obj.r];
             obj.base_vertices=obj.r_vec();
         end

         function Bs_result=forward(obj,n_iter)
             %main function, used to implement the forward solving. n_iter is default to 1.
             if nargin<2
                 n_iter=obj.n_iter;
             end
             obj.b_del=obj.Cable_displacements./n_iter;   %3X1 column vector

             Thetas_ini=zeros(size(obj.Stiffnesses));   %a column vector.
             options = optimset('Display','off','TolFun',1e-7,'TolX',1e-7,'Algorithm','Levenberg-Marquardt','MaxFunEvals',3000,'MaxIter',2000);
             Thetas=zeros(size(Thetas_ini));
             for i=1:obj.num
                   Thetas(((i-1)*3+1):3*i)=fsolve(@(thetas)obj.CirConstraints(obj.Bs(((i-1)*3+1):3*i),thetas,mod(i,2)),Thetas_ini(((i-1)*3+1):3*i),options);
                   %solving for the initial values of thetas according to the values of bs, when cables are not actuated.
             end
             Thetas_vector=Thetas;
             G=repmat(obj.g_straight,obj.num,1); %initial value for the tramsformation vector G for all the modules stacked in series.
             G_vector=G;
             Bs_vector=obj.Bs;

             Fval=obj.potential_energy_addsup(obj.Bs,Thetas); 
             Fval_vector=Fval; %store the values of potential energy for each iterations.
             %the first column of these variables (Thetas_vector, G_vector, Bs_vector, Fval_vector) denote the modules in their initial state.
             %So the numbers of columns of these variables are 1+n_iter
             obj.count=1;
             Proper=true;
             while (Proper==true) && (obj.count<=n_iter)
                 Displacements=obj.b_del*obj.count; %temporary value of cable displacements.    
                 [Bs_result,Fval,Proper,G_new,solution_Thetas]=obj.step_energy(G_vector(:,obj.count),Bs_vector(:,obj.count),Thetas_vector(:,obj.count),Displacements);
                 %The above line is the most important function in the simulation
                 if Proper
                    G_vector=[G_vector,G_new]; %#ok<*PROPLC,*AGROW> 
                    Bs_vector=[Bs_vector,Bs_result];
                    Thetas_vector=[Thetas_vector,solution_Thetas];
                    Fval_vector=[Fval_vector,Fval];                                         
                 end
                 obj.count=obj.count+1; 
             end

             obj.G_vector=G_vector;
             obj.Bs_vector=Bs_vector;
             obj.Thetas_vector=Thetas_vector;
             obj.Fval_vector=Fval_vector;
             if (obj.count-1)==n_iter
                obj.isSuccessful=true;
             end
         end



         function [Bs_result,Fval,Proper,G_new,solution_Thetas]=step_energy(obj,G,Bs,Thetas,Displacements)
             %search for the values of bs, for which the modules possess the minimum potential energy.
             %potential energy is calculated for a given Bs in the function "potential_energy_addsup".
             %The constraint is that, for every iteration, we reduce the cable lengths by b_del (b_del is a 3X1 vector).
             %Note that Thetas here is only used as a convinient initial guess.

             options=optimset('display','off','Algorithm','active-set','MaxFunEvals',10000);              
             A=repmat(eye(3),1,obj.num); %linear inequality constraint.
             b=obj.b_ini*obj.num-Displacements;

             [Bs_result,Fval]=fmincon(@(Bs)obj.potential_energy_addsup(Bs,Thetas),Bs,A,b,[],[],repmat(4,obj.N*obj.num,1),Bs,[],options);
             % search for the values of bs, for which the modules possess the minimum potential energy, while the constraint is satisfied.
             %If the cables are in tension, then A*bs=b. If any cable is slack, then A*bs<b
            
             [Proper,G_new,solution_Thetas]=obj.step(Bs_result,G,Thetas); 
             %according to the new value of bs, update some relavant geometric parameters.
         end 



          function pe_sum=potential_energy_addsup(obj,Bs,Thetas)
              %this function adds up the poetntial energies of all the stacked modules, through the function potential_energy
              %The input argument Thetas would only be used as the initial guess in the fsolve function, to solve for the geometry.
              pe=zeros(obj.num,1);
              if nargin<3
                 Thetas=zeros(size(obj.Stiffnesses)); 
              end              
              for i=1:obj.num
                 pe(i)=obj.potential_energy(Bs(((i-1)*3+1):3*i),obj.Stiffnesses(((i-1)*3+1):3*i),Thetas(((i-1)*3+1):3*i),mod(i,2));
              end
              pe_sum=sum(pe);
          end


          function pe=potential_energy(obj,bs,module_stiffnesses,thetas,direction)
              %note: this function applies to a single module
              % The basic idea is that the geometry of the module can be determined from the values
              % of bs, then based on the geometry of the module, the bending energies of various joints
              % can be determined.
              %bs------------3X1 column vector
              %module_stiffnesses---------------3X1 column vector
              %thetas------------3X1 column vector
              %direction-----------scalar, boolean. direction of the module. 1 (logical true) represents counterclockwise module. 
              % 0 (logical false) represents clockwise.
              %pe------------scalar
              %The input argument thetas would only be used as the initial guess in the fsolve function, to solve for the geometry.
              if nargin<4
                 thetas=[0,0,0]'; 
                 direction=1;
              end
              options = optimset('Display','off','TolFun',1e-7,'TolX',1e-7,'Algorithm','Levenberg-Marquardt','MaxFunEvals',30000,'MaxIter',20000);
              [solution_thetas]=fsolve(@(thetas)obj.CirConstraints(bs,thetas,direction),thetas,options);  
              %The ablove lines are calculating the thetas from bs

              position=top_vertices(obj,bs,solution_thetas,direction);  
              % The above line calculates geometry (coordinates of all the vertices on top plate) from the thetas.              

              %---------------------------------------------------------------------
              %The code below will calculate the bending energy according to the geometry of the module
              b=obj.b_ini/2; %scalar
              bending_angles=1-(bs.^2)/(2*b^2);  %3X1 vector
              bending_angles=acos(bending_angles);
              bending_angles=pi-bending_angles;

              vertical_link=position-obj.r_vec(); %both position and obj.r_vec are 3X3 matrices
              unit_vertical=vertical_link./vecnorm(vertical_link); 
              bot_vertical=dot(unit_vertical,repmat([0,0,1]',1,3));
              bot_vertical=acos(bot_vertical); %angles, in radian, between the vertical links and the bottom plate (1X3 vector)
              
              normal=Modules.get_normal(position);
              top_vertical=dot(unit_vertical,repmat(normal(:),1,3));
              top_vertical=acos(top_vertical);   %angles, in radian, between the vertical links and the top plate (1X3 vector)

              %--------------------
              U_diag=(1/2*module_stiffnesses)'*(bending_angles.^2);              
              stiffness30A=obj.stiffness_plateJoints;  
              U_vert=sum(1/2*stiffness30A*bot_vertical.^2)+sum(1/2*stiffness30A*top_vertical.^2);
              pe=U_diag+U_vert;
          end
    

          function diff_lengths=CirConstraints(obj,bs,thetas,direction)
             %note: this function applies to a single module
             %thetas------------3X1 column vector
             %bs--------------3X1 column vector
             %direction-----------scalar, boolean. direction of the module. 1 (logical true) represents counterclockwise module. 
             % 0 (logical false) represents clockwise.
             %diff_lengths----------------3X1 column vector         
             position=obj.top_vertices(bs,thetas,direction);
             position2=circshift(position,-1,2);
             diff_vertices=position2-position;
             lengths=sqrt(sum(diff_vertices.^2,1)); %a 1X3 row vector
             diff_lengths=lengths(:)-obj.l;   %a 3X1 column vector
         end

         function position=top_vertices(obj,bs,thetas,direction)
             %note: this function applies to a single module
             %bs--------------3X1 column vector
             %thetas------------3X1 column vector, in radian            
             %direction-----------scalar, boolean. direction of the module. 1 (logical true) represents counterclockwise module. 
             % 0 (logical false) represents clockwise.
             position=obj.OM_vector_func(bs,direction)+obj.MQ_vector_func(bs,thetas);
             if ~direction
                 position=circshift(position,1,2);
             end                
         end

         function MQ_vector=MQ_vector_func(obj,bs,thetas)
             %this function applies to a single module
             %used to obtain the unit vectors parallel with the line MQ. (see Fig. 4 in the paper 2022TMECH).
             %bs-----------3X1 column vector
             %thetas------- 3X1 column vector
             %MQ_vector----------3X3 matrix, columns represent different side, rows represent lengths of the unit_vector in xyz dimensions.
             ratio_PM=(1+(obj.as.^2-bs.^2)/obj.l^2)/2;  %a 3X1 column vector. it can be a negative value if angle QP1P2 is obtuse.
             ratio_MQ=sqrt(obj.as.^2/(obj.l^2)-ratio_PM.^2);     %a 3X1 column vector. also applies even when angle QP1P2 is obtuse.
             MQ_length=ratio_MQ.*obj.l;   %a 3X1 column vector        
             MQ_unit_vector=zeros(3,obj.N);
             base_vertices2=circshift(obj.base_vertices,-1,2); %shift by -1 (1 to the left) in the row dimension
             diff_base_vertices=base_vertices2-obj.base_vertices;
             diff_base_vertices=diff_base_vertices./vecnorm(diff_base_vertices); %make the rotation axes unit vectors
             for i=1:obj.N  %i indicates the numbering of face, in a counterclockwise sense.
                MQ_unit_vector(:,i)=Modules.R_axis(diff_base_vertices(:,i),thetas(i))*[0;0;1]; %positive value of thetas means bending outward
             end
             MQ_vector=repmat(MQ_length',3,1).*MQ_unit_vector;
         end

         function base_vertices=r_vec(obj,r)
             if nargin<2
                 r=obj.r;
             end
             %note: this function applies to a single module
             % used to obtain the coordinates of the three vertices on the fixed bottom plate, 3X3 matrix
             % I have changed the i to i-1, which is different from Ben's code,ensuring
             % that the first vertex/vector is along the X axis
             base_vertices=zeros(3,obj.N);  %3X3 matrix, column represents xyz dimension, row represents different vertices
             for i=1:obj.N
                base_vertices(:,i)=r.*(Modules.Rz((i-1)*2*pi/obj.N)*[1;0;0]);
             end
         end

         function OM_vector=OM_vector_func(obj,bs,direction)
             %note: this function applies to a single module
             % used to obtain the vectors OM. (see Fig. 4 in the paper 2022TMECH).
             %bs-----------3X1 column vector
             %direction-----------scalar, boolean. direction of the module. 1 (logical true) represents counterclockwise module. 
             % 0 (logical false) represents clockwise.
             %OM_vector----------3X3 matrix, columns represent different side, rows represent lengths of vector OMs in xyz dimensions.
             ratio_PM=(1+(obj.as.^2-bs.^2)/obj.l^2)/2;  %a 3X1 column vector
             base_vertices2=circshift(obj.base_vertices,-1,2); %shift by -1 (1 to the left) in the row dimension             
             if direction
                OM_vector=repmat(ratio_PM',3,1).*(base_vertices2-obj.base_vertices)+obj.base_vertices;
             else
                OM_vector=repmat(ratio_PM',3,1).*(-base_vertices2+obj.base_vertices)+base_vertices2;
             end             
         end

         function position=g2position(obj,g)
             %note: this function applies to a single module
             %Based on values of obj.base_vertices, solve for the position of the top plate via the value of g, instead of bs and thetas.
             %g-----------6X1 column vector
             %position------------3X3 matrix. each column represents a vertex on the top plate.
             rotation=Modules.R(g(1:3));
             p=g(4:6);
             position=rotation*obj.base_vertices+p;
         end


         function g=bs2g(obj,bs,thetas,g_ini,direction)
             %update value of g based on new bs. Make sure that thetas is solved based on the new bs, not just an initial guess.
             %This function works for one single module. 
             if nargin<4
                 g_ini=obj.g_straight;
             end 
             position=obj.top_vertices(bs,thetas,direction);
             options=optimset('display','off','Algorithm','levenberg-marquardt');
             g=fsolve(@(g)g2position_objFunc(obj,g,position),g_ini,options);             
             function out=g2position_objFunc(obj,g,position)
                 out=obj.g2position(g)-position;
                 %out is a 3X3 matrix
             end
         end

     end 
           


     methods (Static)
         function out = Rz(a)
             %input angle a is in radian
             out = [ cos(a), -sin(a), 0;...
                     sin(a),  cos(a), 0;...
                       0,       0,    1];
         end

         function out = Ry(a)
             %input angle a is in radian
             out = [cos(a), 0, sin(a);...
                       0,   1,    0;...
                    -sin(a), 0, cos(a)];
         end

         function out = Rx(a)
             out = [ 1,    0,     0;...
                     0, cos(a), -sin(a);...
                     0, sin(a), cos(a)];
         end

         function out = R(x)
             % chenzhe: Note the sequence of the three angles 
             out = Modules.Rx(x(1))*Modules.Ry(x(2))*Modules.Rz(x(3));
         end

         function out = R_axis(axis,angle)
             %given the axis x (x should be unit vector defining the rotation direction) and angle a, determine the corresponding rotation matrix
             axis = axis(:);
             out = cos(angle)*eye(3)+sin(angle)*skew(axis)+(1-cos(angle))*(axis*axis'); 
             %the above equation agrees well with the Rodrigue's ratation formula.
             %it's a different form of matrix representation compared to that matrix
             %representation shown in wikipedia.
             function out = skew(x)
                 out =  [0,-x(3),x(2);x(3),0,-x(1);-x(2),x(1),0];
             end
         end
         
         function T=g2T(g)
             %applies to a single module, used to convert the transformation vector to transformation matrix.
             T=[Modules.R(g(1:3)),g(4:6);0,0,0,1];
         end

         function new_position=transform(T,position)
             %applies to a single module, used to transform position by homogeneous transformation matrix T
             new_position=T(1:3,1:3)*position+T(1:3,4);
         end

         function newT=multi_T(T,multiplier)
             %T is a homogeneous transformation matrix T. 
             %multiplier can be either a transformation vector or a transformation matrix
             if size(multiplier,2)==1
                newT=T*Modules.g2T(multiplier);
             else
                newT=T*multiplier;
             end
         end

         function centroid=get_centroid(plane)
             %this function is used to calculate the centroid of a triangular(or quadrilateral) plane.
             %plane-----------------3X3 or 3X4 matrix. each column represents a vertex of the triangle (or quadrilateral).
             %centroid------------------3X1 column vector
             if size(plane,2)==3
                centroid=mean(plane,2); 
             elseif size(plane,2)==4
                centroid123=mean(plane(:,[1,2,3]),2); 
                centroid143=mean(plane(:,[1,4,3]),2);
                r=centroid143-centroid123;
                centroid214=mean(plane(:,[2,1,4]),2); 
                centroid234=mean(plane(:,[2,3,4]),2);
                s=centroid234-centroid214;
                t=norm(cross(centroid214-centroid123,s))/norm(cross(r,s));
                centroid=centroid123+t*r;
             end
         end

        function normal=get_normal(planar_plate)
            %This function is used to get the normal to a triangular planar plate
            %planar_plate--------------3X3 matrix, each column represents a vertex of the triangular plate.
            %normal----------------normal to the plate, 3X1 column vector, of unit length
            vector1=planar_plate(:,2)-planar_plate(:,1);
            vector2=planar_plate(:,3)-planar_plate(:,1);
            normal=cross(vector1,vector2);
            normal=normal/norm(normal);
        end

        function isInside=insideTriangle(triangle,point)
            %This function is used to check if the point is inside the triangle
            %triangle-------------3X3 matrix, each column represents a vectex.
            %point---------------3X1 column vector, the point be checked.
            areaTriangle=Modules.get_area(triangle);
            areaDAB=Modules.get_area([point(:),triangle(:,1:2)]);
            areaDBC=Modules.get_area([point(:),triangle(:,2:3)]);
            areaDCA=Modules.get_area([point(:),triangle(:,[3,1])]);
            if areaTriangle>=(areaDAB+areaDBC+areaDCA-1e-6)
               isInside=true;
            else
               isInside=false;
            end
        end

        function area=get_area(triangle)
           %this function is used to calculate the area of a triangle, with the method of cross product
           %triangle-------------3X3 matrix, each column represents a vectex.
           %area------------------scalar
           area=cross(triangle(:,2)-triangle(:,1),triangle(:,3)-triangle(:,1));
           area=1/2*norm(area);
        end

        isCollision=Ray2FaceCollision(plane,ray)
        isCollision=Face2FaceCollision(plane1,plane2)
     end

end 
    

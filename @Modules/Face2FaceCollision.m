function isCollision=Face2FaceCollision(plane1,plane2)
     %this static class method is used to detect if the plane1 collides with plane2 or not.
     %plane1------------------3X3 or 3X4 matrix. each column represents a vertex of the triangle (or quadrilateral).
     %plane2------------------3X3 or 3X4 matrix. each column represents a vertex of the triangle (or quadrilateral).
     %isCollision-----------------boolean. if true, then collision exists.  

    %test if the two planes are two far away from each other. if yes, then no collision
    cenPlane1=Modules.get_centroid(plane1);  %3X1 vector
    cenPlane2=Modules.get_centroid(plane2);  %3X1 vector
    distance=norm(cenPlane2-cenPlane1);
    Rmax_Plane1= max(vecnorm(plane1-cenPlane1));
    Rmax_Plane2= max(vecnorm(plane2-cenPlane2));
    if distance>Rmax_Plane1+Rmax_Plane2
        isCollision=false;
        return 
    end

    plane2_shift=circshift(plane2,-1,2);   
    for i=1:size(plane2,2)
        isCollision=Modules.Ray2FaceCollision(plane1,[plane2(:,i),plane2_shift(:,i)]);
        if isCollision
           return
        end
    end

    plane1_shift=circshift(plane1,-1,2);
    for i=1:size(plane1,2)
        isCollision=Modules.Ray2FaceCollision(plane2,[plane1(:,i),plane1_shift(:,i)]);
        if isCollision
           return
        end
    end
end
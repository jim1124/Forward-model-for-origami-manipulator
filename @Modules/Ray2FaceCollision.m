function isCollision=Ray2FaceCollision(plane,ray)
    %this static class method is used to detect if the ray (axis, link) collides with the obstacle or not.
    %plane----------------obstacle,3X3 or 3X4 matrix. each column represents a vertex of the triangle (or quadrilateral).
    %ray----------------3X2 matrix. each column represents an end point of the finite ray (or axis, link)
    %isCollision-----------------boolean. if true, then collision exists.
    
    %test if the finite ray and the plane are two far away from each other. if yes, then no collision
    cenPlane=Modules.get_centroid(plane);  %3X1 vector
    cenRay=mean(ray,2);    %3X1 vector
    distance=norm(cenRay-cenPlane);
    Rmax_Plane= max(vecnorm(plane-cenPlane));
    rmax_Ray=1/2*norm(ray(:,2)-ray(:,1));
    if distance>Rmax_Plane+rmax_Ray
        isCollision=false;
        return 
    end


    NormalPlane=Modules.get_normal(plane(:,1:3));


    %test if the starting point and ending point are on the same side of the plane
    if dot(NormalPlane,ray(:,1)-plane(:,1))*dot(NormalPlane,ray(:,2)-plane(:,1))>0
         isCollision=false;
         return
    end

    %now we deal with the situation where the ray intersect the infinite plane.
    V=ray(:,2)-ray(:,1);
    V=V/norm(V);   %unit direction vector of the ray
    if dot(V,NormalPlane)==0
       P=ray;
    else
       t=-(dot(ray(:,1),NormalPlane)-dot(plane(:,1),NormalPlane))/dot(V,NormalPlane);
       P=ray(:,1)+t.*V;  %P is the intersection point of the ray and the plane.
    end

    for i=1:size(P,2)
       if size(plane,2)==3
           isInside=Modules.insideTriangle(plane,P(:,i));
       elseif size(plane,2)==4
           isInside1=Modules.insideTriangle(plane(:,1:3),P(:,i));
           isInside2=Modules.insideTriangle(plane(:,[1,4,3]),P(:,i));
           isInside=isInside1 || isInside2;
       end
       if isInside
           isCollision=true;
           return
       end
    end

    isCollision=false;
end
        







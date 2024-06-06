function generate_animation(obj)
   % The code below would be used to plot the modules and generate an animation, representing the
   % deformation process of the deformed manipulator given a single input (stiffnesses and cable displacements).

   %THe deforming process starts from initial straight shape, when stiffnesses and cable
   %displacements are both set to zeros.
   
   n=size(obj.Bs_vector,2);
   figure
   for i=1:n
     obj.modulesPlot(obj.G_vector(:,i)); 
     pause(1);
     frame=getframe(gcf);
     imind=frame2im(frame);
     [imind,cm] = rgb2ind(imind,256);
     if i==1
          imwrite(imind,cm,'test.gif','gif', 'Loopcount',1,'DelayTime',0.2);
     else
          imwrite(imind,cm,'test.gif','gif','WriteMode','append','DelayTime',0.2);
     end
   end
end
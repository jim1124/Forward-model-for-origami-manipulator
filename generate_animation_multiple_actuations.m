function generate_animation_multiple_actuations(actuation_inputs)
   % 2024.6.5 The code below would be used to plot the modules/manipulator corresponding to each input (stiffnesses and cable displacements) 
   % and generate an animation, representing the morphing history of a manipualator given a series of input (stiffnesses and cable displacements).

   %This code is different from the code "generate_animation.m", in the sense that
   %"generate_animation.m" takes a single input (stiffnesses and cable displacements) and generate
   % an animation of the deforming process of the manipulator starting from initial staright shape.

   % This code does not consider the deforming process for a single given input. 
   % Only the final shape of each input is considered and plotted.

   %actuation_inputs is the series of stiffnesses and displacements.
   %

   n=size(actuation_inputs,2);
   figure
   for i=1:n
     Cable_displacements=actuation_inputs(1:3,i);
     Stiffnesses=actuation_inputs(4:15,i);
     M=Modules(Stiffnesses,Cable_displacements);
     M.forward();
     M.modulesPlot(M.G_vector(:,end)); 
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
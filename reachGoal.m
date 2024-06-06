
Cable_displacements=[97.9748,63.8870,11.1119]';
Stiffnesses=[0.516, 0.817, 1.189, 0.524, 1.205, 1.422, 0.443, 0.234, 0.593, 0.637, 0.848, 1.015]';
Input0=[Stiffnesses;Cable_displacements];
lb=[0.1*ones(size(Stiffnesses));zeros(size(Cable_displacements))];
ub=[5*ones(size(Stiffnesses));150*ones(size(Cable_displacements))];
% options=optimset('display','iter','Algorithm','interior-point','MaxFunEvals',10000); 

options=optimset('display','iter'); 
[Input,fval]=ga(@reach_pose,15,[],[],[],[],lb,ub,[],options);
% [Input,fval]=fminsearch(@reach_pose,Input0,options);



function errPose=reach_pose(Input)
      Input=Input(:);  % a column vector
      Stiffnesses=Input(1:end-3);
      Cable_displacements=Input(end-2:end);
      M=Modules(Stiffnesses,Cable_displacements);
      M.forward();
      errPose=M.postAnalyze();
end
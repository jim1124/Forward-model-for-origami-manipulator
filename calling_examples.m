tic
Cable_displacements=[38.8834,54.3475,17.6192]';
Stiffnesses=[1.8268, 1.2647,0.1951, 0.5570, 1.0938, 1.9150, 1.0298, 0.3152, 1.9412, 1.9143, 0.9708, 1.6006]';
M=Modules(Stiffnesses,Cable_displacements);
B=M.forward();
toc


% tic
% Cable_displacements=[43.39706;43.391045;42.097256];
% Stiffnesses=[32.677837;49.99903;0.21039307;0.10000446;2.658546;49.99465;49.99979;49.753967;49.999752;49.99901;0.10267834;38.863655]';
% M=Modules(Stiffnesses,Cable_displacements);
% B=M.forward();
% toc


% Cable_displacements=rand(3,1)*100;
% Stiffnesses=rand(12,1)*2;
% tic
% Cable_displacements=[97.9748,63.8870,11.1119]';
% Stiffnesses=[0.516, 0.817, 1.189, 0.524, 1.205, 1.422, 0.443, 0.234, 0.593, 0.637, 0.848, 1.015]';
% M=Modules(Stiffnesses,Cable_displacements);
% B=M.forward();
% toc
% M.modulesPlot();
% M.generate_animation()

% Cable_displacements=[97.9748,63.8870,81.1119]';
% Stiffnesses=[2, 2, 2, 0.2, 0.2, 2, 0.2, 0.2, 2, 2, 2, 2]';
% M=Modules(Stiffnesses,Cable_displacements);
% % 不必要每次摧毁、重建 class instance
% B=M.forward(); %12 X 1 column vector
% [collision,errPose]=M.postAnalyze(); 
% M.modulesPlot();



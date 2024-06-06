%Now that we have the forward model (using the function recursive_flow), we will explore the inverse process.
%  that is, given a shape of the arm (bs of all diagonal links), how to determine the input variables of all diagonal links 
%  and the cable lengths on three sides.

%We are using the "Coordinate Descent" algorithm here to update the variable stiffnesses (called S in this script).

% stiffness=[0.34, 0.44, 0.35, 0.22, 0.68, 0.35, 0.35, 0.42, 0.24]';
%2023.9.5 [g_out,jiaodu,bs_vector,fval2_vector]=recursive_flow(17.8979,32/17.8979,10/180*pi,stiffness,[20,30,10]',10);
% The above two lines show the stiffness corresponding to the following T.
% T=18+22*rand(9,1); %random numbers in the interval (18,40)
%T=[38.9741, 30.9287, 42.1302, 34.9613, 42.7524, 43.4027, 36.2725, 26.5269, 34.6750]';
tic
T=[ 43.4027
   43.4027
   43.4027
   5
   5
   43.4027
   15
   15
   15
   5
   43.4027
   43.4027];

% T is a 12X1 column vector, representing the target values of bs.

S=[0.3, 0.4, 0.5, 0.2, 0.38, 0.25, 0.35, 0.42, 0.24,1,1.4,0.28]'; %initial value of the variables, 9X1 column matrix

b_max= 43.4027; %depending on the geometry of the module
num=numel(T)/3; %number of modules

cable_lengths=sum(reshape(T,3,[]),2); %column vector 3X1
final_displacement=b_max*num-cable_lengths;

n_steps=1000;
step=0;
changed=false;


S_sequence=S;
ds=0.1; %note dS is a scalar, representing the distance of new_direction, which is the difference between new_S and S.
dS=zeros(numel(S),1);

Y=recursive_flow(S,final_displacement);
current_error=rmse(T,Y);
disp(current_error);
error_sequence=current_error; 
ds_sequence=[]; 


while step<n_steps
   
   step=step+1;
   disp(['step ', num2str(step)]);
   if ~changed
       ds=0.3*ds;
   end

   changed=false;
   
   for i=1: numel(S)
        dS=zeros(numel(S),1);
        dS(i,1)=ds; %Note dS is a column vector, with just one element being nonzero

        up_Y=recursive_flow(S+dS,final_displacement);
        up_error=rmse(T,up_Y);

        down_Y=recursive_flow(S-dS,final_displacement);
        down_error=rmse(T,down_Y);

        if down_error<current_error
            dS=-dS;
        end


        while 1
            new_S=S+dS;
            new_Y=recursive_flow(new_S,final_displacement);
            
            new_error=rmse(T,new_Y);
            if (new_error>=current_error) || (step>n_steps)
                break
            end
            changed=true;
            S=new_S;
            S_sequence=[S_sequence,S];
            current_error=new_error;
            disp(current_error);
            error_sequence=[error_sequence,current_error];
            ds_sequence=[ds_sequence,ds];
            step=step+1;
            disp(['step ', num2str(step)]);
        end

   end

   if new_error<1e-4
       break
   end
   toc
end

function error=rmse(T,Y)
     error=mean((T-Y).^2);
     error=sqrt(error);
end

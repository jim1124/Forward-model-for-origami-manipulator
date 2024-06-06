function [Proper,G_new,solution_Thetas]=step(obj,Bs_result,G,Thetas)
    % according to the new Bs, update some relevant parameters.   
    %Note that this function applies on all modules as a whole, so Bs_result is num*N X 1 vector
    %Note that Thetas here is only used as a convinient initial guess.
    %vector G here is also ised as a convinient initial guess.
    
    Bs=Bs_result; %for all modules, not a single module
    solution_Thetas=zeros(size(Bs));
    Fval=zeros(size(Bs));
    G_new=zeros(6*obj.num,1); % for every module, its g is a six-element column vector.
    options = optimset('Display','off','TolFun',1e-7,'TolX',1e-7,'Algorithm','Levenberg-Marquardt','MaxFunEvals',3000,'MaxIter',2000);
    for ii=1:obj.num
      [solution_Thetas(((ii-1)*3+1):3*ii),Fval(((ii-1)*3+1):3*ii)]=...
                     fsolve(@(thetas)obj.CirConstraints(Bs(((ii-1)*3+1):3*ii),thetas,mod(ii,2)),Thetas(((ii-1)*3+1):3*ii),options); 
    end    
    constraints=all(abs(Fval)<1e-6,'all'); %to verify that the above fsolve optimization makes sense.
    %used to display the value of CirConstraints
    extension=all(Bs<=obj.b_ini,'all');
    Proper=constraints && extension;
    if (~Proper) && obj.display_error       
       disp('erroneous_Bs=');
       disp(Bs);
    end

    for i=1:obj.num
       [G_new(((i-1)*6+1):6*i)]=bs2g(obj,Bs(((i-1)*3+1):3*i),solution_Thetas(((i-1)*3+1):3*i),G(((i-1)*6+1):6*i),mod(i,2));  
       %Note that g for every module is independent, solely dependent on its corresponding bs
    end
end
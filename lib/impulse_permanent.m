%% 
% @author   Roland Meeks, University of Essex
% @date     2012-10
% @see      impulse_zeros.m
% @purpose  To produce impulse-response functions in which the shocked 
%           variable undergoes a permanent change.
function response = impulse_permanent( By, smat, nstep, shock_index )

[neq,nvar,nlag] = size(By);


% init
response = zeros(nvar,neq,nstep);
response(:,shock_index,1) = smat; % impact response


% iterate up to simulation horizon
for it=2:nstep,
    
    for ilag=1:min(nlag,it-1)

        response(:,:,it) = ...
            response(:,:,it)+By(:,:,ilag)*response(:,:,it-ilag);
        
    end;
    
    response(shock_index,:,it) = response(shock_index,:,it-1);
    
end;

%% end of function
end
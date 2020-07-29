%% 
% @author   Roland Meeks, University of Essex
% @date     2012-10
% @see      impulse.m
% @purpose  To produce impulse-response functions in which one or more
%   variables are 'zeroed out', i.e. are constrained not to respond to the
%   given shock. This 'closes down' some of the transmission channels. For
%   example, if we were to look at output following an oil shock with the
%   monetary policy response shut down, would there be a larger or smaller
%   recession (Bernanke, Gertler, Watson)? The change in the 'policy rule',
%   to assume no response, leaves all other equations unchanged, a
%   violation of the standard cross-equation restrictions. So the preferred
%   interpretation is that there are a sequence of surprises which entail
%   the variable specified remains at exactly zero through the simulation
%   horizon.
function response = impulse_zeros( By, smat, nstep, var_index )

[neq,nvar,nlag] = size(By);


% init
response = zeros(nvar,neq,nstep);
response(:,:,1) = smat'; % impact response
response(var_index,:,1) = 0;


% iterate up to simulation horizon
for it=2:nstep,
    
    for ilag=1:min(nlag,it-1)

        response(:,:,it) = ...
            response(:,:,it)+By(:,:,ilag)*response(:,:,it-ilag);
        
    end;
    
    response(var_index,:,it) = 0;
    
end;

%% end of function
end
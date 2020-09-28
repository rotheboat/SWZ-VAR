function isidentified = CheckGlobalIdentification( mQ, mR )
% Checks whether the sufficient conditions for global identification of the
% structural VAR are satisfied using Rubio-Ramirez, Waggoner, and Zha
% (ReStud, 2010; Thm. 1) condition.
%
% For references to notation see ../tex/LinearRestrictions.pdf
%
% @args: Matrix/matrices containing linear restrictions on A or F
%   <array> mQ = a M x M x M array such that mQ(:,:,j)*A(:,j) = 0
%   <array> mR = a K x K x M array such that mR(:,:,j)*F(:,j) = 0
% @returns: <bool> true if sufficient condition for global identification
%   is met.
% @calledby: SetSimsZhaPrior
% @author: Roland Meeks
% @date: August, 2020
cM = size( mQ, 1 );

for idx = 1:cM
    % Step 1: Construct the matrix of 'stacked' restrictions L
    mL(:,:,idx) = blkdiag( mQ(:,:,idx), mR(:,:,idx) );
    % Step 2: Count the number of restrictions in each equation ('ell')
    ell(idx) = rank( mL(:,:,idx) );
    % Step 3: Construct random A and F matrices that satisfy the
    % restrictions; these matrices are used to check Thm 1. The theorem
    % holds a.s., i.e. the set of parameters for which its conditions fail
    % to hold are measure zero (but that may include some obvious choices,
    % such as all ones).
    %
    % We find the null space of the restriction matrix mQ for each
    % equation, then draw a random vector of the free parameters for that
    % equation.
    mU = null( mQ(:,:,idx) );
    % Note that the number of free parameters is equal to the rank of mU
    mA(:,idx) = mU*rand( rank( mU ), 1 );
    % Same procedure for mR
    mV = null( mR(:,:,idx) );
    mF(:,idx) = mV*rand( rank( mV ), 1 );
end

% Step 4: Reorder equations so the the number of restrictions is decreasing
[~, index] = sort( ell, 'descend' );
mL = mL(:,:,index);
% The function f( A, F ), in the correct order
eff = cat( 1, mA(:,index), mF(:,index) );

for idx = 1:cM
    % Step 5: For each equation, construct the M matrix required for theorem 1.
    mM{idx} = ...
        cat( 1, mL(:,:,idx)*eff, cat( 2, eye(idx), zeros(idx,cM-idx) ) );
    % Step 6: For each M matrix, check rank(M) = m (the dimension of the VAR)
    emm(idx) = rank( mM{idx} );
end

% Step 7: Return flag isidentified
isidentified = all( eq( emm, cM ) );

end

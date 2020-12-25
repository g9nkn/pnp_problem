%% Globally optimal solution to the PnP problem by G. Nakano [1].
% Copyright (c) 2020 NEC Corporation
% This software is released under the NEC Corporation License. See License.txt.
% For commercial use, please contact Gaku Nakano <g-nakano@nec.com>.
%
% USAGE:
%   [R, t, err] = pnp_nakano_bmvc2015(pts3d, pts2d, rot_params)
%
% INPUTS:
%   pts3d - 3xN matrix of 3D points corresponding to the 2D points. (N >= 3)
%           [x1, x2, ..., xN
%            y1, y2, ..., yN
%            z1, z2, ..., zN];
%   pts2d - 2xN matrix of 2D points corresponding to the 3D points. (N >= 3)
%           Each points are normalized by using the intrinsic paramters.
%           [u1, u2, ..., uN
%            v1, v2, ..., vN];
%   rot_params - char to select the paramterization of the PnP problem.
%            'rotation_matrix'   : 3x3 Rotation matrix
%            'quaternion'        : 4x1 Quaternion
%            'cayley' or 'optdls': 3x1 Cayley transformation
%             All solvers were geenrated by using [2].
%
% OUTPUS:
%   R - 3x3xM rotation matrix. R(:,:,i) corresponds to t(:,i) and err(i).
%   t - 3xM translation vector.
%   err - 1xM algebraic cost of the solution.
%
% REFERENCE:
%   [1] Gaku Nakano, "Globally Optimal DLS Method for PnP Problem with
%       Cayley parameterization," BMVC2015.
%   [2] V. Larsson et al., "Efficient Solvers for Minimal Problems by
%       Syzygy-based Reduction," CVPR 2017.
%       http://people.inf.ethz.ch/vlarsson/misc/autogen_v0_5.zip
function [R, t, err] = pnp_nakano_bmvc2015(pts3d, pts2d, rot_params)

    
    %% check arguments
    narginchk(3,3);
    assert( size(pts3d,1)==3 )
    assert( size(pts2d,1)==2 )
    assert( size(pts3d,2)==size(pts2d,2) && size(pts3d,2) >= 3 )

    %% select solver corresponding to the rotational representation
    switch lower(rot_params)
        case 'rotation_matrix'
            solver_func = @solver_pnp_rotmat;
            R_rand      = eye(3);
        case 'quaternion'
            solver_func = @solver_pnp_quat;
            R_rand      = eye(3);
        case {'cayley', 'optdls'}
            solver_func = @solver_pnp_cayley;
            q_rand      = rand(4,1);
            q_rand      = q_rand / norm(q_rand);
            R_rand      = quat2rot( q_rand );
            pts3d       = R_rand * pts3d;
        otherwise
            error('rot_params is wrong. Check valid options by ''help pnp_nakano_bmvc2015''');
    end
    
    %% normalize 3D points
    [pts3d_n, mean3d, var3d] = normalize3dpts(pts3d);

    
    %% call Groebner basis solver
    [M, T] = calcMT(pts3d_n, pts2d);
    R      = solver_func(M);
    
    
    %% recover t and de-normalization
    t   = zeros(3,size(R,3));
    err = zeros(1,size(R,3));
    for i=1:size(R,3)
        R_tmp    = R(:,:,i)';
        r_tmp    = reshape(R_tmp',9,1);
        t_tmp    = T * r_tmp;
        t(:,i)   = t_tmp*var3d - R_tmp*mean3d;
        R(:,:,i) = R_tmp * R_rand;
        err(i)   = r_tmp' * M * r_tmp;
    end
    
        
end


npts  = 100;
planar_scene = false;


% image points
w = 640;
h = 480;
f = 500;
K = [f, 0, w/2
     0, f, h/2
     0, 0,  1];
pts2d = [w*rand(1, npts)
         h*rand(1, npts)
         ones(1, npts)];

% projective depth
min_depth = 5;
max_depth = 10;
if planar_scene
    % 3D plane: normal'*X+p = 0 where X=d*K\m
    theta  = (2*pi/4)*rand(1) - pi/4;
    phi    = pi/3*rand(1);
    normal = -[sin(theta)*cos(phi)
               sin(theta)*sin(phi) 
               cos(phi)];
    p = (max_depth - min_depth)*rand(1);
    d = - p ./ (normal'*(K\pts2d));
else
    d = min_depth + (max_depth - min_depth)*rand(1, npts);
end


% 3D points d*m = K*(R*X+t) -> X = R'*(d*K\m - t)
q = rand(4,1);
R = quat2rot( q/norm(q) );
t = min_depth/2 * rand(3,1);
pts3d = R' * (d.*(K\pts2d) - t);


% add image noise 
sigma = 2;
noisy_pts2d = [pts2d(1:2,:) + sigma*randn(2,npts)
               ones(1,npts)];
noisy_pts2d = K\noisy_pts2d;


% solve pnp problem
[R_r, t_r, algerr_r] = pnp_nakano_bmvc2015(pts3d, noisy_pts2d(1:2,:), 'rotation_matrix');
[R_q, t_q, algerr_q] = pnp_nakano_bmvc2015(pts3d, noisy_pts2d(1:2,:), 'quaternion');
[R_c, t_c, algerr_c] = pnp_nakano_bmvc2015(pts3d, noisy_pts2d(1:2,:), 'cayley');


disp('Estimation errors of rotation matrix are:')
disp(['    9x9 rotmat    : ' num2str(min(calc_R_err(R, R_r)))])
disp(['    quaternion    : ' num2str(min(calc_R_err(R, R_q)))])
disp(['    cayley(optDLS): ' num2str(min(calc_R_err(R, R_c)))])


function R_err = calc_R_err(R_gt, R_est)
    for i=1:size(R_est,3)
        R_err(i) = norm(R_gt'*R_est(:,:,i) - eye(3), 'fro');
    end
end
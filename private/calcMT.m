function [M, T] = calcMT(pts3d, pt2d)
	
	n = size(pts3d,2);
	
	
	%calc M2
	u = pt2d(1,:)'; 
	v = pt2d(2,:)';
	u2 = u.^2;
	v2 = v.^2;

	% M'*M
	sum_u = sum(u);
	sum_v = sum(v);
	MM = [     n,      0,    -sum_u
	           0,      n,    -sum_v
	      -sum_u, -sum_v, sum(u2+v2)];
	invMM = inv(MM);
	
	% N'*N
	X  = pts3d';
	uX = repmat(u,1,3).*X;
	vX = repmat(v,1,3).*X;
	N  = [        -X, zeros(n,3), uX
	      zeros(n,3),         -X, vX];
	NN = N'*N;

	% M'*N
	sum_X   = sum(X);
	sum_uX  = sum(uX);
	sum_vX  = sum(vX);
	sum_u2v2X = sum( (u2+v2).*X );
	
	MN = [     sum_X, zeros(1,3),    -sum_uX
	      zeros(1,3),      sum_X,    -sum_vX
	         -sum_uX,    -sum_vX,  sum_u2v2X];

	M = NN - MN'*invMM*MN;
	
	T = -invMM*MN;
	
return
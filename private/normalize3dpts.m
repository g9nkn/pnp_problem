% shift 3D data so that variance = sqrt(2), mean = 0
function [Xnew, mean3d, var3d] = normalize3dpts( X )
	
	n = size(X,3);
	
    mean3d = mean(X, 2);
    Xnew   = X - mean3d;    
    var3d  = sum( sqrt(sum( Xnew.^2 ) ) )/n;
    Xnew   = 1/var3d*Xnew;
	
return
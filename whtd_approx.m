        function [xhat,errhat,topt,kest] = whtd_approx(y,sig,...
            awhts,bwhts,m,n,k)
%
%                              description:
%
%   Estimates optimal spectral denoiser for data of the form
%
%                               Y = X + G,                          (1)
%
%   where X is a low-rank signal matrix, and G is white noise of 
%   variance sig^2. The method attempts to find the optimal k-by-k topt
%   to minimize the AMSE
%
%                || A (X - \hat U topt \hat V^T) B^T ||_{Fr}^2      (2)
%
%   where A and B are the diagonal weight matrices, and \hat U and \hat V
%   are the empirical singular vectors of Y.
%
%                               input parameters:
%
%   y - the m-by-n matrix of data vectors
%   sig - the noise level
%   awhts - m-dimensional vector of weights
%   bwhts - n-dimensional vector of weights
%   m,n - the dimensionalities
%   k - upper bound on the rank of X
%
%                              output parameters:
%
%   xhat - the estimated matrix
%   topt - estimated k-by-k matrix used in spectral denoiser
%   errhat - an estimate of the Frobenius loss between X and \hat X
%   kest - the rank of xhat
%   
%
        xhat=zeros(m,n);
        topt=zeros(k,k);
        errhat=0;

        gam=m/n;
        tol = 100*eps;

%
%        . . . SVD of data
%
        [uy,sy,vy] = whtd_svdsmart(y,m,n,k);
%
%        rank estimation
%
        bmargin = 0;
%%%        bmargin = 10/n^(2/3);

        kest = whtd_estrank_dumb(sy,gam,sig,k,bmargin);
%%%        prinf('kest=',kest,1);


        if (kest == 0)
%
        return
        end
%
        sy = sy(1:kest);
        uy=uy(:,1:kest);
        vy=vy(:,1:kest);

        [topt,errhat] = whtd_svdtopt(uy,sy,vy,sig,...
            awhts,bwhts,m,n,kest,tol);
%
%        the estimated matrix
%
        xhat = uy * topt * vy';

        end
%

%-------------------------------------------------------------------------------
% Nematic Bingham closure in three dimensions
%
% Input:
%   Q - second moment tensor, struct with fields q11,q12,q13,q22,q23,q33
%   T - rotation tensor T = E + 2*zeta*c*Q, struct with fields t11,t12,t13,t22,t23,t33
%   c11 - Chebyshev coefficients (nu1,nu2) --> ts1111
%   c12 - Chebyshev coefficients (nu1,nu2) --> ts1122
%   c22 - Chebyshev coefficients (nu1,nu2) --> ts2222
%
% Output:
%   ST - contraction S:T, struct with fields st11,st12,st13,st22,st23,st33
%
% Scott Weady, CIMS
% Last updated February 2022
%-------------------------------------------------------------------------------
function ST = bingham3d(Q,T,c11,c12,c22)
      
  % Get components of second moment tensor
  q11 = Q.q11; q12 = Q.q12+1e-15; q13 = Q.q13+1e-15;
  q22 = Q.q22; q23 = Q.q23+1e-15;
  q33 = Q.q33;

  % Get components of rotation tensor
  t11 = T.t11; t12 = T.t12; t13 = T.t13;
  t22 = T.t22; t23 = T.t23;
  t33 = T.t33;
  
  % Transformation between feasible domain and Chebyshev grid
  H1 = @(u,v) (u - 1/3) - (v - 1/3);
  H2 = @(u,v) 2*(u - 1/3) + 4*(v - 1/3);
  G1 = @(u,v) 2*(u + v) - 1;
  G2 = @(u,v) 2*(u./(u+v)) - 1;
  
  % Solve for largest eigenvalue
  tol = 1e-15;  
  mu = (1/3)*ones(size(q11));

  a = q11.*q22 + q11.*q33 + q22.*q33 - (q12.^2 + q13.^2 + q23.^2);
  b = q11.*q23.^2 + q22.*q13.^2 + q33.*q12.^2 - q11.*q22.*q33 - 2*q13.*q12.*q23;

  fnu = mu.^3 - mu.^2 + a.*mu + b;
  maxiter = 50;
  iter = 0;
  
  while max(abs(fnu(:))) > tol && iter < maxiter
    fnu = mu.^3 - mu.^2 + a.*mu + b;
    df = fnu./(3*mu.^2 - 2*mu + a);
    mup1 = mu - df;
    mu = mup1;
    iter = iter + 1;
  end

  nu1 = mu;
  nu2 = (-(mu-1) + sqrt((mu-1).^2 - 4*(a + mu.*(mu-1))))/2;
  nu3 = (-(mu-1) - sqrt((mu-1).^2 - 4*(a + mu.*(mu-1))))/2;

  mu1 = max(nu1,max(nu2,nu3));
  mu3 = min(nu1,min(nu2,nu3));
  mu2 = 1 - (mu1 + mu3);

  % First eigenvector
  O11 = q12.*q23 - q13.*(q22 - mu1);
  O21 = q13.*q12 - (q11-mu1).*q23;
  O31 = (q11-mu1).*(q22-mu1) - q12.*q12;
  m1 = sqrt(O11.^2 + O21.^2 + O31.^2);
  O11 = O11./m1; O21 = O21./m1; O31 = O31./m1;

  % Second eigenvector
  O12 = q12.*q23 - q13.*(q22-mu2);
  O22 = q13.*q12 - (q11-mu2).*q23;
  O32 = (q11-mu2).*(q22-mu2) - q12.*q12;
  m2 = sqrt(O12.^2 + O22.^2 + O32.^2);
  O12 = O12./m2; O22 = O22./m2; O32 = O32./m2;

  % Third eigenvector
  O13 = O21.*O32 - O31.*O22;
  O23 = O31.*O12 - O11.*O32;
  O33 = O11.*O22 - O21.*O12;
  m3 = sqrt(O13.^2 + O23.^2 + O33.^2);
  O13 = O13./m3; O23 = O23./m3; O33 = O33./m3;

  % Improve orthogonality of second eigenvector
  O12 = O21.*O33 - O31.*O23;
  O22 = O31.*O13 - O11.*O33;
  O32 = O11.*O23 - O21.*O13;
  m2 = sqrt(O12.^2 + O22.^2 + O32.^2);
  O12 = O12./m2; O22 = O22./m2; O32 = O32./m2;

  % Map to Chebyshev grid
  nu1 = G1(H1(mu1,mu2),H2(mu1,mu2));
  nu2 = G2(H1(mu1,mu2),H2(mu1,mu2));
        
  % Evaluate Chebyshev interpolants
  Tim2 = 1; Tjm2 = 1;
  Tim1 = nu1; Tjm1 = nu2;

  ts1111 = Tim2.*(c11(1,1).*Tjm2 + c11(1,2).*Tjm1) + ...
           Tim1.*(c11(2,1).*Tjm2 + c11(2,2).*Tjm1);
  ts1122 = Tim2.*(c12(1,1).*Tjm2 + c12(1,2).*Tjm1) + ...
           Tim1.*(c12(2,1).*Tjm2 + c12(2,2).*Tjm1);
  ts2222 = Tim2.*(c22(1,1).*Tjm2 + c22(1,2).*Tjm1) + ...
           Tim1.*(c22(2,1).*Tjm2 + c22(2,2).*Tjm1);
         
  M = length(c11);
  for i = 3:M
    Ti = 2*nu1.*Tim1 - Tim2;
    ts1111 = ts1111 + Ti.*(c11(i,1).*Tjm2 + c11(i,2).*Tjm1);
    ts1122 = ts1122 + Ti.*(c12(i,1).*Tjm2 + c12(i,2).*Tjm1);
    ts2222 = ts2222 + Ti.*(c22(i,1).*Tjm2 + c22(i,2).*Tjm1);
    Tim2 = Tim1;
    Tim1 = Ti;
  end

  for j = 3:M
    Tj = 2*nu2.*Tjm1 - Tjm2;
    Tim2 = 1;
    Tim1 = nu1;
    ts1111 = ts1111 + Tim2.*c11(1,j).*Tj + Tim1.*c11(2,j).*Tj;
    ts1122 = ts1122 + Tim2.*c12(1,j).*Tj + Tim1.*c12(2,j).*Tj;
    ts2222 = ts2222 + Tim2.*c22(1,j).*Tj + Tim1.*c22(2,j).*Tj;
%     for i = 3:(M-(j-1))
    for i = 3:M
      Ti = 2*nu1.*Tim1 - Tim2;
      ts1111 = ts1111 + Ti.*c11(i,j).*Tj;
      ts1122 = ts1122 + Ti.*c12(i,j).*Tj;
      ts2222 = ts2222 + Ti.*c22(i,j).*Tj;
      Tim2 = Tim1;
      Tim1 = Ti;
    end
    Tjm2 = Tjm1;
    Tjm1 = Tj;
  end
  
  ts1133 = mu1 - ts1111 - ts1122;
  ts2233 = mu2 - ts1122 - ts2222;
  ts3333 = mu3 - ts1133 - ts2233;

  disp(abs([ts1111(1) ts1122(1) ts2222(1)]-[1/5 1/15 1/5]))

  % Rotate T to diagonal coordinate system
  tt11 = O11.*(t11.*O11 + 2*t12.*O21 + 2*t13.*O31) + ...
         O21.*(t22.*O21 + 2*t23.*O31) + ...
         t33.*O31.*O31;   

  tt12 = t11.*O11.*O12 + t12.*(O11.*O22 + O21.*O12) + ...
         t22.*O21.*O22 + t13.*(O11.*O32 + O31.*O12) + ...
         t33.*O31.*O32 + t23.*(O21.*O32 + O31.*O22);

  tt13 = t11.*O11.*O13 + t12.*(O11.*O23 + O21.*O13) + ...
         t13.*(O11.*O33 + O31.*O13) + t22.*O21.*O23 + ...
         t23.*(O21.*O33 + O31.*O23) + t33.*O31.*O33;

  tt22 = O12.*(t11.*O12 + 2*t12.*O22 + 2*t13.*O32) + ...
         O22.*(t22.*O22 + 2*t23.*O32) + ...
         t33.*O32.*O32;  

  tt23 = t11.*O12.*O13 + t12.*(O12.*O23 + O22.*O13) + ...
         t13.*(O12.*O33 + O32.*O13) + t22.*O22.*O23 + ...
         t23.*(O22.*O33 + O32.*O23) + t33.*O32.*O33;

  tt33 = O13.*(t11.*O13 + 2*t12.*O23 + 2*t13.*O33) + ...
         O23.*(t22.*O23 + 2*t23.*O33) + ...
         t33.*O33.*O33;  

  % Compute contraction tS:tT
  tst11 = ts1111.*tt11 + ts1122.*tt22 + ts1133.*tt33;
  tst12 = 2*ts1122.*tt12;
  tst13 = 2*ts1133.*tt13;
  tst22 = ts1122.*tt11 + ts2222.*tt22 + ts2233.*tt33;
  tst23 = 2*ts2233.*tt23;
  tst33 = ts1133.*tt11 + ts2233.*tt22 + ts3333.*tt33;

  % Rotate tS:tT to original coordinate system
  st11 = O11.*(O11.*tst11 + O12.*tst12 + O13.*tst13) + ...
                   O12.*(O11.*tst12 + O12.*tst22 + O13.*tst23) + ...
                   O13.*(O11.*tst13 + O12.*tst23 + O13.*tst33);
  st12 = O11.*(O21.*tst11 + O22.*tst12 + O23.*tst13) + ...
                   O12.*(O21.*tst12 + O22.*tst22 + O23.*tst23) + ...
                   O13.*(O21.*tst13 + O22.*tst23 + O23.*tst33);
  st13 = O11.*(O31.*tst11 + O32.*tst12 + O33.*tst13) + ...
                   O12.*(O31.*tst12 + O32.*tst22 + O33.*tst23) + ...
                   O13.*(O31.*tst13 + O32.*tst23 + O33.*tst33);
  st22 = O21.*(O21.*tst11 + O22.*tst12 + O23.*tst13) + ...
                   O22.*(O21.*tst12 + O22.*tst22 + O23.*tst23) + ...
                   O23.*(O21.*tst13 + O22.*tst23 + O23.*tst33);
  st23 = O21.*(O31.*tst11 + O32.*tst12 + O33.*tst13) + ...
                   O22.*(O31.*tst12 + O32.*tst22 + O33.*tst23) + ...
                   O23.*(O31.*tst13 + O32.*tst23 + O33.*tst33);
  st33 = O31.*(O31.*tst11 + O32.*tst12 + O33.*tst13) + ...
                   O32.*(O31.*tst12 + O32.*tst22 + O33.*tst23) + ...
                   O33.*(O31.*tst13 + O32.*tst23 + O33.*tst33);
  % Store fields
  ST = struct('st11',st11,'st12',st12,'st13',st13,...
              'st22',st22,'st23',st23,'st33',st33);

end
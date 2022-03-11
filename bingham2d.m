%-------------------------------------------------------------------------------
% Nematic Bingham closure in two dimensions
%
% Input:
%   Q - second moment tensor, struct with fields q11,q12,q22
%   T - rotation tensor T = E + 2*zeta*c*Q, struct with fields t11,t12,t22
%   cs1111 - Chebyshev coefficients mu1 --> ts1111
%
% Output:
%   ST - contraction S:T, struct with fields st11,st12,st22
%
% Scott Weady, CIMS
% Last updated February 2022
%-------------------------------------------------------------------------------
function ST = bingham2d(Q,T,cs1111)

  % Get components
  q11 = Q.q11; q12 = Q.q12+1e-15; %second moment tensor
  t11 = T.t11; t12 = T.t12; t22 = T.t22; %rotation tensor

  % Eigendecomposition of Q
  mu1 = min(0.5*(1 + sqrt((2*q11-1).^2 + 4*q12.^2)),1.0); %largest eigenvalue
  mu2 = 1 - mu1; %smallest eigenvalue
  ohm = 0.5*atan2(2*q12,2*q11-1); %angle of rotation matrix
  O11 = cos(ohm); O12 = -sin(ohm);
  O21 = -O12; O22 = O11;       

  % Evaluate Chebyshev interpolant
  M = length(cs1111);
  nu = 4*mu1 - 3;
  Tnm1 = 1; Tn = nu;
  ts1111 = cs1111(2)*Tn + cs1111(1)*Tnm1;            

  for k = 3:M
    Tnp1 = 2*nu.*Tn - Tnm1;
    ts1111 = ts1111 + cs1111(k)*Tnp1;
    Tnm1 = Tn; Tn = Tnp1;
  end              
  
  ts1122 = (mu1 - ts1111); 
  ts2222 = (mu2 - ts1122);
  
  % Rotate T to diagonal coordinate system
  tt11 = O11.*O11.*t11 + 2*O11.*O21.*t12 + O21.*O21.*t22;
  tt12 = O11.*O12.*t11 + (O11.*O22 + O21.*O12).*t12 + O21.*O22.*t22;
  tt22 = O12.*O12.*t11 + 2*O12.*O22.*t12 + O22.*O22.*t22;        

  % Compute contraction tS:tT
  tst11 = ts1111.*tt11 + ts1122.*tt22;
  tst12 = 2*ts1122.*tt12;
  tst22 = ts1122.*tt11 + ts2222.*tt22;

  % Rotate tS:tT to original coordinate system
  st11 = O11.*O11.*tst11 + 2.*O11.*O12.*tst12 + O12.*O12.*tst22;
  st12 = O11.*O21.*tst11 + (O11.*O22 + O21.*O12).*tst12 + O12.*O22.*tst22;
  st22 = O21.*O21.*tst11 + 2*O21.*O22.*tst12 + O22.*O22.*tst22;

  % Store fields
  ST = struct('st11',st11,'st12',st12,'st22',st22);

end
%% Simplified model for shaded leaves
% General model Exp2:
% Coefficients (with 95% confidence bounds):
       a =   0.0004953;%  (0.0004952, 0.0004954)
       b =    -0.04727;%  (-0.0473, -0.04723)
       c =   1.922e-06;%  (1.848e-06, 1.996e-06)
       d =     0.04127;%  (0.04063, 0.04191)
x       =   1:60;
     fx = a*exp(b*x) + c*exp(d*x);
% Goodness of fit:
%   SSE: 4.203e-13
%   R-square: 1
%   Adjusted R-square: 1
%   RMSE: 8.663e-08

%% Simplified model for sunlit leaves
%
figure     
subplot(2,1,1)
plot(rad.Pnh,x,fx,x,'x')
xlabel('radiation for shaded leaves [moles m-2 s-1]')
ylabel('layers')
legend('RTMo','simplified')
subplot(2,1,2)
plot((rad.Pnh-fx')./rad.Pnh*100,x,'x')
xlabel('error [%]')

figure
subplot(2,1,1)
plot(rad.Pnh,x,fx,x,'x')
xlabel('radiation for shaded leaves [moles m-2 s-1]')
ylabel('layers')
legend('RTMo','simplified')
subplot(2,1,2)
plot((rad.Pnh-fx')./rad.Pnh*100,x,'x')
xlabel('error [%]')


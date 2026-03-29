
function [Velrange, Vxyzdf, H, f] = leastSquareVel_TSA(vect, satvel, satclkRate, dop, pos, settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function calculates the Least Square Solution
%
% Inputs:
%   obs             - Observations for one epoch  卫星数和doppler
%   sat             - Satellite positions and velocities for one epoch  卫星速度位置
%   allSettings     - receiver settings
%   Vel             - receiver velocity and receiver clock drift 速度  自己
%   pos             - Initial position for the LSE 位置
%
% Outputs:
%   Vel             - receiver velocity and receiver clock drift
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set starting point


nmbOfSatellites = size(satvel, 2);
f = zeros(nmbOfSatellites,1);
H = ones(nmbOfSatellites,4);
Velrange = zeros(1, nmbOfSatellites);
Vxyzdf = zeros(4,1);
nmbOfIterations = 10;


% for iter = 1:nmbOfIterations
%     ind = 0;    
%         for channelNr = 1:nmbOfSatellites
% 
%                 relative_velocity(ind) = dx * sat.(signal).channel(channelNr).Vel(1) + dy * sat.(signal).channel(channelNr).Vel(2) + dz * sat.(signal).channel(channelNr).Vel(3);
%                 relative_velocity(ind) = relative_velocity(ind) / sqrt(dx*dx+dy*dy+dz*dz);
% 
%                 % Observed minus predicted
%                 omp.dRange_rate(ind) = pseudo_range_rate(ind) + relative_velocity(ind) + sat.(signal).channel(channelNr).Vel(4) * SPEED_OF_LIGHT;
% 
%                 Res(ind) = omp.dRange_rate(ind) - vel(3 + signalNr)*SPEED_OF_LIGHT;                
% 
% 
%         end
% end
% 
%     % This is the actual solutions to the LSE optimisation problem
%     clear H;
%     clear dR;    
%     H=sv_matrix;
%     dR=omp.dRange_rate;
%     DeltaVel=(H'*H)^(-1)*H'*dR';
% 
%     % Updating the position with the solution
%     vel(1)=DeltaVel(1);
%     vel(2)=DeltaVel(2);
%     vel(3)=DeltaVel(3);
% 
%     % Update the clock offsets for all systems
%     vel(4:end) = DeltaVel(4:end)/SPEED_OF_LIGHT;
% 
% 


for i = 1:nmbOfSatellites
    f(i) = -dop(i)* settings.c / (1575.42*1000000) - satvel(:,i)'*vect(:,i) + satclkRate(i) * settings.c;
    H(i,1:3) = - vect(:,i)';
end

% Vxyzdf = (H'*H)\H'*f;
Vxyzdf = H \ f;
for i = 1:nmbOfSatellites
    Velrange(i) = (satvel(:,i)-Vxyzdf(1:3))'*vect(:,i) ;
end


end
 
 

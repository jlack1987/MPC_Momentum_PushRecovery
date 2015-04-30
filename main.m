function main()

close all

constants = struct();
robot = struct();
dynamics = struct();
gait = struct();

constants.runningInMatlab = true;

robot.posLeftFootX = 0;
robot.posLeftFootY = 0.15;
robot.posRightFootX = 0;
robot.posRightFootY = -0.15;
robot.footWidth = 0.1;
robot.footLength = 0.15;
constants.T = 0.01;
constants.mass = 145;
constants.initialDoubleSupportDuration = 0.5;
constants.singleSupportDuration = 0.5;
constants.finalDoubleSupportDuration = 1;
constants.initialCenterOfMassHeight = 0.75;
constants.gravity = 9.81;
constants.finalDesiredCenterOfMassHeight = ... 
    constants.initialCenterOfMassHeight + 0.1;

constants.a = 1; % uCmp
constants.b = 1; % COP-COPRef
constants.c = 1e-4; % Ldot
constants.d = 1e-7; %uLdot

gait.comeToStop = true;
gait.nSteps = 1;
gait.stanceSide = 1;
gait.desiredStepLength = 0.25;
gait.desiredStepWidth = abs(robot.posLeftFootY - robot.posRightFootY);
gait.footSteps = cell(gait.nSteps + 2,1);
gait.footSteps{1} = [robot.posLeftFootX,robot.posLeftFootY,0];
gait.footSteps{2} = [robot.posRightFootX,robot.posRightFootY,0];
constants.Nds1 = round(constants.initialDoubleSupportDuration/constants.T);
constants.Nds2 = round(constants.finalDoubleSupportDuration/constants.T);
constants.Nss = round(constants.singleSupportDuration/constants.T);
constants.N = constants.Nds1 + gait.nSteps*(constants.Nds2 + constants.Nss) + ...
    (gait.nSteps-1)*constants.Nds2;
gait.t = 0:constants.T:constants.N*constants.T;

gait = createFootSteps(robot,gait);

gait.COPs = cell(length(gait.footSteps),1);

gait = planCOPs(gait);

% dynamics.initialConditionY = [0;0;0;0;0];
% dynamics.initialConditionX = [0;0;0;0;0];
dynamics.initialConditionY = [0;-0.05;0;0.8;-50];
dynamics.initialConditionX = [0;0.03;0;-1;30];

dynamics = generateComHeightDynamics(dynamics,constants);

dynamics = dynamics_lookahead(dynamics,constants,gait);

gait = setupDomainVectors(gait,constants);
gait = setupPcopRef(gait,constants);
gait.PcomRef = ones(constants.N,2);
gait.PcomRef(:,1) = gait.PcomRef(:,1)*gait.COPs{3}(1);
gait.PcomRef(:,2) = gait.PcomRef(:,2)*gait.COPs{3}(2);
gait = doMPC(constants,dynamics,robot,gait);

make_plots(robot,gait,dynamics,constants);
end

function gait = createFootSteps(robot,gait)
    if(gait.nSteps~=1)
        for i=3:3+gait.nSteps-1
            if(mod(i,2)~=0)
                if(i==3+gait.nSteps-1 && gait.comeToStop)
                    finalX = gait.footSteps{i-1}(1);
                else
                    finalX = robot.posLeftFootX + (i-2)*gait.desiredStepLength; 
                end
                gait.footSteps{i} = ...
                    [finalX,gait.desiredStepWidth/2,0];
            else
                if(i==3+gait.nSteps-1 && gait.comeToStop)
                    finalX = gait.footSteps{i-1}(1);
                else
                    finalX = robot.posRightFootX + (i-2)*gait.desiredStepLength; 
                end
                gait.footSteps{i} = ...
                    [finalX,-gait.desiredStepWidth/2,0];
            end
        end
    else
        gait.footSteps{3} = [gait.footSteps{2}(1)+gait.desiredStepLength,...
                             -sign(gait.footSteps{2}(2))*gait.desiredStepWidth/2,...
                             0];
    end
end

function gait = setupPcopRef(gait,constants)
gait.PcopRef = zeros(constants.N,3);
whichCOP = 2;
iter = 1;
for i = 1:gait.nSteps
    for j = 1:3
        gait.PcopRef = gait.PcopRef + gait.U(:,iter)*gait.COPs{whichCOP};
        if mod(j,2)==0
           whichCOP = whichCOP+1; 
        end
        iter = iter + 1;
    end
end

end

function gait = planCOPs(gait)

for i=1:length(gait.footSteps)
    if(length(gait.footSteps) > 3)
       if i==1
           gait.COPs{i} = (gait.footSteps{i} + gait.footSteps{i+1})/2;
       elseif i==length(gait.footSteps)
          gait.COPs{i} = (gait.footSteps{i} + gait.footSteps{i-1})/2;
       else
           gait.COPs{i} = gait.footSteps{i};
       end
    else 
       if i==1
           gait.COPs{i} = (gait.footSteps{i} + gait.footSteps{i+1})/2;
       elseif i==length(gait.footSteps)
          gait.COPs{i} = (gait.footSteps{i} + gait.footSteps{i-1})/2;
       else
           gait.COPs{i} = gait.footSteps{i};
       end
    end
end

end
function gait = setupDomainVectors(gait,constants)
gait.U = [];
U1 = ones(constants.Nds1,1);
U2 = ones(constants.Nss,1);
U3 = ones(constants.Nds2,1);
for i=1:3*gait.nSteps
   if mod(i,3)==1
       if(i == 1)
           gait.U = blkdiag(gait.U,U1);
       else
           gait.U = blkdiag(gait.U,U3);
       end
   elseif mod(i,3)==2
       gait.U = blkdiag(gait.U,U2);
   else
       gait.U = blkdiag(gait.U,U3);
   end
end
end

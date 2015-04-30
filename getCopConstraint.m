function [Acop,bcop] = getCopConstraint(gait,dynamics,robot,constants)
copConstraintX = zeros(constants.N,2);
copConstraintY = zeros(constants.N,2);

Nds1 = constants.Nds1;
Nss = constants.Nss;
Nds2 = constants.Nds2;
fw = (robot.footWidth/2)*0.8;
fl = (robot.footLength/2)*0.5;

PcopX = dynamics.PcopX*dynamics.initialConditionX;
PcopY = dynamics.PcopX*dynamics.initialConditionY;
PcopU = dynamics.PcopU;

x1 = gait.footSteps{1}(1);
x2 = gait.footSteps{2}(1);
x3 = gait.footSteps{3}(1);
y1 = gait.footSteps{1}(2);
y2 = gait.footSteps{2}(2);
y3 = gait.footSteps{3}(2);

copConstraintX(1:Nds1,1) = ones(Nds1,1)*(max(x1,x2)+fl);
copConstraintX(Nds1+1:Nds1+Nss,1) = ones(Nss,1)*(x2+fl);
copConstraintX(Nds1+Nss+1:Nds1+Nss+Nds2,1) = ones(Nds2,1)*max(x2,x3);
copConstraintX(1:Nds1,2) = ones(Nds1,1)*(-min(x1,x2)+fl);
copConstraintX(Nds1+1:Nds1+Nss,2) = ones(Nss,1)*(-x2+fl);
copConstraintX(Nds1+Nss+1:Nds1+Nss+Nds2,2) = ones(Nds2,1)*(-min(x2,x3));

copConstraintY(1:Nds1,1) = ones(Nds1,1)*(max(y1,y2)+fw);
copConstraintY(Nds1+1:Nds1+Nss,1) = ones(Nss,1)*(y2+fw);
copConstraintY(Nds1+Nss+1:Nds1+Nss+Nds2,1) = ones(Nds2,1)*(max(y2,y3));
copConstraintY(1:Nds1,2) = ones(Nds1,1)*(-min(y1,y2)+fw);
copConstraintY(Nds1+1:Nds1+Nss,2) = ones(Nss,1)*(-y2+fw);
copConstraintY(Nds1+Nss+1:Nds1+Nss+Nds2,2) = ones(Nds2,1)*(-min(y2,y3));

Acopx = [PcopU;-PcopU];
bcopx = [-PcopX + copConstraintX(:,1);
          PcopX + copConstraintX(:,2)];
   
Acopy = [PcopU;-PcopU];
bcopy = [-PcopX + copConstraintY(:,1);
          PcopX + copConstraintY(:,2)];

Acop = blkdiag(Acopx,Acopy);
bcop = [bcopx;bcopy];

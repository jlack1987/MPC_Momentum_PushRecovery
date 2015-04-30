function [Acmp,bcmp] = getCmpConstraint(gait,dynamics,robot,constants)
cmpConstraintX = zeros(constants.N,2);
cmpConstraintY = zeros(constants.N,2);

Nds1 = constants.Nds1;
Nss = constants.Nss;
Nds2 = constants.Nds2;
fw = (robot.footWidth/2)*0.8;
fl = (robot.footLength/2)*0.5;

PcopX = dynamics.PcopX*dynamics.initialConditionX;
PcopU = dynamics.PcopU;

LdotX = dynamics.LdotX*dynamics.initialConditionX;
LdotU = dynamics.LdotU;

LdotY = dynamics.LdotX*dynamics.initialConditionY;

zddot = dynamics.zddot(2:end);

for i = 1:size(LdotU,1)
   for j = 1:size(LdotU,2)
      LdotU(i,j) = LdotU(i,j)/(constants.mass*(zddot(i) + constants.gravity));
   end
   
   for j = 1:size(LdotX,2)
       LdotX(i,j) = LdotX(i,j)/(constants.mass*(zddot(i) + constants.gravity));
       LdotY(i,j) = LdotY(i,j)/(constants.mass*(zddot(i) + constants.gravity));
   end
end

PcmpU = PcopU + LdotU;
PcmpX = PcopX + LdotY;

PcmpY = PcopX + LdotX;

x1 = gait.footSteps{1}(1);
x2 = gait.footSteps{2}(1);
x3 = gait.footSteps{3}(1);
y1 = gait.footSteps{1}(2);
y2 = gait.footSteps{2}(2);
y3 = gait.footSteps{3}(2);

cmpConstraintX(1:Nds1,1) = ones(Nds1,1)*(max(x1,x2)+fl);
cmpConstraintX(Nds1+1:Nds1+Nss,1) = ones(Nss,1)*(x2+fl);
cmpConstraintX(Nds1+Nss+1:Nds1+Nss+Nds2,1) = ones(Nds2,1)*max(x2,x3);
cmpConstraintX(1:Nds1,2) = ones(Nds1,1)*(-min(x1,x2)+fl);
cmpConstraintX(Nds1+1:Nds1+Nss,2) = ones(Nss,1)*(-x2+fl);
cmpConstraintX(Nds1+Nss+1:Nds1+Nss+Nds2,2) = ones(Nds2,1)*(-min(x2,x3));

cmpConstraintY(1:Nds1,1) = ones(Nds1,1)*(max(y1,y2)+fw);
cmpConstraintY(Nds1+1:Nds1+Nss,1) = ones(Nss,1)*(y2+fw);
cmpConstraintY(Nds1+Nss+1:Nds1+Nss+Nds2,1) = ones(Nds2,1)*(max(y2,y3));
cmpConstraintY(1:Nds1,2) = ones(Nds1,1)*(-min(y1,y2)+fw);
cmpConstraintY(Nds1+1:Nds1+Nss,2) = ones(Nss,1)*(-y2+fw);
cmpConstraintY(Nds1+Nss+1:Nds1+Nss+Nds2,2) = ones(Nds2,1)*(-min(y2,y3));

Acmpx = [PcmpU;-PcmpU];
bcmpx = [-PcmpX + cmpConstraintX(:,1);
          PcmpX + cmpConstraintX(:,2)];
   
Acmpy = [PcmpU;-PcmpU];
bcmpy = [-PcmpY + cmpConstraintY(:,1);
          PcmpY + cmpConstraintY(:,2)];

Acmp = blkdiag(Acmpx,Acmpy);
bcmp = [bcmpx;bcmpy];
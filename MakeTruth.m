function [x,y,u] = MakeTruth(tStart,tEnd,numSteps,system,x0,uNom)
deltaT = (tEnd-tStart)/numSteps;
u_Nom = uNom;
u_hold = uNom;
x_Nom = [0,0,20,0,0,0]';
Qtrue = system.Qtrue;
Klin = system.Klin;
Rtrue = system.Rtrue;
tmp = x0;
NL_data = x0;

for k = 1:numSteps
    tSpan = [((k-1)*deltaT) k*deltaT];
    u(:,k) = u_Nom - Klin*(tmp-x_Nom);
    u_hold = u_Nom - Klin*(tmp-x_Nom);
    w = mvnrnd([0 0 0]',Qtrue,1)';
    [~,x_NL] = ode45(@(t,x)RHS(t,x,u_hold,w),tSpan,tmp);
    
    tmp = x_NL(end,:)';
    NL_data = [NL_data tmp];
    
    y_tmp = RHS(0,tmp,u_hold,w);
    v = mvnrnd([0 0 0 0]',Rtrue);
    y_data(:,k) = [tmp(1); tmp(3); tmp(6); y_tmp(2)]+v';
    
    
end
y = y_data;
x = NL_data;
end


close all;
clear;
tic
%********NOTE******** DE's must be SECOND ORDER

%Turn on adaptive convergence? (0 False or 1 True)
adaptive = 1;
%Convergence Rate
da=1;
noTimes = 5;
%initial Position
X1=-1000;
Y1=-500;
t1=0;
%Final Position
X2=200;
Y2=-100;
t2=20;
%How many points?
noPoints = 25;
t=linspace(t1,t2,noPoints);
%How many cycles to apply convergence?
noCycles = 100000;
endCycle = 0;
dt=t(2)-t(1);
disp(dt)
%X=(X2-X1)/(t2-t1)*(t-t1)+X1-0*(X2-X1)*(rand(1,noPoints)-0.5);
%Y=(Y2-Y1)/(t2-t1)*(t-t1)+Y1-0*(Y2-Y1)*(rand(1,noPoints)-0.5);
X=linspace(X1,X2,noPoints);
Y=linspace(Y1,Y2,noPoints);
plot(X,Y,"b*-")
figure
X(1)=X1;
X(end)=X2;
Y(1)=Y1;
Y(end)=Y2;
dx=0.00001;
dy=0.00001;
%trackda = zeros(noCycles,1); DEBUG
K1=zeros(length(Y)-2,length(Y)-2+length(X)-2);
K2=zeros(length(Y)-2,length(Y)-2+length(X)-2);
totalerror = zeros(noCycles,1);
times=0;
%Run for noCycles amount of times
f0 = zeros(length(Y)-2,1);
%Populate K1 matrix (jacobian)
for functionIndex = (1:(length(Y)-2))+1
    for varIndex = (1:(length(Y)-2+length(X)-2))+1
        if(varIndex-1 <= length(X)-2)
            xValIndex = varIndex;
            yValIndex=-100;
        elseif(varIndex-1 > length(X)-2)
            xValIndex=-100;
            yValIndex = varIndex-(length(X)-2);
        end
        %NOte i didn't change dx and dy
        K1(functionIndex-1,varIndex-1) = (getDerivFunction1Val(functionIndex,xValIndex,yValIndex,dx,dy,t,X,Y,dt)-getFunction1Val(functionIndex,t,X,Y,dt))/dy;
    end
end
for functionIndex = (1:(length(Y)-2))+1
    for varIndex = (1:(length(Y)-2+length(X)-2))+1
        if(varIndex-1 <= length(X)-2)
            xValIndex = varIndex;
            yValIndex=-100;
        elseif(varIndex-1 > length(X)-2)
            xValIndex=-100;
            yValIndex = varIndex-(length(X)-2);
        end
        K2(functionIndex-1,varIndex-1) = (getDerivFunction2Val(functionIndex,xValIndex,yValIndex,dx,dy,t,X,Y,dt)-getFunction2Val(functionIndex,t,X,Y,dt))/dy;
    end
end
%disp(size(K1))
%disp(size(K2))
K=[K1;K2];
%Populate function output vector
for functionIndex = 2:(length(Y)-1)
    f0(functionIndex-1) = getFunction1Val(functionIndex,t,X,Y,dt);
end
for functionIndex = length(Y):((length(Y)-2+length(X)-2)+1)
    f0(functionIndex-1) = getFunction2Val(functionIndex-(length(Y)-2),t,X,Y,dt);
end
%disp(K)
invK = inv(K);
%disp(det(K))
%disp(invK)
dY=-K\f0;
dY = dY*da/max(abs(dY));
%Seperate dY into the respective delta values
dX = dY(1:(length(X)-2));
dY = dY(length(X)-1:end);
%Update X and Y
x0 = X(2:end-1);
x1 = [X1;X(2:end-1)'+dX;X2]';
y0 = Y(2:end-1);
y1 = [Y1;Y(2:end-1)'+dY;Y2]';
f1 = zeros(length(Y)-2,1);
for functionIndex = 2:(length(Y)-1)
    f1(functionIndex-1) = getFunction1Val(functionIndex,t,x1,y1,dt);
end
for functionIndex = length(Y):((length(Y)-2+length(X)-2)+1)
    f1(functionIndex-1) = getFunction2Val(functionIndex-(length(Y)-2),t,x1,y1,dt);
end
x01 = [x0,y0];
x11 = [x1(2:end-1),y1(2:end-1)];
s=x11'-x01';
y=f1-f0;
invK = invK+((s-invK*y)*s'*invK)/(s'*invK*y);
x0=x1;
y0=y1;
f0=f1;
totalerror(1) = sum(abs(f0));
minError = totalerror(1);
bestSolutionx = 0;
bestSolutiony = 0;
for k = 2:noCycles
    %make sure K is invertible
    dY=-invK*f0;
    dY = dY*da/max(abs(dY));
    %Seperate dY into the respective delta values
    dX = dY(1:(length(X)-2));
    dY = dY(length(X)-1:end);
 %   disp(dX)
%disp(dY)
    %Update X and Y
    x1 = [X1;x0(2:end-1)'+dX;X2]';
    y1 = [Y1;y0(2:end-1)'+dY;Y2]';
    f1 = zeros(length(Y)-2,1);
    for functionIndex = 2:(length(Y)-1)
        f1(functionIndex-1) = getFunction1Val(functionIndex,t,x1,y1,dt);
    end
    for functionIndex = length(Y):((length(Y)-2+length(X)-2)+1)
        f1(functionIndex-1) = getFunction2Val(functionIndex-(length(Y)-2),t,x1,y1,dt);
    end
    x01 = [x0(2:end-1),y0(2:end-1)];
    x11 = [x1(2:end-1),y1(2:end-1)];
    %disp(size(x01));
    %disp(size(x11));
    s=x11'-x01';
    y=f1-f0;
    if((s'*invK*y) == 0)
        endCycle = k;
        break;
    end
    a=invK*y;
    invK = invK+((s-a)*s'*invK)/(s'*a);
    x0=x1;
    y0=y1;
    f0=f1;
    %Track error
    totalerror(k) = sum(abs(f0));
    if(totalerror(k) < minError)
        minError = totalerror(k);
        bestSolutionx = x0;
        bestSolutiony = y0;
    end
    %trackda(k) = da; DEBUG
    %Adaptive Convergence
    if(adaptive)
        if(k > 1)
            %If error is increasing, lower convergence rate
            if((totalerror(k)-totalerror(k-1))/totalerror(k-1) > 0.001)
                da = da*0.8;
            elseif((totalerror(k)-totalerror(k-1))/totalerror(k-1) < 0.001)
                %If error is decreasing steadily, increase convergence rate
                times = times+1;
                if(times == noTimes)
                    da = da*1.1;
                    times = 0;
                end
            end
        end
    end
end
disp("done!")
%Plot data
%disp(f0);
disp(sum(abs(f0)));
disp(sum(sqrt((diff(x0)/dt-x0(1:end-1)).^2+(diff(y0)/dt-y0(1:end-1)).^2))*dt);
disp(sum(sqrt((diff(bestSolutionx)/dt-bestSolutionx(1:end-1)).^2+(diff(bestSolutiony)/dt-bestSolutiony(1:end-1)).^2))*dt);
plot(x0,y0,"b*-");
xlabel("x")
ylabel("y")
title("Final Solution equation y(x)")
figure
plot(bestSolutionx,bestSolutiony,"b*-");
xlabel("x")
ylabel("y")
title("Best Solution equation y(x)")
figure
if(endCycle == 0)
    if(floor(noCycles/200) ~= 0)
        plot(1:floor(noCycles/200):noCycles,totalerror(1:floor(noCycles/200):end)/noPoints)
        %figure DEBUG
        %plot(1:floor(noCycles/200):noCycles,trackda(1:floor(noCycles/200):end))
    else
        plot(1:noCycles,totalerror/noPoints)
        %figure DEBUG
        %plot(1:noCycles,trackda)
    end
else
    if(floor(endCycle/200) ~= 0)
        plot(1:floor(endCycle/200):endCycle,totalerror(1:floor(endCycle/200):endCycle)/noPoints)
        %figure DEBUG
        %plot(1:floor(noCycles/200):noCycles,trackda(1:floor(noCycles/200):end))
    else
        plot(1:endCycle,totalerror/noPoints)
        %figure DEBUG
        %plot(1:noCycles,trackda)
    end
end
xlabel("Cycle")
ylabel("Average Error")
title("Average Error vs. Cycle")
toc
%figure
%plot(X(3:end-2),f0(2:end-1))
%xlabel("x")
%ylabel("Error")
%title("Error vs. x")


%Get value of f_i where i is functionIndex
function f1 = getFunction1Val(functionIndex,T,X,Y,dt)
y=getY(functionIndex,Y);
yp=getYPrime(functionIndex,Y,dt);
y2p=getY2Prime(functionIndex,Y,dt);
x=getX(functionIndex,X);
xp=getXPrime(functionIndex,X,dt);
x2p=getX2Prime(functionIndex,X,dt);
t=T(functionIndex);

v=sqrt((xp-x)^2+(yp-y)^2);

f1 = -(xp-x)/v-((x2p-xp)/v-(x'-x)*((xp-x)*(x2p-xp)+(yp-y)*(y2p-yp))/v^3);                                        %%Change Differential Equation Here (only allows y, y', and y'')

if(isnan(f1) || isinf(f1))
    f = 0;
end
%f = getY(functionIndex,Y)-getYPrime(functionIndex,Y,dt)/10;
end
%Get value of f_i where i in function index and we find how it changes when
%we increase y_j where j is yValIndex some amount. (To find d(f_i)/d(y_j),
% we need (f_i(y_j+deltaY)-f_i(y_j))/deltaY)
function f1 = getDerivFunction1Val(functionIndex,xValIndex,yValIndex,deltaX,deltaY,T,X,Y,dt)
y=getDerivY(functionIndex,yValIndex,deltaY,Y);
yp=getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt);
y2p=getDerivY2Prime(functionIndex,yValIndex,deltaY,Y,dt);
x=getDerivX(functionIndex,xValIndex,deltaX,X);
xp=getDerivXPrime(functionIndex,xValIndex,deltaX,X,dt);
x2p=getDerivX2Prime(functionIndex,xValIndex,deltaX,X,dt);
t=T(functionIndex);
v=sqrt((xp-x)^2+(yp-y)^2);
f1 = -(xp-x)/v-((x2p-xp)/v-(x'-x)*((xp-x)*(x2p-xp)+(yp-y)*(y2p-yp))/v^3);                                            %%Change Differential Equation Here must be the same as above

if(isnan(f1) || isinf(f1))
    f1 = 0;
end
%f = getDerivY(functionIndex,yValIndex,deltaY,Y)-getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt)/10;
end
%Get value of f_i where i is functionIndex
function f2 = getFunction2Val(functionIndex,T,X,Y,dt)
y=getY(functionIndex,Y);
yp=getYPrime(functionIndex,Y,dt);
y2p=getY2Prime(functionIndex,Y,dt);
x=getX(functionIndex,X);
xp=getXPrime(functionIndex,X,dt);
x2p=getX2Prime(functionIndex,X,dt);
t=T(functionIndex);
v=sqrt((xp-x)^2+(yp-y)^2);
f2 = -(yp-y)/v-((y2p-yp)/v-(y'-y)*((xp-x)*(x2p-xp)+(yp-y)*(y2p-yp))/v^3);                                        %%Change Differential Equation Here (only allows y, y', and y'')

if(isnan(f2) || isinf(f2))
    f2 = 0;
end
%f = getY(functionIndex,Y)-getYPrime(functionIndex,Y,dt)/10;
end
%Get value of f_i where i in function index and we find how it changes when
%we increase y_j where j is yValIndex some amount. (To find d(f_i)/d(y_j),
% we need (f_i(y_j+deltaY)-f_i(y_j))/deltaY)
function f2 = getDerivFunction2Val(functionIndex,xValIndex,yValIndex,deltaX,deltaY,T,X,Y,dt)
y=getDerivY(functionIndex,yValIndex,deltaY,Y);
yp=getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt);
y2p=getDerivY2Prime(functionIndex,yValIndex,deltaY,Y,dt);
x=getDerivX(functionIndex,xValIndex,deltaX,X);
xp=getDerivXPrime(functionIndex,xValIndex,deltaX,X,dt);
x2p=getDerivX2Prime(functionIndex,xValIndex,deltaX,X,dt);
t=T(functionIndex);
v=sqrt((xp-x)^2+(yp-y)^2);
f2 = -(yp-y)/v-((y2p-yp)/v-(y'-y)*((xp-x)*(x2p-xp)+(yp-y)*(y2p-yp))/v^3);                                              %%Change Differential Equation Here must be the same as above

if(isnan(f2) || isinf(f2))
    f2 = 0;
end
%f = getDerivY(functionIndex,yValIndex,deltaY,Y)-getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt)/10;
end

%Get y_i
function x = getX(functionIndex,X)
x = X(functionIndex);
end
%get y_j+deltaY
function xDeriv = getDerivX(functionIndex,xValIndex,deltaX,X)
xDeriv = X(functionIndex);
if(functionIndex == xValIndex)
    xDeriv = xDeriv+deltaX;
end
end
%Get y_i'
function xPrime = getXPrime(functionIndex,X,dt)
xi1=X(functionIndex+1);
xin1=X(functionIndex-1);
xPrime=(xi1-xin1)./(2*dt);
end
%Get y_j'
function xPrime = getDerivXPrime(functionIndex,xValIndex,deltaX,X,dt)
xi1=X(functionIndex+1);
xin1=X(functionIndex-1);
if(functionIndex+1==xValIndex)
    xi1 = xi1+deltaX;
elseif(functionIndex-1==xValIndex)
    xin1 = xin1+deltaX;
end
xPrime=(xi1-xin1)./(2*dt);
end
%Get y_i''
function x2Prime = getX2Prime(functionIndex,X,dt)
xi1=X(functionIndex+1);
xi=X(functionIndex);
xin1=X(functionIndex-1);
x2Prime=(xi1-2*xi+xin1)./(dt*dt);
end
%Get y_j''
function x2Prime = getDerivX2Prime(functionIndex,xValIndex,deltaX,X,dt)
xi1=X(functionIndex+1);
xi=X(functionIndex);
xin1=X(functionIndex-1);
if(functionIndex+1==xValIndex)
    xi1 = xi1+deltaX;
elseif(functionIndex-1==xValIndex)
    xin1 = xin1+deltaX;
elseif(functionIndex == xValIndex)
    xi=xi+deltaX;
end
x2Prime=(xi1-2*xi+xin1)./(dt*dt);
end


%Get y_i
function y = getY(functionIndex,Y)
y = Y(functionIndex);
end
%get y_j+deltaY
function yDeriv = getDerivY(functionIndex,yValIndex,deltaY,Y)
yDeriv = Y(functionIndex);
if(functionIndex == yValIndex)
    yDeriv = yDeriv+deltaY;
end
end
%Get y_i'
function yPrime = getYPrime(functionIndex,Y,dt)
yi1=Y(functionIndex+1);
yin1=Y(functionIndex-1);
yPrime=(yi1-yin1)./(2*dt);
end
%Get y_j'
function yPrime = getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt)
yi1=Y(functionIndex+1);
yin1=Y(functionIndex-1);
if(functionIndex+1==yValIndex)
    yi1 = yi1+deltaY;
elseif(functionIndex-1==yValIndex)
    yin1 = yin1+deltaY;
end
yPrime=(yi1-yin1)./(2*dt);
end
%Get y_i''
function y2Prime = getY2Prime(functionIndex,Y,dt)
yi1=Y(functionIndex+1);
yi=Y(functionIndex);
yin1=Y(functionIndex-1);
y2Prime=(yi1-2*yi+yin1)./(dt*dt);
end
%Get y_j''
function y2Prime = getDerivY2Prime(functionIndex,yValIndex,deltaY,Y,dt)
yi1=Y(functionIndex+1);
yi=Y(functionIndex);
yin1=Y(functionIndex-1);
if(functionIndex+1==yValIndex)
    yi1 = yi1+deltaY;
elseif(functionIndex-1==yValIndex)
    yin1 = yin1+deltaY;
elseif(functionIndex == yValIndex)
    yi=yi+deltaY;
end
y2Prime=(yi1-2*yi+yin1)./(dt*dt);
end


close all;
clear;
tic
%********NOTE******** DE's must be SECOND ORDER

%Turn on adaptive convergence? (0 False or 1 True)
adaptive = 1;
%Convergence Rate
da=0.0001;
noTimes = 5;
%How many points?
noPoints = 25;
%How many cycles to apply convergence?
noCycles = 10000;
endCycle = 0;
%Initial Time
t1=0;
%Final Time
t2=1;
%Create Time vector
t=linspace(t1,t2,noPoints);
dt=t(2)-t(1);
disp(dt)

%initial Position
X1=0.2;
Y1=0;

%Initial Velocity
Vx1=0;
Vy1=1;

X11 = X1-Vx1*dt;
Y11 = Y1-Vy1*dt;
%Final Position
X2=0;
Y2=1;
%Final Velocity
Vx2=-0.2;
Vy2=0;

X22 = X2+Vx2*dt;
Y22 = Y2+Vy2*dt;

%X=(X2-X1)/(t2-t1)*(t-t1)+X1-0*(X2-X1)*(rand(1,noPoints)-0.5);
%Y=(Y2-Y1)/(t2-t1)*(t-t1)+Y1-0*(Y2-Y1)*(rand(1,noPoints)-0.5);
Ax=(Vx1-(X2-X1)/(t2-t1))/(t1-t2);
Bx=3*(Vx2-(2*(X2-X1)/(t2-t1)-Vx1))/(t2-t1)^2;
Cx=-2*Bx*(t2-t1)/3;
X=(X2-X1)/(t2-t1)*(t-t1)+X1+Ax*((t-(t1+t2)/2).^2-(t1-t2)^2/4)+Bx*(t-t1).^3/3+Cx*(t-t1).^2/2;
Ay=(Vy1-(Y2-Y1)/(t2-t1))/(t1-t2);
By=3*(Vy2-(2*(Y2-Y1)/(t2-t1)-Vy1))/(t2-t1)^2;
Cy=-2*By*(t2-t1)/3;
Y=(Y2-Y1)/(t2-t1)*(t-t1)+Y1+Ay*((t-(t1+t2)/2).^2-(t1-t2)^2/4)+By*(t-t1).^3/3+Cy*(t-t1).^2/2;
X=[X11,X,X22];
Y=[Y11,Y,Y22];
%X=[X11,linspace(X1,X2,noPoints),X22];
%Y=[Y11,linspace(Y1,Y2,noPoints),Y22];
plot(X,Y,"b*-")
figure
dx=0.00001;
dy=0.00001;
%trackda = zeros(noCycles,1); DEBUG
K1=zeros(length(Y)-4,length(Y)-4+length(X)-4);
K2=zeros(length(Y)-4,length(Y)-4+length(X)-4);
totalerror = zeros(noCycles,1);
times=0;
%Run for noCycles amount of times
f0 = zeros(length(Y)-4,1);
%Populate K1 matrix (jacobian)
for functionIndex = (1:(length(Y)-4))+2
    for varIndex = (1:(length(Y)-4+length(X)-4))+2
        if(varIndex-2 <= length(X)-4)
            xValIndex = varIndex;
            yValIndex=-100;
        elseif(varIndex-2 > length(X)-4)
            xValIndex=-100;
            yValIndex = varIndex-(length(X)-4);
        end
        %NOte i didn't change dx and dy
        %disp(functionIndex+","+varIndex+","+xValIndex+","+yValIndex);
        K1(functionIndex-2,varIndex-2) = (getDerivFunction1Val(functionIndex,xValIndex,yValIndex,dx,dy,t,X,Y,dt)-getFunction1Val(functionIndex,dx,t,X,Y,dt))/dy;
    end
end
for functionIndex = (1:(length(Y)-4))+2
    for varIndex = (1:(length(Y)-4+length(X)-4))+2
        if(varIndex-2 <= length(X)-4)
            xValIndex = varIndex;
            yValIndex=-100;
        elseif(varIndex-2 > length(X)-4)
            xValIndex=-100;
            yValIndex = varIndex-(length(X)-4);
        end
        K2(functionIndex-2,varIndex-2) = (getDerivFunction2Val(functionIndex,xValIndex,yValIndex,dx,dy,t,X,Y,dt)-getFunction2Val(functionIndex,dy,t,X,Y,dt))/dy;
    end
end
%disp(size(K1))
%disp(size(K2))
K=[K1;K2];
%disp(K1)
%disp(K2)
%Populate function output vector
for functionIndex = (1:(length(Y)-4))+2
    f0(functionIndex-2) = getFunction1Val(functionIndex,dx,t,X,Y,dt);
end
for functionIndex = length(Y):((length(Y)-4+length(X)-4)+2)
    f0(functionIndex-2) = getFunction2Val(functionIndex-(length(Y)-4),dy,t,X,Y,dt);
end
%disp(K)
invK = inv(K);
%disp(det(K))
%disp(invK)
dY=-K\f0;
dY = dY*da/max(abs(dY));
%Seperate dY into the respective delta values
dX = dY(1:(length(X)-4));
dY = dY(length(X)-3:end);
%Update X and Y
x0 = X(3:end-2);
x1 = [X11;X1;X(3:end-2)'+dX;X2;X22]';
y0 = Y(3:end-2);
y1 = [Y11;Y1;Y(3:end-2)'+dY;Y2;Y22]';
f1 = zeros(length(Y)-4,1);
for functionIndex = (1:(length(Y)-4))+2
    f1(functionIndex-2) = getFunction1Val(functionIndex,dx,t,x1,y1,dt);
end
for functionIndex = length(Y):((length(Y)-4+length(X)-4)+2)
    f1(functionIndex-2) = getFunction2Val(functionIndex-(length(Y)-4),dy,t,x1,y1,dt);
end
x01 = [x0,y0];
x11 = [x1(3:end-2),y1(3:end-2)];
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
    dX = dY(1:(length(X)-4));
    dY = dY(length(X)-3:end);
 %   disp(dX)
%disp(dY)
    %Update X and Y
    x1 = [X11;X1;x0(3:end-2)'+dX;X2;X22]';
    y1 = [Y11;Y1;y0(3:end-2)'+dY;Y2;Y22]';
    f1 = zeros(length(Y)-4,1);
    for functionIndex = (1:(length(Y)-4))+2
        f1(functionIndex-2) = getFunction1Val(functionIndex,dx,t,x1,y1,dt);
    end
    for functionIndex = length(Y):((length(Y)-4+length(X)-4)+2)
        f1(functionIndex-2) = getFunction2Val(functionIndex-(length(Y)-4),dy,t,x1,y1,dt);
    end
    x01 = [x0(3:end-2),y0(3:end-2)];
    x11 = [x1(3:end-2),y1(3:end-2)];
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
x2p=(x0(3:end)-2*x0(2:end-1)+x0(1:end-2))/dt;
y2p=(y0(3:end)-2*y0(2:end-1)+y0(1:end-2))/dt;
G=6.67*10^-11;
M=6.67*10^11;
Fx=-G*M*(x0(2:end-1)./sqrt(x0(2:end-1).^2+y0(2:end-1).^2).^3);
Fy=-G*M*(y0(2:end-1)./sqrt(x0(2:end-1).^2+y0(2:end-1).^2).^3);
disp(sum(sqrt((x2p-Fx).^2+(y2p-Fy).^2)*dt));
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
function f1 = getFunction1Val(functionIndex,dx,T,X,Y,dt)

Px=getPxVal(functionIndex,X,Y,dt);
Py=getPyVal(functionIndex,X,Y,dt);

PxDot=getPxDotVal(functionIndex,X,Y,dt);
PyDot=getPyDotVal(functionIndex,X,Y,dt);

Px2Dot=getPx2DotVal(functionIndex,X,Y,dt);
Py2Dot=getPy2DotVal(functionIndex,X,Y,dt);

dPx_dx=getXDerivPxVal(functionIndex,dx,X,Y,dt);
dPy_dx=getXDerivPyVal(functionIndex,dx,X,Y,dt);

P=sqrt(Px^2+Py^2);
f1 = (Px*dPx_dx+Py*dPy_dx)/P+Px2Dot/P-2*PxDot*(Px*PxDot+Py*PyDot)/P^3-Px*((Px*Px2Dot+PxDot^2+Py*Py2Dot+PyDot^2)/P^3-3*(Px*PxDot+Py*PyDot)^2/P^5);                                        %%Change Differential Equation Here (only allows y, y', and y'')

if(isnan(f1) || isinf(f1))
    f1 = 0;
end
%f = getY(functionIndex,Y)-getYPrime(functionIndex,Y,dt)/10;
end
%Get value of f_i where i in function index and we find how it changes when
%we increase y_j where j is yValIndex some amount. (To find d(f_i)/d(y_j),
% we need (f_i(y_j+deltaY)-f_i(y_j))/deltaY)
function f1 = getDerivFunction1Val(functionIndex,xValIndex,yValIndex,deltaX,deltaY,T,X,Y,dt)
Px=getDerivPxVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
Py=getDerivPyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
%disp(Px-getPxVal(functionIndex,X,Y,dt));
%disp(Py-getPyVal(functionIndex,X,Y,dt));

PxDot=getDerivPxDotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
PyDot=getDerivPyDotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
%disp(PxDot-getPxDotVal(functionIndex,X,Y,dt));
%disp(PyDot-getPyDotVal(functionIndex,X,Y,dt));

Px2Dot=getDerivPx2DotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
Py2Dot=getDerivPy2DotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);

dPx_dx=getDerivXDerivPxVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
dPy_dx=getDerivXDerivPyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);

P=sqrt(Px^2+Py^2);
f1 = (Px*dPx_dx+Py*dPy_dx)/P+Px2Dot/P-2*PxDot*(Px*PxDot+Py*PyDot)/P^3-Px*((Px*Px2Dot+PxDot^2+Py*Py2Dot+PyDot^2)/P^3-3*(Px*PxDot+Py*PyDot)^2/P^5);                                        %%Change Differential Equation Here (only allows y, y', and y'')

if(isnan(f1) || isinf(f1))
    f1 = 0;
end
%f = getDerivY(functionIndex,yValIndex,deltaY,Y)-getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt)/10;
end
%Get value of f_i where i is functionIndex
function f2 = getFunction2Val(functionIndex,dy,T,X,Y,dt)
Px=getPxVal(functionIndex,X,Y,dt);
Py=getPyVal(functionIndex,X,Y,dt);

PxDot=getPxDotVal(functionIndex,X,Y,dt);
PyDot=getPyDotVal(functionIndex,X,Y,dt);

Px2Dot=getPx2DotVal(functionIndex,X,Y,dt);
Py2Dot=getPy2DotVal(functionIndex,X,Y,dt);

dPx_dy=getYDerivPxVal(functionIndex,dy,X,Y,dt);
dPy_dy=getYDerivPyVal(functionIndex,dy,X,Y,dt);

P=sqrt(Px^2+Py^2);
f2 = (Px*dPx_dy+Py*dPy_dy)/P+Py2Dot/P-2*PyDot*(Px*PxDot+Py*PyDot)/P^3-Py*((Px*Px2Dot+PxDot^2+Py*Py2Dot+PyDot^2)/P^3-3*(Px*PxDot+Py*PyDot)^2/P^5);                                        %%Change Differential Equation Here (only allows y, y', and y'')

if(isnan(f2) || isinf(f2))
    f2 = 0;
end
%f = getY(functionIndex,Y)-getYPrime(functionIndex,Y,dt)/10;
end
%Get value of f_i where i in function index and we find how it changes when
%we increase y_j where j is yValIndex some amount. (To find d(f_i)/d(y_j),
% we need (f_i(y_j+deltaY)-f_i(y_j))/deltaY)
function f2 = getDerivFunction2Val(functionIndex,xValIndex,yValIndex,deltaX,deltaY,T,X,Y,dt)
Px=getDerivPxVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
Py=getDerivPyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);

PxDot=getDerivPxDotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
PyDot=getDerivPyDotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);

Px2Dot=getDerivPx2DotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
Py2Dot=getDerivPy2DotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);

dPx_dy=getDerivYDerivPxVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);
dPy_dy=getDerivYDerivPyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,deltaX,deltaY);

P=sqrt(Px^2+Py^2);
f2 = (Px*dPx_dy+Py*dPy_dy)/P+Py2Dot/P-2*PyDot*(Px*PxDot+Py*PyDot)/P^3-Py*((Px*Px2Dot+PxDot^2+Py*Py2Dot+PyDot^2)/P^3-3*(Px*PxDot+Py*PyDot)^2/P^5);                                        %%Change Differential Equation Here (only allows y, y', and y'')

if(isnan(f2) || isinf(f2))
    f2 = 0;
end
%f = getDerivY(functionIndex,yValIndex,deltaY,Y)-getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt)/10;
end


function Px = getPxVal(functionIndex,X,Y,dt)
x2p=getX2Prime(functionIndex,X,dt);
Fx=getFxVal(functionIndex,X,Y);
Px = x2p-Fx;
end
function dPx_dx = getXDerivPxVal(functionIndex,dx,X,Y,dt)
Fx=getXDerivFxVal(functionIndex,X,dx,Y,dt);
dPx_dx = -Fx;
end
function dPx_dy = getYDerivPxVal(functionIndex,dy,X,Y,dt)
Fx=getYDerivFxVal(functionIndex,X,dy,Y,dt);
dPx_dy = -Fx;
end
function Px = getPxDotVal(functionIndex,X,Y,dt)
x3p=getX3Prime(functionIndex,X,dt);
FxDot=getFxDotVal(functionIndex,X,Y,dt);
Px = x3p-FxDot;
end
function Px = getPx2DotVal(functionIndex,X,Y,dt)
x4p=getX4Prime(functionIndex,X,dt);
Fx2Dot=getFx2DotVal(functionIndex,X,Y,dt);
Px = x4p-Fx2Dot;
end

function Px = getDerivPxVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
x2p=getDerivX2Prime(functionIndex,xValIndex,dx,X,dt);
Fx=getDerivFxVal(functionIndex,X,Y,xValIndex,yValIndex,dx,dy);
Px = x2p-Fx;
end
function dPx_dx = getDerivXDerivPxVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
Fx=getDerivXDerivFxVal(functionIndex,X,dx,Y,dt,xValIndex,yValIndex,dy);
dPx_dx = -Fx;
end
function dPx_dy = getDerivYDerivPxVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
Fx=getDerivYDerivFxVal(functionIndex,X,dx,Y,dt,xValIndex,yValIndex,dy);
dPx_dy = -Fx;
end
function Px = getDerivPxDotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
x3p=getDerivX3Prime(functionIndex,xValIndex,dx,X,dt);
%disp(x3p-getX3Prime(functionIndex,X,dt));
FxDot=getDerivFxDotVal(functionIndex,X,Y,dt,dx,dy,xValIndex,yValIndex);
Px = x3p-FxDot;
end
function Px = getDerivPx2DotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
x4p=getDerivX4Prime(functionIndex,xValIndex,dx,X,dt);
Fx2Dot=getDerivFx2DotVal(functionIndex,X,Y,dt,xValIndex,yValIndex,dx,dy);
Px = x4p-Fx2Dot;
end


function Py = getPyVal(functionIndex,X,Y,dt)
y2p=getY2Prime(functionIndex,Y,dt);
Fy=getFyVal(functionIndex,X,Y);
Py = y2p-Fy;
end
function dPy_dx = getXDerivPyVal(functionIndex,dx,X,Y,dt)
Fy=getXDerivFyVal(functionIndex,X,dx,Y,dt);
dPy_dx = -Fy;
end
function dPy_dy = getYDerivPyVal(functionIndex,dy,X,Y,dt)
Fy=getYDerivFyVal(functionIndex,X,dy,Y,dt);
dPy_dy = -Fy;
end
function Py = getPyDotVal(functionIndex,X,Y,dt)
y3p=getY3Prime(functionIndex,Y,dt);
FyDot=getFyDotVal(functionIndex,X,Y,dt);
Py = y3p-FyDot;
end
function Py = getPy2DotVal(functionIndex,X,Y,dt)
y4p=getY4Prime(functionIndex,Y,dt);
Fy2Dot=getFy2DotVal(functionIndex,X,Y,dt);
Py = y4p-Fy2Dot;
end

function Py = getDerivPyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
y2p=getDerivY2Prime(functionIndex,yValIndex,dy,Y,dt);
Fy=getDerivFyVal(functionIndex,X,Y,xValIndex,yValIndex,dx,dy);
Py = y2p-Fy;
end
function dPy_dx = getDerivXDerivPyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
Fy=getDerivXDerivFyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy);
dPy_dx = -Fy;
end
function dPy_dy = getDerivYDerivPyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
Fy=getDerivYDerivFyVal(functionIndex,X,dx,Y,dt,xValIndex,yValIndex,dy);
dPy_dy = -Fy;
end
function Py = getDerivPyDotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
y3p=getDerivY3Prime(functionIndex,yValIndex,dy,Y,dt);
%disp(y3p-getY3Prime(functionIndex,Y,dt));
FyDot=getDerivFyDotVal(functionIndex,X,Y,dt,dx,dy,xValIndex,yValIndex);
Py = y3p-FyDot;
end
function Py = getDerivPy2DotVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
y4p=getDerivY4Prime(functionIndex,yValIndex,dy,Y,dt);
Fy2Dot=getDerivFy2DotVal(functionIndex,X,Y,dt,xValIndex,yValIndex,dx,dy);
Py = y4p-Fy2Dot;
end

function Fx = getFxVal(functionIndex,X,Y)
G = 6.67*10^-11;
M = 6.67*10^11;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
Fx = -G*M*x/sqrt(x^2+y^2)^3;
end
function Fx = getDerivFxVal(functionIndex,X,Y,xValIndex,yValIndex,dx,dy)
G = 6.67*10^-11;
M = 6.67*10^11;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
if(functionIndex == xValIndex)
    x = x+dx;
end
if(functionIndex == yValIndex)
    y = y+dy;
end
Fx = -G*M*x/sqrt(x^2+y^2)^3;
end
function dFx_dx = getXDerivFxVal(functionIndex,X,dx,Y,dt)
G = 6.67*10^-11;
M = 6.67*10^11;
xpdx=getX(functionIndex,X)+dx;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
dFx_dx = -G*M*(xpdx/sqrt(xpdx^2+y^2)^3-x/sqrt(x^2+y^2)^3)/dx;
end
function dFx_dy = getYDerivFxVal(functionIndex,X,dy,Y,dt)
G = 6.67*10^-11;
M = 6.67*10^11;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
ypdy=getY(functionIndex,Y)+dy;
dFx_dy = -G*M*(x/sqrt(x^2+ypdy^2)^3-x/sqrt(x^2+y^2)^3)/dy;
end
function dFx_dx = getDerivXDerivFxVal(functionIndex,X,dx,Y,dt,xValIndex,yValIndex,dy)
G = 6.67*10^-11;
M = 6.67*10^11;
xpdx=getX(functionIndex,X)+dx;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
if(functionIndex == xValIndex)
    xpdx = xpdx+dx;
    x=x+dx;
end
if(functionIndex == yValIndex)
y=y+dy;
end
dFx_dx = -G*M*(xpdx/sqrt(xpdx^2+y^2)^3-x/sqrt(x^2+y^2)^3)/dx;
end
function dFx_dy = getDerivYDerivFxVal(functionIndex,X,dx,Y,dt,xValIndex,yValIndex,dy)
G = 6.67*10^-11;
M = 6.67*10^11;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
ypdy=getY(functionIndex,Y)+dy;
if(functionIndex == xValIndex)
x=x+dx;
end
if(functionIndex == yValIndex)
y=y+dy;
ypdy=ypdy+dy;
end
dFx_dy = -G*M*(x/sqrt(x^2+ypdy^2)^3-x/sqrt(x^2+y^2)^3)/dy;
end
function dFx_dt = getFxDotVal(functionIndex,X,Y,dt)
dFx_dt = 0;%(getFxVal(functionIndex+1,X,Y)-getFxVal(functionIndex-1,X,Y))/(2*dt);
end
function dFx_dt = getDerivFxDotVal(functionIndex,X,Y,dt,dx,dy,xValIndex,yValIndex)
dFx_dt = 0;%(getDerivFxVal(functionIndex+1,X,Y,xValIndex,yValIndex,dx,dy)-getDerivFxVal(functionIndex-1,X,Y,xValIndex,yValIndex,dx,dy))/(2*dt);
end
function d2Fx_dt2 = getFx2DotVal(functionIndex,X,Y,dt)
d2Fx_dt2 = 0;%(getFxVal(functionIndex+1,X,Y)-2*getFxVal(functionIndex,X,Y)+getFxVal(functionIndex-1,X,Y))/(dt^2);
end
function d2Fx_dt2 = getDerivFx2DotVal(functionIndex,X,Y,dt,xValIndex,yValIndex,dx,dy)
d2Fx_dt2 = 0;%(getDerivFxVal(functionIndex+1,X,Y,xValIndex,yValIndex,dx,dy)-2*getDerivFxVal(functionIndex,X,Y,xValIndex,yValIndex,dx,dy)+getDerivFxVal(functionIndex-1,X,Y,xValIndex,yValIndex,dx,dy))/(dt^2);
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
function x3Prime = getX3Prime(functionIndex,X,dt)
x2pi1 = getX2Prime(functionIndex+1,X,dt);
x2pin1 = getX2Prime(functionIndex-1,X,dt);
x3Prime = (x2pi1-x2pin1)/(2*dt);
end
%Get y_j'''
function x3Prime = getDerivX3Prime(functionIndex,xValIndex,deltaX,X,dt)
x2pi1 = getDerivX2Prime(functionIndex+1,xValIndex,deltaX,X,dt);
x2pin1 = getDerivX2Prime(functionIndex-1,xValIndex,deltaX,X,dt);
x3Prime = (x2pi1-x2pin1)/(2*dt);
end
%Get y_i''''
function x4Prime = getX4Prime(functionIndex,X,dt)
x2pi1 = getX2Prime(functionIndex+1,X,dt);
x2p = getX2Prime(functionIndex,X,dt);
x2pin1 = getX2Prime(functionIndex-1,X,dt);
x4Prime = (x2pi1-2*x2p+x2pin1)/(dt^2);
end
%Get y_j''''
function x4Prime = getDerivX4Prime(functionIndex,xValIndex,deltaX,X,dt)
x2pi1 = getDerivX2Prime(functionIndex+1,xValIndex,deltaX,X,dt);
x2p = getDerivX2Prime(functionIndex,xValIndex,deltaX,X,dt);
x2pin1 = getDerivX2Prime(functionIndex-1,xValIndex,deltaX,X,dt);
x4Prime = (x2pi1-2*x2p+x2pin1)/(dt^2);
end

function Fy = getFyVal(functionIndex,X,Y)
G = 6.67*10^-11;
M =6.67*10^11;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
Fy = -G*M*y/sqrt(x^2+y^2)^3;
end
function Fy = getDerivFyVal(functionIndex,X,Y,xValIndex,yValIndex,dx,dy)
G = 6.67*10^-11;
M = 6.67*10^11;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
if(functionIndex == xValIndex)
    x = x+dx;
end
if(functionIndex == yValIndex)
    y = y+dy;
end
Fy = -G*M*y/sqrt(x^2+y^2)^3;
end
function dFy_dy = getYDerivFyVal(functionIndex,X,dy,Y,dt)
G = 6.67*10^-11;
M =6.67*10^11;
x=getX(functionIndex,X);
ypdy=getY(functionIndex,Y)+dy;
y=getY(functionIndex,Y);
dFy_dy = -G*M*(ypdy/sqrt(x^2+ypdy^2)^3-y/sqrt(x^2+y^2)^3)/dy;
end
function dFy_dx = getXDerivFyVal(functionIndex,X,dx,Y,dt)
G = 6.67*10^-11;
M = 6.67*10^11;
xpdx=getX(functionIndex,X)+dx;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
dFy_dx = -G*M*(y/sqrt(xpdx^2+y^2)^3-y/sqrt(x^2+y^2)^3)/dx;
end
function dFy_dy = getDerivYDerivFyVal(functionIndex,X,dx,Y,dt,xValIndex,yValIndex,dy)
G = 6.67*10^-11;
M =6.67*10^11;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
ypdy=getY(functionIndex,Y)+dy;
if(functionIndex == xValIndex)
    x=x+dx;
end
if(functionIndex == yValIndex)
    ypdy=ypdy+dy;
    y=y+dy;
end
dFy_dy = -G*M*(ypdy/sqrt(x^2+ypdy^2)^3-y/sqrt(x^2+y^2)^3)/dx;
end
function dFy_dx = getDerivXDerivFyVal(functionIndex,xValIndex,yValIndex,X,Y,dt,dx,dy)
G = 6.67*10^-11;
M = 6.67*10^11;
xpdx=getX(functionIndex,X)+dx;
x=getX(functionIndex,X);
y=getY(functionIndex,Y);
if(functionIndex == xValIndex)
    x=x+dx;
    xpdx=xpdx+dx;
end
if(functionIndex == yValIndex)
    y=y+dy;
end
dFy_dx = -G*M*(y/sqrt(xpdx^2+y^2)^3-y/sqrt(x^2+y^2)^3)/dx;
end
function dFy_dt = getFyDotVal(functionIndex,X,Y,dt)
dFy_dt = 0;%(getFyVal(functionIndex+1,X,Y)-getFyVal(functionIndex-1,X,Y))/(2*dt);
end
function dFy_dt = getDerivFyDotVal(functionIndex,X,Y,dt,dx,dy,xValIndex,yValIndex)
dFy_dt = 0;%(getDerivFyVal(functionIndex+1,X,Y,xValIndex,yValIndex,dx,dy)-getDerivFyVal(functionIndex-1,X,Y,xValIndex,yValIndex,dx,dy))/(2*dt);
end
function d2Fy_dt2 = getFy2DotVal(functionIndex,X,Y,dt)
d2Fy_dt2 = 0;%(getFyVal(functionIndex+1,X,Y)-2*getFyVal(functionIndex,X,Y)+getFyVal(functionIndex-1,X,Y))/(dt^2);
end
function d2Fy_dt2 = getDerivFy2DotVal(functionIndex,X,Y,dt,xValIndex,yValIndex,dx,dy)
d2Fy_dt2 = 0;%(getDerivFyVal(functionIndex+1,X,Y,xValIndex,yValIndex,dx,dy)-2*getDerivFyVal(functionIndex,X,Y,xValIndex,yValIndex,dx,dy)+getDerivFyVal(functionIndex-1,X,Y,xValIndex,yValIndex,dx,dy))/(dt^2);
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
%Get y_i'''
function y3Prime = getY3Prime(functionIndex,Y,dt)
y2pi1 = getY2Prime(functionIndex+1,Y,dt);
y2pin1 = getY2Prime(functionIndex-1,Y,dt);
y3Prime = (y2pi1-y2pin1)/(2*dt);
end
%Get y_j'''
function y3Prime = getDerivY3Prime(functionIndex,yValIndex,deltaY,Y,dt)
y2pi1 = getDerivY2Prime(functionIndex+1,yValIndex,deltaY,Y,dt);
y2pin1 = getDerivY2Prime(functionIndex-1,yValIndex,deltaY,Y,dt);
y3Prime = (y2pi1-y2pin1)/(2*dt);
end
%Get y_i''''
function y4Prime = getY4Prime(functionIndex,Y,dt)
y2pi1 = getY2Prime(functionIndex+1,Y,dt);
y2p = getY2Prime(functionIndex,Y,dt);
y2pin1 = getY2Prime(functionIndex-1,Y,dt);
y4Prime = (y2pi1-2*y2p+y2pin1)/(dt^2);
end
%Get y_j''''
function y4Prime = getDerivY4Prime(functionIndex,yValIndex,deltaY,Y,dt)
y2pi1 = getDerivY2Prime(functionIndex+1,yValIndex,deltaY,Y,dt);
y2p = getDerivY2Prime(functionIndex,yValIndex,deltaY,Y,dt);
y2pin1 = getDerivY2Prime(functionIndex-1,yValIndex,deltaY,Y,dt);
y4Prime = (y2pi1-2*y2p+y2pin1)/(dt^2);
end

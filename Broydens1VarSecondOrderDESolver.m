close all;
clear;
%********NOTE******** DE's must be SECOND ORDER

%Turn on adaptive convergence? (0 False or 1 True)
adaptive = 0; 
%Convergence Rate
da=0.001;
%initial Position
Y1=0;
X1=0;
%Final Position
X2=1;
Y2=0;
%How many points?
noPoints = 100;
X=linspace(X1,X2,noPoints);
%How many cycles to apply convergence?
noCycles = 1000;
endCycle = 0;
dt=X(2)-X(1);
Y=(Y2-Y1)/(X2-X1)*(X-X1)+Y1-0*(Y2-Y1)*(rand()-0.5);
Y(1)=Y1;
Y(end)=Y2;
dy=0.00001;


%trackda = zeros(noCycles,1); DEBUG
K=zeros(length(Y)-2,length(Y)-2);
totalerror = zeros(noCycles,1);
times=0;
%Run for noCycles amount of times
f0 = zeros(length(Y)-2,1);
%Populate K matrix (jacobian)
for functionIndex = (1:(length(Y)-2))+1
    for yValIndex = (1:(length(Y)-2))+1
        K(functionIndex-1,yValIndex-1) = (getDerivFunctionVal(functionIndex,yValIndex,dy,X,Y,dt)-getFunctionVal(functionIndex,X,Y,dt))/dy;
    end
end
%Populate function output vector
for functionIndex = 2:(length(Y)-1)
    f0(functionIndex-1) = getFunctionVal(functionIndex,X,Y,dt);
end
invK = inv(K);
dY=-K\f0;
dY = dY*da/max(abs(dY));
%Update Y
x0 = Y(2:end-1);
x1 = [Y1;Y(2:end-1)'+dY;Y2]';
f1 = zeros(length(Y)-2,1);
for functionIndex = 2:(length(Y)-1)
    f1(functionIndex-1) = getFunctionVal(functionIndex,X,x1,dt);
end
s=x1(2:end-1)'-x0';
y=f1-f0;
invK = invK+((s-invK*y)*s'*invK)/(s'*invK*y);
x0=x1;
f0=f1;
totalerror(1) = sum(abs(f0));
minError = totalerror(1);
bestSolution = 0;
for k = 2:noCycles
    %make sure K is invertible
    dY=-invK*f0;
    dY = dY*da/max(abs(dY));
    %Update x1
    x1 = [Y1;x0(2:end-1)'+dY;Y2]';
    %update f1
    for functionIndex = 2:(length(x0)-1)
        f1(functionIndex-1) = getFunctionVal(functionIndex,X,x1,dt);
    end
    s=(x1(2:end-1)-x0(2:end-1))';
    y=f1-f0;
    if((s'*invK*y) == 0)
        endCycle = k;
        break;
    end
    a = invK*y;
    invK = invK+((s-a)*s'*invK)/(s'*a);
    x0=x1;
    f0=f1;
    %Track error
    totalerror(k) = sum(abs(f0));
    if(totalerror(k) < minError)
        minError = totalerror(k);
        bestSolution = x0;
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
                if(times == 3)
                    da = da*1.1;
                    times = 0;
                end
            end
        end
    end
end
disp("done!")
%Plot data
plot(X,x0);
xlabel("x")
ylabel("y")
title("Final Solution equation y(x)")
figure
plot(X,bestSolution);
xlabel("x")
ylabel("y")
title("Best Solution equation y(x)")
figure
if(endCycle == 0)
    if(floor(noCycles/200) ~= 0)
        plot(1:floor(noCycles/200):noCycles,totalerror(1:floor(noCycles/200):end)/length(X))
        %figure DEBUG
        %plot(1:floor(noCycles/200):noCycles,trackda(1:floor(noCycles/200):end))
    else
        plot(1:noCycles,totalerror)
        %figure DEBUG
        %plot(1:noCycles,trackda)
    end
else
    if(floor(endCycle/200) ~= 0)
        plot(1:floor(endCycle/200):endCycle,totalerror(1:floor(endCycle/200):endCycle)/length(X))
        %figure DEBUG
        %plot(1:floor(noCycles/200):noCycles,trackda(1:floor(noCycles/200):end))
    else
        plot(1:endCycle,totalerror)
        %figure DEBUG
        %plot(1:noCycles,trackda)
    end
end
xlabel("Cycle")
ylabel("Average Error")
title("Average Error vs. Cycle")
figure
plot(X(3:end-2),f0(2:end-1))
xlabel("x")
ylabel("Error")
title("Error vs. x")


%Get value of f_i where i is functionIndex
function f = getFunctionVal(functionIndex,X,Y,dt)
y=getY(functionIndex,Y);
yp=getYPrime(functionIndex,Y,dt);
y2p=getY2Prime(functionIndex,Y,dt);
x=X(functionIndex);
f = sqrt(1-yp^2)-y2p;                                        %%Change Differential Equation Here (only allows y, y', and y'')
if(isnan(f) || isinf(f))
    f = 0;
end
%f = getY(functionIndex,Y)-getYPrime(functionIndex,Y,dt)/10;
end
%Get value of f_i where i in function index and we find how it changes when
%we increase y_j where j is yValIndex some amount. (To find d(f_i)/d(y_j), 
% we need (f_i(y_j+deltaY)-f_i(y_j))/deltaY)
function f = getDerivFunctionVal(functionIndex,yValIndex,deltaY,X,Y,dt)
y=getDerivY(functionIndex,yValIndex,deltaY,Y);
yp=getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt);
y2p=getDerivY2Prime(functionIndex,yValIndex,deltaY,Y,dt);
x=X(functionIndex);
f = sqrt(1-yp^2)-y2p;                                            %%Change Differential Equation Here must be the same as above
if(isnan(f) || isinf(f))
    f = 0;
end
%f = getDerivY(functionIndex,yValIndex,deltaY,Y)-getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt)/10;
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


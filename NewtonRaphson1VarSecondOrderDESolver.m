close all;
clear;
%********NOTE******** DE's must be SECOND ORDER

%Turn on adaptive convergence? (0 False or 1 True)
adaptive = 0; 
%Convergence Rate
da=0.1;
%initial Position
Y1=0;
X1=0;
%Final Position
X2=4;
Y2=4;
%How many points?
noPoints = 100;
X=linspace(X1,X2,noPoints);
%How many cycles to apply convergence?
noCycles = 100;

dt=X(2)-X(1);
Y=(Y2-Y1)/(X2-X1)*(X-X1)+Y1+(Y2-Y1)*(rand()-0.5);
Y(1)=Y1;
Y(end)=Y2;
dy=0.01;

%trackda = zeros(noCycles,1); DEBUG
K=zeros(length(Y)-2,length(Y)-2);
totalerror = zeros(noCycles,1);
times=0;
%Run for noCycles amount of times
for k = 1:noCycles
    %initialize f
    f = zeros(length(Y)-2,1);
    %Populate K matrix (jacobian)
    for functionIndex = (1:(length(Y)-2))+1
        for yValIndex = (1:(length(Y)-2))+1
            K(functionIndex-1,yValIndex-1) = (getDerivFunctionVal(functionIndex,yValIndex,dy,Y,dt)-getFunctionVal(functionIndex,Y,dt))/dy;
        end
    end
    %Populate function output vector
    for functionIndex = 2:(length(Y)-1)
        f(functionIndex-1) = getFunctionVal(functionIndex,Y,dt);
    end
    %make sure K is invertible
    if(det(K) ~= 0)
        dY=-K\f; %-inv(K)*f;
        dY = dY*da/max(abs(dY));
        %Update Y
        Y = [Y1;Y(2:end-1)'+dY;Y2]';
    else
        %If K is not invertible, we can't do anything
        break;
    end
    %Track error
    totalerror(k) = sum(abs(f));
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
%Plot data
plot(X,Y);
xlabel("x")
ylabel("y")
title("Solution equation y(x)")
figure
if(floor(noCycles/200) ~= 0)
    plot(1:floor(noCycles/200):noCycles,totalerror(1:floor(noCycles/200):end)/length(X))
    %figure DEBUG
    %plot(1:floor(noCycles/200):noCycles,trackda(1:floor(noCycles/200):end))
else
    plot(1:noCycles,totalerror)
    %figure DEBUG
    %plot(1:noCycles,trackda)
end
xlabel("Cycle")
ylabel("Average Error")
title("Average Error vs. Cycle")
figure
plot(X(3:end-2),f(2:end-1))
xlabel("x")
ylabel("Error")
title("Error vs. x")


%Get value of f_i where i is functionIndex
function f = getFunctionVal(functionIndex,Y,dt)
y=getY(functionIndex,Y);
yp=getYPrime(functionIndex,Y,dt);
y2p=getY2Prime(functionIndex,Y,dt);
f = 20*y+y2p;                                           %%Change Differential Equation Here (only allows y, y', and y'')
if(isnan(f) || isinf(f))
    f = 0;
end
%f = getY(functionIndex,Y)-getYPrime(functionIndex,Y,dt)/10;
end
%Get value of f_i where i in function index and we find how it changes when
%we increase y_j where j is yValIndex some amount. (To find d(f_i)/d(y_j), 
% we need (f_i(y_j+deltaY)-f_i(y_j))/deltaY)
function f = getDerivFunctionVal(functionIndex,yValIndex,deltaY,Y,dt)
y=getDerivY(functionIndex,yValIndex,deltaY,Y);
yp=getDerivYPrime(functionIndex,yValIndex,deltaY,Y,dt);
y2p=getDerivY2Prime(functionIndex,yValIndex,deltaY,Y,dt);
f = 20*y+y2p;                                            %%Change Differential Equation Here must be the same as above
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


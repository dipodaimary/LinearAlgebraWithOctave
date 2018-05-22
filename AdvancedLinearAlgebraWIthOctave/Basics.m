% f = input signal (1xa)
% h = impulse response (1xb)
% g = output signal (1x(a + b-1))
function [g] = convm1d(f,h)
    a = length(h);
    b = length(f);
    n = a + b-1;
    he = [h zeros(1,n-a)];
    fe = [f zeros(1,n-b)]';
    h1 = [he(1) he(n:-1:2)];
    H = h1;
    for i = 1:n-1
        h1 = [h1(end) h1(1:end-1)];
        H = [H;h1];
    end
    g = H*fe;
    g = g'
end

A = [1 2; 3 4; 5 6]
size(A)
size(A,1)
size(A,2)
v = [1 2 3 4]
length(v)
length(A)
length([1;2;3;4])
%load('featuresX.dat') or load featuresX.dat
clear featuresX;
v = priceY(1:10) %first 10 elements
clear
whos
whos
load hello.mat
save hello.txt v -ascii
A(3,2)
A(2,:)
A(:,2)
A([1,3],:)
A(:,2)
A(:,2) = [10;11;12]
A = [A,[101,102]];
A(:)
A = [1 2; 3 4; 5 6]
C = [A B]
A = [A; B]
%Calculations
A.*B
1./A
log(v)
exp(v)
abs(v)
Y+1
find(a<3)
a<3
magic(3)
floor(u)
ceil(a)
rand(3)
max(rand(3),rand(3))
max(A,[],1) %col wise maximum
max(A,[],2) %row wise maximum
max(A) 
max(A(:))
A = magic(9)
sum(a,1) %col wise sum
sum(A,2) %row wise sum
sum(A.*eye(9)) %diagonal sum
flipupd(eye(9)) %upside down
pinv(magic(3)) %pseudo inverse
image.sc(magic(15)), colorbar, colormap grey

% control statements
for i = 1:10,
    v(i)  = 2^i; %disp(i)
end;

i=1
while i<=5,
    v(i) = 100;
    i = i+1;
end

while true,
    v(i) = 100
    i = i+1
    if(i==6),
        break;
    end
end

v(1) = 2
if v(1)==1;
    disp('The value is one')
elseif v(1)==2;
    disp('The value is two')
else
    disp('The value is not one or two')
%octave search path
addPath('/home/..')

%functions
function [y1,y2] = squareAndCubeThisNumber(x)
    y1 = x^2
    y2 = x^3

X = [1 1; 1 2; 1 3]
y = [1;2;3]
theta = [0,3]

%CostFunctionJ.m
function J = CostFunctionJ(X,y,theta)
    m = size(X,1)
    predictions = X*theta
    %examples
    sqrErrors = (predictions-y).^2  %squared errors 
    J = 1/(2*m)*sum(sqrErrors)

%plots
t = [0:0.01:0.98];
y1 = sin(2*pi*4*t);
y2 = cos(2*pi*4*t);
plot(t,y1);
hold on ;
plot(t,y2,'r');
legend('sin','cos');
title('myplot');
print -dpng '/tmp/pylot.png';
figure(1); plot(t,y1)
figure(2); plot(t,y2)

subplot(1,2,1) %divide the plot to two grids
plot(t,y1)
subplot(1,2,2)
plot(t,y2)
axis([0.5 1 -1 1])

A = magic(15)
imagesc(A), colorbar, colormap, gray;

%Control statements
for i=1:10,
    v(i) = 2^i;
end;
v
%
i =1
while(i<=5),
    v(i)=100;
    i=i+1;
end;
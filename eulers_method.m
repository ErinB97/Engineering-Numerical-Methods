% Euler's Method
% A quick demonstration of Euler's method to solve differntial equations
% based on https://www.khanacademy.org/math/ap-calculus-bc/bc-differential-equations-new/bc-7-5/v/eulers-method

clear all
close all


%differential to be solved
dydx = @(y) y; 
%known analytical solution
ysol = @(x) exp(x);

%IV's and structures
x0 = 0; %IV
delx = 0.5; %step-size
xf = 3;

X = x0:delx:xf; %series of X values to solve across
Y = zeros(1,length(X)); %results
Y(1) = 1; %IV

%error
err = zeros(1,length(X));
err(1) = ( ysol( X(1) ) - Y(1) )/ ysol(X(1));

%Euler's method
for i = 2:1:length(X)
    
    Y(i) = Y(i-1) + dydx( Y(i-1) )*delx; %euler's method
    
    err(i) = ( ysol( X(i) ) - Y(i) )/ ysol(X(i)); %calculating error
    
end

%displaying the results
figure; title("e^x vs. Euler's approximation");hold on;xlabel('X');ylabel('Y');grid on
scatter(X,Y,'b','x');scatter([x0:delx:xf],exp(x0:delx:xf),'r','x');

fprintf("X step-size is %.3f \n", delx);
fprintf("Maximum error is %.3f \n", max(err));



%Newton-Raphson Method University Assignment
clear all 

a = 3; %From assignment

%Delcaring the function F and diffierntial dF/dx
func = @(x) (x.^5 - 16.5.*x.^4 + 102.85.*x.^3 - 299.475.*x.^2 + 401.1634.*x -193.5612 - (a/10) ) ; 
dfunc = @(x) (5.*x.^4 - 66.*x.^3 + (6171.*x.^2)/20 - (11979.*x)/20 +  401.1634);

disp('This is a Newton-Rasphon calculator for the function');
disp(func);
fprintf('Where a = %d \n',a);

[imax,emax,dmin,xi] = GetUser(); %User inputs: Maximum iterations, maximum
% error, minimum gradient (for catching zero division), and their initial
% guess for the root location.


%Newton-Raphson Method
%!!! Give this more step-by-step commenting and better variable names!

while i < imax   %i:current iteration, loops until i = imax.
    
   f = func(xi);
   f1 = dfunc(xi);
   
   if abs(f1) < dmin %Stops the loop when the minimum gradient condition is met.
        break        
   end   
    
   xipp = xi - f/f1; %Calculates the next value for Xroot
   
   deltax = abs(xipp-xi/xipp);
   
   if deltax <= emax      
       break    %Stops the loop when the error condition is met.
   end
   
   xi = xipp;
   
   i = i+1;    
end

%Returning the results and calculation steps
fprintf('The root was calculated to be at Xr = %f4 \n',xi'); 
fprintf('The root was calculated in %d iterations \n',i'); 
fprintf('The final error was err = %f1 \n',deltax'); 


function [imax,emax,dmin,xi] = GetUser()
%This function simply obtains the user input for max iteration, maximum
%error, minimum gradient (for catching division by zero), and user's initial
%guess of X. 

    pass = 0; %used to keep user in loop if they enter a negative value 
    
    %User inputs maximum solve iterations
    while pass == 0
      imax = input('Maximum Iteration: ');  
        if imax <= 0 
         disp('Enter a +ve value.');
        else
         pass = 1;
        end

    end

    pass = 0;
    %User inputs maximum error tolerance
    while pass == 0
      dmin = input('Maximum Error:  ');  
        if dmin <= 0 
         disp('Enter a +ve value.');
        else 
         pass = 1;   
        end    
    end

    pass = 0;
    %User imputs the mimum graident allowed
    while pass == 0
      emax = input('Minmum Gradient (-> 0 fails): ');  
        if emax <= 0 
         disp('Enter a +ve value.'.');
        else 
         pass = 1;
        end  
    end

    pass = 0;
    %User enters their intiial guess for the X-value of the function root 
    while pass == 0
      xi = input('Initial root of function (X) : ');  
        if dmin <= 0 
         disp('Enter a +ve real value.');
        else 
         pass = 1;   
        end    
    end


end


% Gauss Seidel Method GB00365 URN: 6411353
clear all
close all
global C_Err ES nmax  nc h

%-----------------------------------------------------------
% Matrix and vector input by user, which is then checked for problems
 disp('Enter the co-efficients of your system of equations as a square matrix: ');
 A = input('');
 
 h = size(A,1); %The height of the matrix = n

 
 if  CheckZeros(A) == 0  %See error for description
    error('You entered an equation with a diagonal element of 0, it cannot be re-arranged to a diagonal dominant form. Try entering the equations in a different order, or double check you entered them correctly.'); 
 end     

disp('Enter the results of your system of equations as a column vector: ');
B = input('');

 if CheckSquare(A,B) == 0 %See error for description
   error('You have not entered the same number of equations as variables to solve. You must enter a square matrix, and its size must equal the length of the vector B'); 
 end    



%-----------------------------------------------------------
%Check that Matrix is diagonally dominant, if not, re-arrange.
if CheckDiagB(A) == 0
 disp('The matrix you entered was not diagonally dominant.');   
 [A,B] = MakeDiag(A,B); 

end 

%
%-----------------------------------------------------------
%User inputs initial guesses, error, and maximum iterations. 
lab = [1:h]; %Labels for X1,X2 etc..
guide = ('X%1.0f'); %Produces a format for X1,X2,...,Xh

X = GetGuess(); %User inputs initial values
 
[nmax, ES] = GetErrIt();

C_Err = [100 100 100 100]; %Initial error, holds no important value.
  
%-----------------------------------------------------------
%Finds solutions, then displays the resultsa and solution information.

X = Gauss_Seidel(X,A,B);
Results = X(size(X,1),:);

disp('The system of equations has been solved to produce the results:');


head = ('    X%1.0f    '); %Displaying the results in a table form
fprintf(head,lab);
fprintf('\n');
disp(Results);

maxErr =  max(C_Err(nc-1));

fprintf('The equations were solved in %1.0f iterations \n ', nc);


%-----------------------------------------------------------
%Functions used in the program.

function [X] = Gauss_Seidel(X,A,B)
%Represents the Gauss-Seidel formula in function form
 global C_Err ES nmax  nc h
    
 nc = 2; %nc: iteration count
 
 while nc < nmax 
    
   if   max(C_Err(nc-1)) < ES   
     break    
   end
     
     for i = 1:4    
     [old_vals ,new_vals] = Edge_Detect(X,A,i,h,nc); %See relevant function
     X(nc,i) = ( B(i) - (old_vals + new_vals)  )/( A(i,i) ); %Gauss-Seidel formula   
     end
     
     for i = 1:4
     C_Err(nc,i) = abs((abs(X(nc,i)) -abs(X(nc-1,i)))/abs(X(nc,i))); %Calculating error
     end    
   nc = nc+1;  
 end    
   
 
end

function [old_vals,new_vals] = Edge_Detect(X,A,i,h,nc)
  %The general equation for Gauss-Seidel (shown in document) has
  %conflictions with MATLAB.In particular regarding which X values use new
  %current values, and which use the values from the previous iterations.
  %When not calculating the first and last X values, the equation holds
  %true. However when at the first and last X values i.e X1 and X4 in a 4x4
  %matrix, the summations go in the range of 0 to 1, and 4 to 5 j terms.
  %When doing this by hand, a person can recognize to skip these terms, but
  %a case statement is needed in MATLAB to spot this. (No shorter way to
  %explain this!)
   switch i
    case h
     old_vals = 0;
     new_vals = A(i,1:i-1)*X(nc,1:i-1)';   
    case 1
     old_vals = A(i,i+1:h)*X(nc-1,i+1:h)';
     new_vals = 0;
    otherwise
     old_vals = A(i,i+1:h)*X(nc-1,i+1:h)';
     new_vals = A(i,1:i-1)*X(nc,1:i-1)'; 
   end
    
end

function [bin] = CheckDiagB(in)
%Checks if a matrix is diagonally dominant, returns a value of 0 or 1.
%Returns 1 if it is diagonally dominant, 0 if it isnt.
i = 1;
n = 1;

while i == 1
    if n > size(in,1)
        break
    end
    
    if abs(in(n,n)) > sum(abs(in(n,n+1:size(in,1)))) %Checks that the value
        i =1;  %is position a(n,n) is larger than the other values in the row.
        
    else
        i = 0;
    end     
    n = n+1;
end

bin = i;

end

function [in,B_in] = MakeDiag(in,B_in)

%This function will make the input matrix diagonally dominant. 
%The input matrix is considered as A(n,m), where m = n.
 
i = 1; % i : iteration count
h = size(in,1);
pass = 0; %pass: Does the matrix pass the diag dom test? 1 = yes, 0 = no

 while pass == 0

   if i == h*(2^h) %Allows for all combinations to be found.
     error('Matrix cannot be made diagonally dominant. Only diagonally dominant matrices can be solved using the Gauss-Seidel method. Dobule check that you have entered the matrix correctly.');  
   end  %Stops the entire program if the matrix cannot be made diagonally dominant
                                       
    
    for n = 1:h-1 %For loop allows function to work for any size of square matrix.
         maxn = max(abs(in(n,:)));
    
     if find(abs(in(n,:))==maxn,2) ~= n %Checks that the max value of row n is in a position of n,n
             

       hold = in(n,:);
       in(n,:) = in(n+1,:); %Row swapping in matrix A
       in(n+1,:) = hold;
       
       B_hold = B_in(n);
       B_in(n) = B_in(n+1); %Swapping the corresponding B rows
       B_in(n+1) = B_hold;
       
       clear hold maxn   %Avoiding potential row duplication
     end
    
     pass = CheckDiagB(in);    
    
    end        
    
    i = i+1; %Increasing the iteration count
       
 end    
  disp('The matrix has been made diagonally dominant');


end

function [bin] = CheckSquare(A_in,B_in)
%Simple function to check that the matrix of co-efficiencts is square.
bin = 0;
 if size(A_in,1) == size(A_in,2)
    if size(A_in,1) == size(B_in,1)
       bin = 1; 
    end
 end   

end

function [bin] = CheckZeros(A_in)
%Simple function to check for zeros in the diagonal of the input matrix A.
%If there are zeros on the diagonal input, the matrix cannot be re-arranged
%by MakeDiag().
  global h
bin = 1;
 for i = 1:h
    if A_in(i,i) ==0 
      bin = 0;  
    end    
 end

end

function [X] = GetGuess()
%User inputs the initial guessed values for X1, X2.. Adjusts to the size of
%the input matrix, i.e a 6x6 matrix will require 6 initial guesses.
 global h 
 disp('Please enter the initial guesses for ');
 guide = ('X%1.0f : \n');

 for i = 1:h
     fprintf(guide,i);
    X(1,i) = input('');
 end    

end

function [nmax,ES] = GetErrIt()
%User enters values for maximum iterations and maximum error, checks for
%nonon-positive values.
pass = 0;
while pass == 0
 n  = input('Please enter the maximum number of iterations you wish to run: ');
 if n <= 0
    disp('Please enter a positive non-zero value');
 else
    pass = 1;
 end

end

    nmax = n+1;

pass = 0;

while pass ==0

ES = input('Please enter the maximum error in the solutions you wish to have: ');

 if ES <= 0                 
    disp('Please enter a positive non-zero value');
 else
    pass = 1;
 end

end


end

function [f,g] = linear_regression(theta, X,y)

  %
  % Arguments:
  %   theta - A vector containing the parameter values to optimize.
  %   X - The examples stored in a matrix.
  %       X(i,j) is the i'th coordinate of the j'th example.
  %   y - The target value for each example.  y(j) is the target for example j.
  %
  
  m=size(X,2);
  n=size(X,1);

  f=0;
  g=zeros(size(theta));

  %
  % TODO:  Compute the linear regression objective efficiently and store the objective
  %        function value in 'f'.
  %
  % TODO:  Compute the gradient of the objective with respect to theta efficiently and  
  %        store the computed gradient in 'g'.
  
%%% YOUR CODE HERE %%%

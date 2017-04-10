format long
A=mmread('unsym.mtx');
b=transpose(size(A,1):-1:1);
x=A\b

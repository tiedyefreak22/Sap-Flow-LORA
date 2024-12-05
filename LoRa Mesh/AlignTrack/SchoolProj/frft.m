function x2 = frft(X, T, U, a)
% Function to Calculate CONTIONOUS FRACTIONAL FOURIER TRANSFORM
% Inputs
% T: 
% X: Input function
% U: Range of frequency
% a: Rotation angle (can't be 0)
% Outputs
% x2:

p=0;
q=0;
x1 = [];
for t=-T:0.01:T
    p=p+1;
    for u=-U:0.1:U
        q=q+1;
%         x1(p,q)=(sqrt((1-j*cot(a))/(2*pi)))*exp(j*(((0.5*(u^2)+0.5*(t^2))*cot(a))-u*t*(1/sin(a))));
        x1(p,q)=(sqrt(1-1i*cot(a*pi/2)))*exp(1i*pi*((u^2)*cot(a*pi/2)+(t^2)*cot(a*pi/2)-2*u*t*(1/sin(a*pi/2))));
    end
    q=0;
end

x2=0.01*(X*x1);

% x3=20*log10(abs(x2)./max(x2));
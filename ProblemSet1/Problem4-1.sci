function y = f(x)
y = ((log(x)-(1-(1/x)))/(log(x)+((1-(1/x))/((5/3)-1))))-0.3
endfunction

function y = fprime(x)
num = -(5/(2*x*x)) - (5*log(x)/(2*x*x)) + (5/(2*x));
den = log(x)*log(x) +3*log(x)+ (9/(4*x*x)) - (((12*log(x))+18)/(4*x)) + (9/4);
y = num/den;
endfunction

//function y = fprime2(x)
//low = log(x) - (3/(2*x)) + (3/2);
//high = (29/(20*x)) + ((7*log(x))/10) - (29/20);
//dlow = (1/x) + (3/(2*x*x));
//dhigh = (-29/(20*x*x)) + (7/(10*x));
//y = ((low*dhigh) - (high*dlow))/(low*low);
//endfunction

xold = 6; // Initial x value

//NEWTON'S METHOD
while abs((f(xold)))>=(10^(-6))
    xold
    f(xold)
    xnew = xold - (f(xold)/fprime(xold))
    xold = xnew;
end

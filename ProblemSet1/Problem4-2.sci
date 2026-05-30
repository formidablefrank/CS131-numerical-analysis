function y = f(x)
   y = cosh(x)*cos(x) + 1
endfunction

//implementation of the secant method with 2 initial points as parameters
function y = secantmethod(xn1,xn2)
    while abs(f(xn1))>(10^(-6))
        xn = xn1 - ((f(xn1)*(xn1-xn2))/(f(xn1)-f(xn2)));
        xn2 = xn1;
        xn1 = xn;
    end
    y = xn;
endfunction

function y = fsubi(betai)
    y = sqrt(((betai*betai*betai*betai*200*0.025*0.0025*0.0025*0.0025)/12)/(4*%pi*%pi*7850*0.9*0.025*0.0025*0.9*0.9*0.9))
endfunction

//(xn-1 and xn-2)
//solving for beta using the function secantmethod
beta1 = secantmethod(-5,-5.5) //leftmost root
beta2 = secantmethod(-2,-2.5) //2nd from the left root
beta3 = secantmethod(1.5,1) //3rd from the left root
beta4 = secantmethod(4.5,4) //rightmost root

//solving for the f sub i's 
fsub1 = fsubi(beta1) 
fsub2 = fsubi(beta2)
fsub3 = fsubi(beta3)
fsub4 = fsubi(beta4)
fsub5 = fsubi(1)

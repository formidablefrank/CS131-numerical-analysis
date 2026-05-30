function[F]=X(t,x,y,z,a)
    F=(-4*(%pi^2))log(((x^2)+(y^2)+(z^2))^(1/2))+a)
endfunction

function[Yi]=RK4(f,t,y,h)
    k1=f(t,y)*h;
    k2=f(t+(h/2),y+(k1/2))*h;
    k3=f(t+(h/2),y+(k2/2))*h;
    k4=f(t+h,y+k3)*h;
    Yi=y+((1/6)*(k1+(2*k2)+(2*k3)+k4));
endfunction

function[]=adapt(f,t,y,h)
    stop=0;
    p=1;//first-order accuracy
    tol =10^(-6);
    emin = tol/(2^(p+1));
    emax = tol;
    while stop==0
        Y1=RK4(f,t,y,h);
        Y2=RK4(f,t,y,h/2);
        Y2=RK4(f,t,Y2,h/2);
        E=abs((Y2-Y1)/((2^p)-1));
        if emin>E then
            Y=Y2;
            h=2*h;
        elseif emax<E then h=h/2
        else 
            Y=Y2;
            stop=1;
        end
    
    end
endfunction

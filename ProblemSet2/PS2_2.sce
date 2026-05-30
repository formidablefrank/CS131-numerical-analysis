function [F] = f(x)
   F=exp(-abs(x));
   //mprintf("%f\n",F); 
endfunction

function [R]=rom(p,f)
    F=0;
    tol=10^(-6)
    h=p*2;
    R(1,1)=f(h);
    h=h/2;
    R(1,2)=f(h);
    R(2,2)=(R(1,1) - (4*R(1,2)))/(-3);
    i=2;
    j=2;
    while (abs(R(i,j)-R(i-1,j-1))>tol)
        h=h/2;
        i=0;
        j=j+1;
        while i<j
            i=i+1;
            if i==1 then R(i,j)=f(h);
            else
                R(i,j)=(R(i-1,j-1) - ((4^(i-1))*R(i-1,j)))/(1-4^(i-1));
            end
            //mprintf("%d %d\n",i,j);
        end
    end

    R=R(i,j);

endfunction

function[O]=Ao(p,rom,f)
    O=(1/p)*rom(p,f);
endfunction

function [R]=romA(p,f,i)
    F=0;
    tol=10^(-6)
    h=p*2;
    R(1,1)=f(h)*cos((i*%pi*h)/p);
    h=h/2;
    R(1,2)=f(h)*cos((i*%pi*h)/p);
    R(2,2)=(R(1,1) - (4*R(1,2)))/(-3);
    i=2;
    j=2;
    while (abs(R(i,j)-R(i-1,j-1))>tol)
        h=h/2;
        i=0;
        j=j+1;
        while i<j
            i=i+1;
            if i==1 then R(i,j)=f(h)*cos((i*%pi*h)/p);
            else
                R(i,j)=(R(i-1,j-1) - ((4^(i-1))*R(i-1,j)))/(1-4^(i-1));
            end
            //mprintf("%d %d\n",i,j);
        end
    end

    R=R(i,j);

endfunction

function[I]=Ai(p,romA,f,i)
    I=(1/p)*romA(p,f,i);
endfunction

function [R]=romB(p,f,i)
    F=0;
    tol=10^(-6)
    h=p*2;
    R(1,1)=f(h)*sin((i*%pi*h)/p);
    h=h/2;
    R(1,2)=f(h)*sin((i*%pi*h)/p);
    R(2,2)=(R(1,1) - (4*R(1,2)))/(-3);
    i=2;
    j=2;
    while (abs(R(i,j)-R(i-1,j-1))>tol)
        h=h/2;
        i=0;
        j=j+1;
        while i<j
            i=i+1;
            if i==1 then R(i,j)=f(h)*sin((i*%pi*h)/p);
            else
                R(i,j)=(R(i-1,j-1) - ((4^(i-1))*R(i-1,j)))/(1-4^(i-1));
            end
            //mprintf("%d %d\n",i,j);
        end
    end

    R=R(i,j);

endfunction

function[I]=Bi(p,romB,f,i)
    I=(1/p)*romB(p,f,i);
endfunction


for n=1:100
    A=0;
    B=0;
    p=1;
    x=1;
    for i=1:n
        A=A+(Ai(p,romA,f,i)*cos((i*%pi*x)/p));
        B=B+(Bi(p,romB,f,i)*sin((i*%pi*x)/p));
    end
    O=Ao(p,rom,f)/p;
    F=O+A+B;
    mprintf("%d %f\n",n,F);
end

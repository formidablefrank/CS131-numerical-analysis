funcprot(0)

//Gaussian function
function y = p1(x)
    y = 1/sqrt(6*%pi)*%e^(-x^2/6)
endfunction

//t-distribution function
function y = p2(x)
    y = 1/(0.886226925453*sqrt(3*%pi))*(1+x^2/3)^(-2)
endfunction

//definite integral of Gaussian function for z >= 0
//Simpson's 1/3 rule with 5 panels, 6 points
function y = P1(x)
    h = x/5
    y = h/3*p1(0) + 4*h/3*p1(x/5) + 2*h/3*p1(2*x/5) + 4*h/3*p1(3*x/5) + 2*h/3*p1(4*x/5) + h/3*p1(x)
endfunction

//definite integral of t-distribution function for z >= 0
//Simpson's 1/3 rule with 5 panels, 6 points
function y = P2(x)
    h = x/5
    y = h/3*p2(0) + 4*h/3*p2(x/5) + 2*h/3*p2(2*x/5) + 4*h/3*p2(3*x/5) + 2*h/3*p2(4*x/5) + h/3*p2(x)
endfunction

//root-finding function for Gaussian function
function y = Q1(z, a)
    y = 0.5 + P1(z) - a
endfunction

//root-finding function for t-distribution function
function y = Q2(z, a)
    y = 0.5 + P2(z) - a
endfunction

//regula-falsi method for Gaussian function
//value at risk at alpha
function m = RF1(alpha)
    a = 1
    b = 40
    if(Q1(a, alpha)*Q1(b, alpha) > 0) then
        error('redefine bracket')
    else
        while 1
            m = (a * Q1(b, alpha) - b * Q1(a, alpha)) / (Q1(b, alpha) - Q1(a, alpha))
            if (Q1(m, alpha) < 1D-9) then
                return
                disp('nice')
            else
                if(Q1(a, alpha) * Q1(m, alpha) > 0) then
                    a = m
                else
                    b = m
                end
            end
        end
    end
endfunction

//regula-falsi method for t distribution function
//value at risk at alpha
function m = RF2(alpha)
    a = 1
    b = 40
    if(Q2(a, alpha)*Q2(b, alpha) > 0) then
        error('redefine bracket')
    else
        while 1
            m = (a * Q2(b, alpha) - b * Q2(a, alpha)) / (Q2(b, alpha) - Q2(a, alpha))
            if (Q2(m, alpha) < 1D-9) then
                return
            else
                if(Q2(a, alpha) * Q2(m, alpha) > 0) then
                    a = m
                else
                    b = m
                end
            end
        end
    end
endfunction

function myGaussianDriver()
    x = []
    y = []
    for val=0.8:0.01:0.99
        x = [x val]
        y = [y RF1(val)]
    end
    disp('Gaussian:')
    disp('  alpha  value-at-risk')
    disp([x' y'])
    plot(x, y, 'r')
endfunction

function myTdistDriver()
    x = []
    y = []
    for val=0.8:0.01:0.99
        x = [x val]
        y = [y RF2(val)]
    end
    disp('T-dist')
    disp('  alpha  value-at-risk')
    disp([x' y'])
    plot(x, y, 'g')
    hl = legend(['Gaussian','t distribution'])
    xtitle('Value-at-Risk', 'Level of confidence (alpha)', 'Value-at-Risk (z)')
endfunction

myGaussianDriver()
myTdistDriver()

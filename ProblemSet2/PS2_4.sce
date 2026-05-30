funcprot(0)
//Given a matrix M, it returns its reduced row-echelon form.
function X = myGJR(M)
    [m n] = size(M)
    for k = 1:m
        pivot = argmax(M, k)
        M = swap(M, k, pivot)
        if abs(M(k, k)) <= 1D-10 then error('singular matrix!'); end
        M(k, k:$) = M(k, k:$) / M(k, k)
        for i = k+1:m
            for j = n:-1:k
                M(i, j) = M(i, j) - M(k, j) * (M(i, k) / M(k, k))
            end
            M(i, k) = 0
        end
    end
    for k = m:-1:1
        for i = k-1:-1:1
            for j = n:-1:k
                if M(k, k) ~= 0 then
                    M(i, j) = M(i, j) - M(k, j) * (M(i, k) / M(k, k))
                end
            end
        end
    end
    X = M
endfunction

//Swap rows a and b of a matrix A.
function X = swap(A, a, b)
    temp = A(a, :)
    A(a, :) = A(b, :)
    A(b, :) = temp
    X = A
endfunction

//Returns a index number of a row that has the maximum absolute value of an element
//For a matrix M in column k
function X = argmax(M, k)
    [m n] = size(M)
    X = k
    temp = 0
    for i=k:m
        if abs(M(i, k)) > temp then
            temp = abs(M(i, k))
            X = i
        end
    end
endfunction

//logistic model
function y = l(t, x, P)
    y = x(1)*P*(x(2)-P)
endfunction

//RK4 solution of logistic model
function y = RK4L(t, x, P, dt)
    x1 = l(t, x, P)
    x2 = l(t + dt/2, x, P + x1*dt/2)
    x3 = l(t + dt/2, x, P + x2*dt/2)
    x4 = l(t + dt, x, P + x3*dt)
    X = (x1 + 2*x2 + 2*x3 + x4)/6
    y = P + X*dt
endfunction

//gompertz model
function y = g(t, x, P)
    y = x(1)*P*log(x(2)-P)
endfunction

//RK4 solution of Gompertz model
function y = RK4G(t, x, P, dt)
    x1 = g(t, x, P)
    x2 = g(t + dt/2, x, P + x1*dt/2)
    x3 = g(t + dt/2, x, P + x2*dt/2)
    x4 = g(t + dt, x, P + x3*dt)
    X = (x1 + 2*x2 + 2*x3 + x4)/6
    y = P + X*dt
endfunction

//x is the parameter vector that contains parameters k and C, P is the independent variable. this is the function we want to optimize.
//x(1)=k; x(2)=C
function y = fL(x, P)
    y = P + (x(1)/6)*P*(x(2)-P) + 2*x(1)/6*(x(2)-P-(x(1)/2)*P*(x(2)-P))*(P+(x(1)/2)*P*(x(2)-P)) + (2*x(1)/6)*(x(2)-P-(x(1)/2)*(x(2)-P-(x(1)/2)*P*(x(2)-P))*(P+(x(2)/2)*P*(x(2)-P)))*(P+(x(1)/2)*(x(2)-P-(x(1)/2)*P*(x(2)-P))*(P+(x(2)/2)*P*(x(2)-P))) + (x(1)/6)*(x(2)-P-x(1)*(x(2)-P-(x(1)/2)*(x(2)-P-(x(1)/2)*P*(x(2)-P))*(P+x(1)/2*P*(x(2)-P)))*(P+(x(1)/2)*(x(2)-P-(x(1)/2)*P*(x(2)-P))*(P+x(1)/2*P*(x(2)-P))))*(P+x(1)*(x(2)-P-(x(1)/2)*(x(2)-P-(x(1)/2)*P*(x(2)-P))*(P+x(1)/2*P*(x(2)-P)))*(P+(x(1)/2)*(x(2)-P-(x(1)/2)*P*(x(2)-P))*(P+x(1)/2*P*(x(2)-P)))) 
endfunction

function y = fG(x, P)
    y = P + (x(1)/6)*P*log(x(2)-P) + 2*x(1)/6*log(x(2)-P-(x(1)/2)*P*log(x(2)-P))*(P+(x(1)/2)*P*log(x(2)-P)) + (2*x(1)/6)*log(x(2)-P-(x(1)/2)*log(x(2)-P-(x(1)/2)*P*log(x(2)-P))*(P+(x(2)/2)*P*log(x(2)-P)))*(P+(x(1)/2)*log(x(2)-P-(x(1)/2)*P*log(x(2)-P))*(P+(x(2)/2)*P*log(x(2)-P))) + (x(1)/6)*log(x(2)-P-x(1)*log(x(2)-P-(x(1)/2)*log(x(2)-P-(x(1)/2)*P*log(x(2)-P))*(P+x(1)/2*P*log(x(2)-P)))*(P+(x(1)/2)*log(x(2)-P-(x(1)/2)*P*log(x(2)-P))*(P+x(1)/2*P*log(x(2)-P))))*(P+x(1)*log(x(2)-P-(x(1)/2)*log(x(2)-P-(x(1)/2)*P*log(x(2)-P))*(P+x(1)/2*P*log(x(2)-P)))*(P+(x(1)/2)*log(x(2)-P-(x(1)/2)*P*log(x(2)-P))*(P+x(1)/2*P*log(x(2)-P)))) 
endfunction

//population
P = [1.03519 1.06937 1.12796 1.17191 1.19633 1.21586 1.23051 1.25004 1.26469 1.28422 1.29398 1.30375 1.30863 1.31352 1.32328 1.32328 1.32328 1.32817 1.32817 1.33305 1.33305]

function myLogisticDriver()
    //P = previous population, x = column vector of parameters (x_old)
    x_new = [1 1]'
    disp('Running Gauss-Newton optimization for logistic model...')
    while 1
        J = []
        r = []
        x = x_new
        for val = [1:20]
            r(val,:) = fL(x, P(val)) - P(val+1)
            J(val,:) = derivative(list(fL, P(val)), x, order = 4)
        end
        //J is the Jacobian matrix of residual matrix r.
        temp = myGJR([J'*J J'*r])
        delta = temp(:, $)
        x_new = x - delta
        if(abs(norm(x_new - x))< 1D-9) then
            disp('optimal parameter <k, C>:')
            disp(x_new)
            break
        end
    end
    s = [0:0.05:20]
    M = [P(1)]
    for val = [1:400]
        M = [M RK4L(s(val), x_new, M(val), 0.05)]
    end
    disp('   time     population')
    disp([s' M'])
endfunction

function myGompertzDriver()
    //P = previous population, x = column vector of parameters (x_old)
    x_new = [1 1]'
    disp('Running Gauss-Newton optimization for Gompertz model...')
    while 1
        J = []
        r = []
        x = x_new
        for val = [1:20]
            r(val,:) = fG(x, P(val)) - P(val+1)
            J(val,:) = derivative(list(fG, P(val)), x, order = 4)
        end
        //J is the Jacobian matrix of residual matrix r.
        temp = myGJR([J'*J J'*r])
        delta = temp(:, $)
        x_new = x - delta
        if(abs(norm(x_new - x))< 1D-9) then
            disp('optimal parameter <k, C>:')
            //x_new = real(x_new)
            disp(x_new)
            break
        end
    end
    s = [0:0.05:20]
    M = [P(1)]
    for val = [1:400]
        M = [M RK4G(s(val), x_new, M(val), 0.05)]
    end
    disp('   time     population')
    disp([s' M'])
endfunction

//myLogisticDriver()
myGompertzDriver()

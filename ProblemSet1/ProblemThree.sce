funcprot(0)

function M = myDivDiff(X, Y)
    [m n] = size(X)
    M = zeros(n, n)
    M(:, 1) = Y'
    for col = 2:n
        for row = 2:n
            if col > row then
                continue
            else
                M(row, col) = (M(row, col-1) - M(col-1, col-1)) / (X(row)-X(col-1))
            end
        end
    end
endfunction

function Z = myNewtonInterpolate(M, X, val)
    [m n] = size(X)
    for i = 0:n-1
        if i == 0 then
            y = M(n, n)
        else
            k = n-i
            y = M(k, k) + (val - X(1, k)) * temp
        end
        temp = y
    end
    Z = temp
endfunction

function myNewtonDriver(X, Y)
    disp('Newton Interpolation:')
    M = myDivDiff(X, Y)
    disp('Divided differences matrix:')
    disp(M)
    newX = []
    newY = []
    for val = X(1, 1):10:X(1, $)
        newX = [newX val]
        newY = [newY myNewtonInterpolate(M, X, val)]
    end
    disp('Temp (K)   Specific Heat (cp) [Newton interpolation at intervals of 10K]')
    disp([newX' newY'])
    disp('Generating plot...')
    scf(0)
    //trueX = [77 173 273 373 573 773]
    //trueY = [0.336 0.743 0.880 0.937 1.021 1.130]
    plot(newX, newY, 'r')
    //plot(trueX, trueY)
    //hl = legend(['Interpolated values','True values'])
    //xtitle( 'Newton Interpolation for the specific heat of aluminum', 'Temperature (K)', 'Specific Heat (J/gk)')
    Xnot = myNewtonInterpolate(M, X, 273)
    disp(msprintf('At T = 273 K, the interpolated specific heat is %f J/gK while the true specific heat is 0.880 J/gK', Xnot))
    disp(msprintf('Absolute error: %f, relative error: %f', abs(Xnot-0.88), Xnot/0.88-1))
    disp(msprintf('Using regula falsi from the interval [150, 300] with tol 1 x 10^-5, 1.0 J/gk is in %f K', myRegulaFalsiNewton(150, 300, 1D-5, M, X)))
    disp(msprintf('--------------------------------------------------------\n'))
endfunction

function R = myRegulaFalsiNewton(a, b, tol, M, X)
    if (myNewtonInterpolate(M, X, a) - 1) * (myNewtonInterpolate(M, X, b) - 1) > 0 then
        error('Define a new interval')
    else
        while 1
            R = (a * (myNewtonInterpolate(M, X, b) - 1) - b * (myNewtonInterpolate(M, X, a) - 1)) / ((myNewtonInterpolate(M, X, b) - 1) - (myNewtonInterpolate(M, X, a) - 1))
            if (myNewtonInterpolate(M, X, R) - 1) < tol then
                return
            else
                if (myNewtonInterpolate(M, X, a) - 1) * (myNewtonInterpolate(M, X, R) - 1) < 0 then
                    b = R
                else
                    a = R
                end
            end
        end
    end
endfunction

function K = myCubicSplineK(X, Y)
    [m n] = size(X)
    A = zeros(n-2, n-2)
    b = []
    for i = 2:n-1
        b(1, i-1) = 6*((Y(i+1)-Y(i))/(X(i+1)-X(i))-(Y(i)-Y(i-1))/(X(i)-X(i-1)))
        r = X(1, i) - X(1, i-1)
        s = 2*(X(1, i+1) - X(1, i-1))
        t = X(1, i+1) - X(1, i)
        if i == 2 then
            A(1, 1) = s
            A(1, 2) = t
        elseif i == n-1
            A(n-2, n-3) = r
            A(n-2, n-2) = s
        else
            A(i-1, i-2) = r
            A(i-1, i-1) = s
            A(i-1, i) = t
        end
    end
    disp('Tridiagonal matrix:')
    disp(A)
    disp('Coefficient matrix:')
    disp(b)
    K = myGJR([A b'])
    K = [0 K(:, $)' 0]
endfunction

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

function X = swap(A, a, b)
    temp = A(a, :)
    A(a, :) = A(b, :)
    A(b, :) = temp
    X = A
endfunction

function Z = myCubicInterpolate(K, X, Y, val)
    [m n] = size(X)
    v = 0
    for i = 1:n-1
        if val >= X(1, i) then
            if val < X(1, i+1) then v = i; end
        end
    end
    v = 5
    Z = K(1, v)*(X(1, v+1)-val)^3 / 6 / (X(1, v+1)-X(1, v)) + K(1, v+1)*(val-X(1, v))^3/6/ (X(1, v+1)-X(1, v))+ (X(1, v+1)-val)*(Y(1, v)/(X(1, v+1)-X(1, v))-K(1, v)*(X(1, v+1)-X(1, v))/6) + (val-X(1, v))*(Y(1, v+1)/(X(1, v+1)-X(1, v)) - K(1, v+1)/6*(X(1, v+1)-X(1, v)))
endfunction

function R = myRegulaFalsiCubic(a, b, tol, K, X, Y)
    if (myCubicInterpolate(K, X, Y, a) - 1) * (myCubicInterpolate(K, X, Y, b) - 1) > 0 then
        error('Define a new interval')
    else
        while 1
            R = (a * (myCubicInterpolate(K, X, Y, b) - 1) - b * (myCubicInterpolate(K, X, Y, a) - 1)) / ((myCubicInterpolate(K, X, Y, b) - 1) - (myCubicInterpolate(K, X, Y, a) - 1))
            if (myCubicInterpolate(K, X, Y, R) - 1) < tol then
                return
            else
                if (myCubicInterpolate(K, X, Y, a) - 1) * (myCubicInterpolate(K, X, Y, R) - 1) < 0 then
                    b = R
                else
                    a = R
                end
            end
        end
    end
endfunction

function myCubicSplineDriver(X, Y)
    disp('Cubic Spline Interpolation:')
    K = myCubicSplineK(X, Y)
    disp('k:')
    disp(K)
    newX = []
    newY = []
    for val = X(1, 1):10:X(1, $)
        newX = [newX val]
        newY = [newY myCubicInterpolate(K, X, Y, val)]
    end
    disp('Temp (K)   Specific Heat (cp) [Cubic spline interpolation at intervals of 10K]')
    disp([newX' newY'])
    disp('Generating plot...')
    //scf(1)
    plot(newX, newY, 'g')
    trueX = [77 173 273 373 573 773]
    trueY = [0.336 0.743 0.880 0.937 1.021 1.130]
    plot(trueX, trueY)
    hl = legend(['Newton','Cubic Spline','True values'])
    xtitle('Interpolation for the specific heat of aluminum', 'Temperature (K)', 'Specific Heat (J/gk)')
    Xnot = myCubicInterpolate(K, X, Y, 273)
    disp(msprintf('At T = 273 K, the interpolated specific heat is %f J/gK while the true specific heat is 0.880 J/gK', Xnot))
    disp(msprintf('Absolute error: %f, relative error: %f', abs(Xnot-0.88), Xnot/0.88-1))
    disp(msprintf('Using regula falsi from the interval [150, 300] with tol 1 x 10^-5, 1.0 J/gk is in %f K', myRegulaFalsiCubic(150, 300, 1D-5, K, X, Y)))
    disp(msprintf('--------------------------------------------------------\n'))
endfunction

X = [-250 -200 -100 0 100 300]
Y = [0.0163 0.318 0.699 0.870 0.941 1.04]
myNewtonDriver(X, Y)
myCubicSplineDriver(X, Y)

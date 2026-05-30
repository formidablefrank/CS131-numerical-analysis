funcprot(0)

//Gauss-Jordan Reduction method for matrix inverse.
//Returns the inverse of the matrix
function X = myGJRdriver(M)
    [m n] = size(M)
    A = myGJR([M eye(m,n)])
    X = A(:,n+1:$)
endfunction

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

//Returns the upper triangular form of matrix
function S = uppert(M)
    S = M
    [m n] = size(M)
    for i=1:m
        for j=1:n
            if j < i 
               S(i,j)=0
            end
        end
    end
endfunction

//Returns the lower triangular form of matrix
function T = lowert(M)
    T = M
    [m n] = size(M)
    for i=1:m
        for j=1:n
            if j >= i
               T(i,j)=0
            end
        end
    end
endfunction

//Returns the lower and upper matrices obtained using Crout's decomposition
function [L, U]= myLU(M)
    [y z]=size(M)
    for i=1:z-1
        M(i+1:z,i) = M(i+1:z,i) / M(i,i)
        M(i+1:z,i+1:z) = M(i+1:z,i+1:z) - M(i+1:z,i) * M(i,i+1:z)
    end
    
    L=eye(z,z)+lowert(M)
    U=uppert(M)
endfunction

function y = ForwardSub(L,b)
    [m n] = size(L)
    y = zeros(m, 1)
    y(1) = b(1) / L(1, 1)
    for i = 2:m
        y(i) = b(i) - L(i, 1:i-1) * y(1:i-1)
        y(i) = y(i) / L(i, i)
    end
endfunction

function x = BackSub(U, y)
    [m n] = size(U)
    x = zeros(m, 1)
    x(m) = y(m) / U(m,m)
    for i = m-1:-1:1
        x(i) = y(i) - U(i, i+1:m) * x(i+1:m)
        x(i) = x(i) / U(i, i)
    end
endfunction

//Given matrix M, returns its inverse using LU Factorization
function A = myLUdriver(M)
    [L U] = myLU(M)
    [m n] = size(M)
    I = eye(m, n)
    B = []
    X = []
    for i = 1:n
        B = [B ForwardSub(L, I(:, i))]
        X = [X BackSub(U, B(:, i))]
    end
    A = X
endfunction

function X = myOneNorm(A)
    [m n] = size(A)
    maxi = 0
    for i = 1:m
        temp = 0
        for j = 1:n
            temp = temp + abs(A(j, i))
        end
        if temp > maxi then
            maxi = temp
        end
    end
    X = maxi
endfunction

function X = myInfNorm(A)
    [m n] = size(A)
    maxi = 0
    for i = 1:m
        temp = 0
        for j = 1:n
            temp = temp + abs(A(i, j))
        end
        if temp > maxi then
            maxi = temp
        end
    end
    X = maxi
endfunction

//Iterative method for matrix inverse, initial guess by Heath
function X = myIter(M)
    [m n] = size(M)
    tol = 1D-3
    X_old = M'./(myOneNorm(M) * myInfNorm(M) + 1)
    i = 1
    while 1
        X_new = X_old*(2*eye(m, n) - M*X_old)
        if myOneNorm(X_new - X_old) < tol then
            disp(sprintf("Took %d iterations", i))
            break
        end
        X_old = X_new
        i = i + 1
    end
    X = X_new
endfunction

//Iterative method for matrix inverse, initial guess from the AMS Journal
function X = myIter2(M)
    [m n] = size(M)
    tol = 1D-3
    X_old = M'./max(real(spec(M*M')))
    i = 1
    while 1
        X_new = X_old*(2*eye(m, n) - M*X_old)
        if myOneNorm(X_new - X_old) < tol then
            disp(sprintf("Took %d iterations", i))
            break
        end
        X_old = X_new
        i = i + 1
    end
    X = X_new
endfunction

//Given a matrix X, returns its reduced form with its row r and column c removed.
function A = reduce(X, r, c)
    [m n] = size(X)
    Y = zeros(m-1, n-1)
    j = 1
    for h = 1:n
        if h == r then
            continue
        end
        k = 1
        for i = 1:n
            if i == c then
                continue
            end
            Y(j, k) = X(h, i)
            k = k + 1
        end
        j = j + 1
    end
    A = Y
endfunction

//Computes the determinant of a matrix: base case is when matrix is 3x3.
//Use only for small matrices, else, recursive calls are deep
function A = myDet(X)
    [m n] = size(X)
    if m ~= n then
        error('not a square matrix!')
    end
    maxi = m
    ret = 0
    if maxi == 3 then
        A = X(1, 1) * X(2, 2) * X(3, 3) + X(1, 2) * X(2, 3) * X(3, 1) + X(1, 3) * X(2, 1) * X(3, 2) - X(1, 3) * X(2, 2) * X(3, 1) - X(1, 2) * X(2, 1) * X(3, 3) - X(1, 1) * X(2, 3) * X(3, 2)
    else
        Y = []
        for h = 1:maxi
            if X(h, 1) == 0 then
                continue
            end
            Y = reduce(X, h, 1)
            if modulo(h, 2) == 0 then
                ret = ret - myDet(Y) * X(h, 1)
            else
                ret = ret + myDet(Y) * X(h, 1)
            end
        end
        A = ret
    end
endfunction

//Given a matrix X, returns if it is nice
function B = nice(X)
    B = 1
    try
        M = myGJR(X)
        [M N] = myLU(X)
    catch
        B = 0
    end
endfunction

function testMethod(method, M)
    disp(msprintf('Inverse using %s:', method))
    tic()
    if method == 'GJRM' then
        M_inv = myGJRdriver(M)
    elseif method == 'LU' then
        M_inv = myLUdriver(M)
    elseif method == 'Iter' then
        M_inv = myIter(M)
    end
    time = toc()
    disp(M_inv)
    disp(msprintf('Time elapsed: %f seconds', time))
endfunction

M = [0]
while nice(M) == 0
    order = grand(1, 1, "uin", 5, 20)
    M = grand(order, order, "uin", -10, 10)
end
disp('Random matrix:')
disp(M)
testMethod('GJRM', M)
testMethod('LU', M)
testMethod('Iter', M)

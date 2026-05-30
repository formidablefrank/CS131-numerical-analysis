function A = GE(M)
    A = M
    [m n] = size(A)
    for col = 1:n
        A(col,col:n) = A(col,col:n)/A(col,col)
        A(2,1:$) = A(2,1:$)-A(2,1)*A(1,1:$)
//        for row = col+1:n
//            disp(row)
//            disp(col)
//            A(row,col:n) = A(row,col:n) - A(row:col)*A(col,col:n)
//        end
    end
endfunction

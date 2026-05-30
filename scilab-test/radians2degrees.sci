function [degrees] = radians2degrees(radians)
    degrees = radians * (180 / %pi)
endfunction

deff('[radians] = degrees2radians(degrees)','radians = degrees * %pi / 180')

function [fact] = factorial(n)
    if n < 2 then
        fact = 1
    else
        fact = n * factorial(n-1)
    end
endfunction

function [ans] = sample(x)
    r = 2.220D-16
    ans = log10(r/(10^x) + 10^x)
endfunction

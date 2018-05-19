% f = input signal (1xa)
% h = impulse response (1xb)
% g = output signal (1x(a + b-1))
function [g] = convm1d(f,h)
    a = length(h);
    b = length(f);
    n = a + b-1;
    he = [h zeros(1,n-a)];
    fe = [f zeros(1,n-b)]';
    h1 = [he(1) he(n:-1:2)];
    H = h1;
    for i = 1:n-1
        h1 = [h1(end) h1(1:end-1)];
        H = [H;h1];
    end
    g = H*fe;
    g = g'
end
function outbuf = hermnu(inptr, points)

%points are (non-integer) index values into inptr that are to be
%interpolated

outsamps=0;
outbuf = zeros(size(points));
for i=1:length(points)
    pos = floor(points(i)+1);     %shift by 1 to match c zero indexing
    
    x3 = inptr(pos-3);
    x2 = inptr(pos-2);
    x1 = inptr(pos-1);
    x0 = inptr(pos);
    
    frac = points(i)+1 -floor(points(i)+1);

    %/* 4-tap Hermite, using Farrow structure */
    acc0 = (3 * (x2 - x1) + x0 - x3) / 2;
    acc0 = (acc0 * frac);
    acc0 = acc0 + 2 * x1 + x3 - ((5 * x2 + x0) / 2);
    acc0 = (acc0 * frac);
    acc0 = acc0 + (x1 - x3) / 2;
    acc0 = (acc0 * frac);
    acc0 = acc0 + x2;
    
    outsamps = outsamps+1;   
    outbuf(outsamps) = acc0;
end


function convert16to4chan(in, out)

%f = fopen('warpedFrame_16chan.bin');
f = fopen(in);
d = fread(f, 'uint16=>uint16');
fclose(f);
d2 = reshape(d,16, 6656, numel(d)/6656/16);
%f2 = fopen('warpedFrame_4chan.bin', 'w');
f2 = fopen(out, 'w');
d4=d2(1:4,:,:);
fwrite(f2, d4, 'uint16');
fclose(f2)

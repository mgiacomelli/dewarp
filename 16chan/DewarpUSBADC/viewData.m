function [test, d] = viewData(filename, xpix, ypix)

if(nargin < 2)
    xpix = 1024;
    ypix = 1024;
end

f=fopen(filename);
data=fread(f, xpix*ypix*16, 'int16');
fclose(f);
test=reshape(data,xpix, ypix, 16);
%figure;imagesc(test(700:790,700:790,1))
%figure;imagesc(test(1000:1100,1+1000:2:1100,1))
figure;imagesc(test(:,:,1))
axis square

max(max(test(:,:,2)))
min(min(test(:,:,2)))

f=fopen('myDebugData10.bin');
d = fread(f, 'float32');
fclose(f);
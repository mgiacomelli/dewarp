function [taps, fc, barray,barrayh, remap_index, shifted_cutoff, linedata, filtered_data] = loadDewarpBins()
path = '.\';

f=fopen([path 'myDebugData102.bin']);
taps = fread(f, 'int32');
fclose (f);

%cutoff frequencies
f=fopen([path 'myDebugData501.bin']);
fc = fread(f, 'float');
fclose (f);

%the actual filter tap array
f=fopen([path 'myDebugData100.bin']);
barray = fread(f, 'float');
fclose (f);

f=fopen([path 'myDebugData99.bin']);
barrayh = fread(f, 'uint16=>uint16');   %no native fp16 support in fread
barrayh = half.typecast(barrayh);
fclose (f);

f=fopen([path 'myDebugData10.bin']);
remap_index = fread(f, 'float');
fclose (f);

f=fopen([path 'myDebugData103.bin']);
shifted_cutoff = fread(f, 'double');
fclose (f);

f=fopen([path 'myDebugData401.bin']);
linedata = fread(f, 'uint16=>uint16');  
fclose (f);

f=fopen([path 'myDebugData400.bin']);
filtered_data = fread(f, 'uint16=>uint16');   %no native fp16 support in fread
filtered_data = half.typecast(filtered_data);
fclose (f);


%dewarp a raw frame from the alazar using cubic hermite interpolation
%this code tries to mimick the c implementation as much as possible at the expense of performance
%it also avoids using built in matlab functions which would not be available in c
%slowLines:  how many slow axis positions are in the frame
%fastPixels: how many fast axis pixels to reconstruct
%fir:  use FIR filtering
%type:  type of interpoloation, 0: nearest, 1: linear, 2: hermite
function out = dewarpFrame(path, slowLines, fastPixels, channels, fir, type)

if(nargin <4)
    channels = 16;
    fir = 1;
    type = 2;
end

samples = 6700; %samples per res cycle
dataLen = 6656; %samples recorded (may be less if there is ADC rearm time)

fraction_nyquist = 0.7;  %how oversampled the image is at the center of the FOV (1.0 = no oversampling)

%change the taps every four samples?  Needed in the 4 channel case so that
%we can pack 4 samples of 4 channels into a 16 element YMM register
tap_by_four = false;

%calculate the interpolation table
table = zeros(fastPixels*2,1);
m=samples;
n=fastPixels;
%for standard test image
shift = 34;

%shift L by 0.5 so that we index to the center of each bin rather than the left edge 
for L = 0:fastPixels-1      
    table(L+1) = m/2/pi*acos(1-2*(L+0.5)/n) - shift;
end


for L = fastPixels:2*fastPixels-1
    table(3*n-L)=m/2+m/2/pi*acos(2/n*(L+0.5-n)-1)  - shift;
end


%clamp ends for interpolation
table(table<3)=3;
table(table>dataLen-1)=dataLen-1;   %-1 since we are going to shift over to match the zero index c code


%FIR tap calculation
outtaps = zeros(fastPixels,1);


%calculate the cutoff frequency assuming that the image is Nyquist
%sampled at samples/4 (center FOV)
oversample = diff(table(1:fastPixels));
%index 1 and end are missing/wrong due to diff, so duplicate 2 and end-1
oversample(1)=oversample(2);
oversample(end+1) = oversample(end);


%cutoff frequencies, note symetric so can use for both directions
cutoff = 1./oversample;
cutoff = cutoff*fraction_nyquist;
shifted_cutoff = zeros(size(cutoff));

%calculate the number of taps for each pixel assuming the transition
%band will be 3.5/Ntaps wide
for i=1:fastPixels

    %TODO:  should replace this with a smoothed function that changes the
    %taps more slowly!
    ntaps = 3*ceil(oversample(i));

    %calculate the transistion band half width
    tband = 1.75/ntaps;
     shifted_cutoff2(i) = cutoff(i)+tband;
    %make sure the tband isn't too large to be useful
    newtaps = ntaps;

	%0.25 is a good tradeoff since the computational complexity goes up rapidly for very small improvement in pass band
    while(tband> 0.25)

        newtaps = newtaps+1;
        tband = 1.75/newtaps;
    end

%     %optional feature to keep ntaps constant for groups of 4 samples to
%     %enable more efficient vector processing in the 4 channel case
%    if(tap_by_four &&  mod(i,4) ~= 1)
%        newtaps=outtaps(i-1);
%    end
%        

    %first shift the cut off by half the transition band
    newcutoff = cutoff(i)+tband;


   
    shifted_cutoff(i) =newcutoff;
    %now round up oversample to nearest even integer

    %always use even n for integer group delay (fir1 adds 1)
    
    if(mod(newtaps, 2) == 1)
        newtaps = newtaps+1;
    end
    outtaps(i) = newtaps;




end

%clamp cutoff at 1
shifted_cutoff(shifted_cutoff>1) = 1;

%nearest neighbor the taps/cutoff to the sample domain from pixel domain
ii = table(1);
tablepos = 1;

cutoff_samp(1:ceil(table(1))) = shifted_cutoff(1);    %get any samples on the edge
taps_samp(1:ceil(table(1))) = outtaps(1);


%should this be floor?  we double count but i think its ok?
for i=floor(table(1)):samples/2

    if(i>=ii)
        tablepos=tablepos+1;
        ii = table(tablepos);
    end
    if(tablepos>length(shifted_cutoff))
        %points after the last table value aren't used anyway
        cutoff_samp(i) = shifted_cutoff(end);
        taps_samp(i) = outtaps(end);
    else
        cutoff_samp(i) = shifted_cutoff(tablepos);
        %optional feature to update taps in groups of 4 for easier
        %vectorization
        if(tap_by_four && mod(i,4)~=1)
            taps_samp(i)=taps_samp(i-1);
        else
            taps_samp(i) = outtaps(tablepos);
        end
    end
end

%cutoff_samp = smooth(cutoff_samp, 100);
%calculate the filter coefficients


barray = zeros(dataLen/2,max(taps_samp));

%flattened (1D) version of the barray
%barray2 = zeros(dataLen/2,max(taps_samp));
bpos = 1;
for i=1:dataLen/2
   
    if(cutoff_samp(i) > 0.99)     %can't believe fir1 fails at this
        taps_samp(i) = 1;
        barray(i,1) = 1;
        
        barray2(bpos) = 1;
         bpos = bpos+1;
    else
		%can use the built in matlab fir1 function for comparison
        %b = fir1(taps_samp(i),cutoff_samp(i), 'low');
		%windowed sinc is nice because its simple to implement in c
       b = wsfiltgen(taps_samp(i)+2, cutoff_samp(i));
		
		%can zero the taps corresponding to missing data at the edges,
		%but this is a bad idea since it will result in non-constant delay at edges
       if(0)
		    if (i <= taps_samp(i) / 2)
		    
			    %handle the boundary condition	    
                b(1:taps_samp(i) / 2) = 0;
                %renormalize so DC is correct
                b=b./sum(b);
    
           end
       end
      
        barray(i,1:taps_samp(i)+1) = b;
        
        %1D ARRAY AS NEEDED FOR C OPTIMIZATION
        barray2(bpos:bpos+taps_samp(i)) = b;

        bpos = bpos + taps_samp(i)+1;
    end
    
end

%read the data
f=fopen(path);
d=fread(f, 'uint16');
fclose(f);

d2 = reshape(d,channels, dataLen, slowLines);

%write out the raw data for reference
if(0)
    img = squeeze(d2(1:3,:,:));
    img = img - mean(mean(mean(d(3,:,:))));
    img= -1*img;
    test = zeros(dataLen, slowLines, 3);
    
    test(:,:,1) = 4*squeeze(img(1,:,:))./max(max(squeeze(img(1,2:end,2:end))));
    test(:,:,2) = squeeze(img(2,:,:))./max(max(squeeze(img(2,:,:))));
    imwrite(test(:,16:end,:), 'raw.png');
end
%d2filt = zeros(size(d2));
out = zeros(channels, fastPixels, slowLines);

%d2(1,:,1)=1:6656;

%dewarp the data
tic
for i=1:slowLines
    for j=1:channels
    %for j=1:4
        %filter the entire sample data
        if(fir)
            if(0)   %impulse response testing
                d2(j,:,i) = zeros(dataLen,1);
                d2(j,1700,i) = 100000;
            end

			%this is very slow but mimicks the way the c code accesses arrays for testing, could be optimized if anyone wants the matlab code
            if(i>1 && i< slowLines)

                threeline = [squeeze(d2(j,:, i-1)) squeeze(d2(j,:, i)) squeeze(d2(j,:, i+1))];
            elseif(i==1)

                threeline = [zeros(1,dataLen) squeeze(d2(j,:, i)) squeeze(d2(j,:, i+1))];
            else
                threeline = [squeeze(d2(j,:, i-1)) squeeze(d2(j,:, i)) fliplr(squeeze(d2(j,:, i)))];
            end
            filtered = zeros(dataLen,1);
            bindx = 1;
            for k=1:dataLen/2 %sample loop

                ntaps = taps_samp(k);
                if(ntaps>1)


                    for n=1:ntaps+1   %FIR loop 

                        startindex = dataLen+k-1;


                        filtered(k) = filtered(k) + threeline(startindex-ntaps/2+n)*barray2(bindx);
    
                        %barray is symmetric, so can use the same coefficient for
                        %the backwards sweep pixel

                        filtered(dataLen/2+k) = filtered(dataLen/2+k) + threeline(dataLen/2+startindex-ntaps/2+n)*barray2(bindx);

                        bindx=bindx+1;
                    end
                else
                    filtered(dataLen/2+k) = squeeze(d2(j,dataLen/2+k, i));
                    filtered(k) = squeeze(d2(j,k, i));
                    bindx=bindx+1;
                end
                

            end
            
        else
            filtered = squeeze(d2(j,:, i));
        end
    
        switch type
            case 0
                out(j, :, i*2-1) =interp1(filtered,-2+table(1:fastPixels),'nearest');    
                out(j, :, i*2) = flipud (interp1(filtered, (-2+table(fastPixels+1:fastPixels*2)), 'nearest'));
            case 1
                out(j, :, i*2-1) =interp1(filtered,-2+table(1:fastPixels),'linear');    
                out(j, :, i*2) = flipud (interp1(filtered, (-2+table(fastPixels+1:fastPixels*2)), 'linear'));
            case 2
                out(j, :, i*2-1) =hermnu(filtered, table(1:fastPixels));
                out(j, :, i*2) = flipud (hermnu(filtered, (table(fastPixels+1:fastPixels*2))));
        end


    end
end
toc



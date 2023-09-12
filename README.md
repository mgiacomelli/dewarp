# dewarp
AVX accelerated resonant scanner dewarping e.g. for high speed two photon or confocal imaging.  Compared to other methods used, the approach here is optimal in terms of maximizing shot-noise limited SNR.  Higher order interpolation is also used to preserve higher spatial frequencies.  This code is forked off of our microscope software, and is designed to be run real-time with minimum processing latency (few tens of milliseconds).  Code is self-contained and should be straightforward to interface into other microscope software.  Although designed around an Alazar digitizer, we have also used it on our custom USB digitizers as well and it should work on any 4 or 16 channel digitizer returning 10, 12, 14 or 16 bit sample data.

# algorithm
Publication forthcoming. 

# versions
The algorithm is implemented in matlab and c using Intel's AVX intrinsics.  See:  https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html 

The c implementation has two versions, one optimized for 4 channel processing and another optimized for 16 channel processing.  Both require the AVX512-FP16 extensions, which are currently officially supported on only the Intel Sapphire Rapids processors, although they can be unofficially enabled on the 12th gen Alder Lake.  They will likely also be supported on the upcoming AMD Zen 5 processors and then most new CPUs thereafter.  I have mostly working versions that work on older processors (although at much lower performance), but have not yet had time to clean them up.

# compiling
VS2022 project files are provided.  Initially the code was targeted at compilation using clang (either version included in VS or the free Intel fork) due to bugs in the VS compiler, but the VS team responded quickly and as of August 2023, they were able to get the VS compiler working (although it is not as fast as clang).  

# performance
With full FIR filtering enabled to maximize sensitivity, an Intel 12700k w/ DDR4-3000 gives 1400 MP/s in clang and 1250 MP/s in MSOC for 16 channel dewarping.  Performance scales fairly well out to about 3 threads before memory bandwidth limitations become a factor, so roughly 4000 MP/s multithreaded.  With a 12 KHz resonant scanner, you're limited to about 50 MP/s per channel, or about 800 MP/s across 16 channels, so this is about 5x real-time.  Disabling FIR filtering will increase performance roughly 50%.  

# viewing the dewarped data
The library is intended to plugin into OpenGL or similar graphics libraries, so it generates one buffer (texture) per spectral channel which is intended to be loaded into a GPU for display.  For testing the library saves these to disk.  An example matlab script (viewData.m) is provided for viewing them. 

# forward/backward scan alignment
A comnmon problem with bidirectional scanning is the alignment between the forward and backward scans, which are ideally at zero delay relative to one another (meaning eactly half of each buffer is forward scan, half backwards scan).  In reality since each sample is only a few nanoseconds long there will be a few samples of misalignment due to the trigger signal not being perfectly aligned with the galvo turn around, scanner frequency drift, length of cables in the system, etc.  The "shift" variable in the script is a floating point value that represents how many samples of delay should be added.  As a floating point value, fraction sample shifts are allowed, and in general this number will not be an integer.  It can be calculated exactly by cross-correlation of the forward and backward scans.  

# example data
Data files from two different microscopes are provided in the releases section.  There is also a script providied for turning 16 channel data into 4 channel data (converter16to4chan.m).  


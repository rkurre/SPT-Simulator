# SPT-Simulator

Matlab software for simulation of two-dimensional single-color single-molecule tracking as well as 
sptPALM and PAINT eperiments. Single emitter signals are randomly distributed 
and can diffuse over time based on 2D random walks. Simulation accounts for using 
an EMCCD camera (e.g. Andor iXon Ultra) for single-molecule imaging which 
typically introduces excess noise at high EMGAIN (in our experiments: EMGAIN=300 
for iXon 897 Ultra). Excess noise is simulated by an additional noise factor of sqrt(2). 
Single-emitter simulation is based on Poisson statistics for excitation and pixelated
detection of point-spread functions. The program accounts for additional noise contributions
by photonelectron conversion, background fluorescence noise as well as experimentally determined camera noise using dark frames. 

This program generates a tif-stack of simulated data and a ground-truth table of each single emitter as csv-file. A Matlab m-file containing all variables of the simulation is also exported.

# System Requirements (tested environment)
- Operating System: Microsoft Windows 11 Home Version 23H2 (Build 22631.3527)
- MATLAB Version: 9.9.0.1467703 (R2020b) 
- Statistics and Machine Learning Toolbox, Version 12.0 (R2020b)

# Containing files:
- SPT_Simulator.fig: graphical user interface.
- SPT_Simulator.m: graphical user interface (Matlab source code).
- simFrames.m: source code for simulating single-molecule imaging stacks.
- example data output

# Installation guide
- Install Matlab and required toolboxes (installation time < 1 h)
- extract ZIP-file 'SPT_Simulator.zip'

# Instructions to run:
1. Start Matlab 2020b
2. In Matlab browser, navigate to software directory (\...\SPT_Simulator\)
3. run 'SPT_Simulator.m' in Matlab
4. Set parameters (details below) and push "Simulate" to generate tif-stacks and ground-truth data.

Expected run time with default parameters less than 1 minute, using a computer with AMD Ryzen 7 4700U or similar

# Input parameters:
## Imaging parameters:
- Exposure time: time interval of image acquisition and simulation (typical 5-50 ms).
- Imaging width: width of camera field of view in pixel (typical 64-512).
- Imaging Heigth: height of camera field of view in pixel (typical 64-512).
- Pixel size: pixel size of final simulated image in units of nanometer (typical 100 nm).
- Activation cycles: Cycles of photoactivation (sptPALM), for PAINT equals total number of frames; for SPT equals 1.
- Frames per cycle: Number of frames per photoactivation cycle. For PAINT equals 1 and for SPT equals total number of frames.
- Detection wavelength: wavelength of emitter fluorescence in units of micrometer.
- Numerical aperture: Numerical aperture of objective lens.
- Quantum efficiency: quantum efficiency of detector.
- Electron conversion: Photoelectron conversion factor of detector in photoelectrons per count (here: 0.015 electron per count).
- Camera Offset: camera offset of detector in units of digital counts. Can be measured by mean of dark frames.
- Camera noise: camera offset of detector in units of digital counts. Can be measured by standard deviation of dark frames.
- Illumination profile: Gaussian illumination profile with defined sigma or flat-top illumination profile.
- Axial penetration depth of evanescent field.

## Emitter parameters:
- Emitter intensity: Integrated signal of single emitter in units of photons (typical: 50-500 photons)
- Background intensity: Mean pixel background intensity (typical: 0-5 photons under TIRF conditions)
- Particle density: Global particle density of simulation. (typical: 0.1-0.5 for single-molecule localization)
- Lifetime: Lifetime of single-emitter fluorescence in units of frames based on single exponential decay (typical: 200-2000 frames).
- Oligomerization state: defines oligomer state, 1: monomer, 2: dimer, 3:trimer 4, tetramer, etc..
- Diffusion monomer: Diffusion constant of monomers in units of µm^2/s (typical 0.1-0.2 for transmembrane proteins)
- Diffusion oligomer:Diffusion constant of oligomer in units of µm^2/s (typical 0.05-0.1 for transmembrane proteins)
- Mean axial z-position of simulated particles above coverslip (e.g. position of plasma membrane).
- Standard deviation of z position variations along trajectories.
- Number of generated independent data sets based on same parameters.
- FileName (preffix): Name of the tif-stack.

# Output:
1. Single Tif-stack of simulated data
2. Ground-truth csv table: header = "frame","x [nm]","y [nm]","z [nm]","sigma [nm]","intensity [photon]","background [photon]';
3. Matlab M-file containing all parameters and main variables.
    
# Author
Rainer Kurre, PhD 

Email: rainer.kurre@uos.de

Center for Cellular Nanoanalytics and Division of Biophysics

Osnabrueck University, Barabastr. 11, D-49076 Osnabrück

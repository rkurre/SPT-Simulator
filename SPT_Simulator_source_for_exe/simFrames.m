function simFrames(dt, wAreaPxl, hAreaPxl, pxlSzIm, nCycles, nCycleFrames,...
    lambda, NA, QE, ph2cnts, camOffset, camNoise, useGauss, sigmaGauss,...
    pd, sigE, sigBG, pDen, ltDye, oligState, frOlig, Dm, Do, zAvg, zStd,...
    nSets, fileName, outputPath, handles)

    %Simulation of sptPALM/PAINT/SPT TIRFM eperiments. Single emitter 
    %signals are randomly distributed and can diffuse over time based on 
    %random walks. Axial position for each trajectory can vary over time 
    %representing e.g. plasma membrane ripples and is modelled by a 
    %normally-distributed random number generator. Please, note, that we 
    %don't modell the membrane, but we only introduce an additional
    %variation on single trajectory level to account for signal intensity
    %variations due to the expontential decay of the evanescent wave. 
    %To obtain smooth trajectories, we define equidistant nodes for random 
    %axial positions and interpolate remaining z positions by a spline 
    %function. Simulation is based on using an EMCCD camera (e.g. Andor 
    %iXon Ultra) for single-molecule imaging which typically introduces 
    %excess noise at high EMGAIN (in our experiments: EMGAIN=300 for iXon 
    %897 Ultra). Excess noise is simulated by an additional noise factor 
    %of 1.4. Illumination profile can be chosen as a predefined Gaussian
    %illumination (classical TIRFM) or a flat-top profile. The function 
    %generates a tif-stack and a ground-truth table as Matlab m-file as 
    %well as a csv-file.
    %
    %Input parameters
    %   dt: time interval of image acquisition and simulation.
    %   wAreaPxl: width of imaging area in pixel.
    %   hAreaPxl: height of imaging area in pixel.
    %   pxlSzIm: pixel size of final simulated image in units of micrometer.
    %   nCycles: Cycles of photoactivation (sptPALM), for PAINT equals total
    %            number of frames.
    %   nCycleFrame: Number of frames by photoactivation cycle. For PAINT
    %                nCycleFrame = 1 and for SPT equals total number of frames.
    %   lambda: wavelength of emitter in units of micrometer.
    %   NA: Numerical aperture of objective.
    %   QE: quantum efficiency of detector.
    %   ph2cnts: photonelectron conversion factor in photonelectrons per count.
    %   camOffset: camera offset in units of digital counts.
    %   camNoise: camera noise in units of digital counts.
    %   useGauss: True/(1): Gaussian illumination profile with defined sigma
    %             (sigmaGauss) or false/(0): Flattop illumination profile.
    %   sigmaGauss: Standard deviation sigma of Gaussian intensity profile.
    %   pd: Penetration depth of evanescent wave.
    %   sigE: Number of photons per emitter.
    %   sigBG: Number of photons generating continuous background.
    %   pDen: particle density in units of particles/µm^2.
    %   ltDye: characteristic lifetime of dye (assumption single exponential
    %          decay).
    %   oligState: Defines oligomer state, 1: monomer, 2: dimer, 3:trimer 4,
    %              tetramer, etc..
    %   frOlig: Fraction of particles being an oligomer. Set this
    %           parameter to 0 to simulate monomers only.
    %   Dm: %diffusion constant of monomers in units of µm^2/s.
    %   Do: %diffusion constant of oligomer in units of µm^2/s.
    %   zAvg: Mean axial z-position of simulated particles above coverslip
    %         (e.g. position of plasma membrane).
    %   zStd: Standard deviation of z position variations along
    %         trajectories.
    %   nSets: Number of generated independent datasets based on same 
    %          parameters
    %   fileName: Name of the tif file.
    %   outputPath: user defined output filepath
    %   handles: Handles of GUI for updating status and image in GUI.

    %Date: 05/01/2024
    %written by Rainer Kurre (rainer.kurre@uos.de)
    %Center for Cellular Nanoanalytics and division of biophysics
    %Osnabrueck University, Barabastr. 11, D-49076 Osnabrück
    
    
    rng('shuffle')
    
    exNoise = sqrt(2); %excess noise of EMCCD multiplication register
    %1D mean displacement of a free diffusing monomeric particle in µm
    mdm = sqrt(2*Dm*dt); 
    %1D mean displacement of a free diffusing dimeric particle in µm
    mdo = sqrt(2*Do*dt);
    %super resolution pixel size in units of micrometer
    pxlSzSR = 0.005; 
    %width and heigth of image in microns
    wArea = wAreaPxl*pxlSzIm;
    hArea = hAreaPxl*pxlSzIm;
    %mean number of particles per frame
    nPar = round(pDen*wArea*hArea);
    %Definition of PSF of single emitters
    psfFWHM = lambda/(2*NA);
    psfWin = 4*psfFWHM;
    psfWinSz = round(psfWin/pxlSzSR); %size of PSF in units of SR pixels
    if mod(psfWinSz,2) ~= 0
        psfWinSz = psfWinSz +1;
    end
    %calc standard deviation for PSF model
    sigx = 1.3*psfFWHM/2.355;
    sigy = sigx;
    %function of PSF modelling
    gaussPsf = @(x,y) exp(-x.^2./(2*sigx^2)-y.^2./(2*sigy^2));
    
    mImSR = ceil(hArea/pxlSzSR);
    nImSR = ceil(wArea/pxlSzSR);

    if mod(mImSR,2)~=0
        mImSR = mImSR + 1;
    end
    if mod(nImSR,2)~=0
        nImSR = nImSR + 1;
    end

    bin = round(pxlSzIm/pxlSzSR);
    mImSR = mImSR - mod(mImSR, bin);
    nImSR = nImSR - mod(nImSR, bin);
    
    %Definition of PSF grid
    psfWinHalfSz = round(psfWinSz/2);
    [psfx, psfy] = meshgrid(-psfWinHalfSz:psfWinHalfSz,...
        -psfWinHalfSz:psfWinHalfSz);
    psfx = pxlSzSR*psfx;
    psfy = pxlSzSR*psfy;
    
    %Definition of Gaussian illumination profile
    %size of illumination in units of SR pixels
    [illGaussX, illGaussY] = meshgrid(-mImSR/2:mImSR/2, -nImSR/2:nImSR/2);
    sigmaG = sigmaGauss/pxlSzIm*bin;
    illGauss = @(x,y) exp(-x.^2./(2*sigmaG^2)-y.^2./(2*sigmaG^2));
    illProfileGauss = illGauss(illGaussX, illGaussY);
    illProfileGaussAvg = mean(illProfileGauss(:));
    illCenterX = nImSR/2;
    illCenterY = mImSR/2;
    
    %Correction of emitter intensity by axial position
    sigEcorr =  sigE/exp(-zAvg/pd);
    
    %Definition of distance between nodes for z-position spline function
    dNodesM = round(2/mdm);
    dNodesO = round(2/mdo);

    
    for iSet=1:nSets
        %generate struct array containing tracks
        tracks = struct('xTr', {}, 'yTr', {}, 'zTr', {}, 'frame', {});
        nFluoroTotal = 0;
        maxTau = 0;
        nTr = 0;
        set(handles.status, 'String', ['Data set ',int2str(iSet),...
            ': Generation of tracks.']);
        drawnow;
        for n=1:nCycles
            %oligomerization
            nOlig = round(frOlig*nPar); %fraction of particles being an oligomer
            nMono = nPar-nOlig;
            nFluoro = nMono+oligState*nOlig;
            nFluoroTotal = nFluoroTotal + nFluoro;
            x0 = repmat(rand(nOlig,1)*wArea,oligState,1);
            y0 = repmat(rand(nOlig,1)*hArea,oligState,1);
            x0 = vertcat(x0, rand(nMono,1)*wArea);
            y0 = vertcat(y0, rand(nMono,1)*hArea);
            for i=1:nOlig
                nTr = nTr + 1;
                tau = sort(round(random('exp', ltDye,oligState,1)), 'descend');
                if (tau(1) > maxTau)
                    maxTau = tau(1);
                end
                dxLast = randn(tau(1),1)*mdo;
                dyLast = randn(tau(1),1)*mdo;
                xTr = x0(i)+cumsum(dxLast);
                yTr = y0(i)+cumsum(dyLast);
                nPoints = length(xTr);
                if nPoints<=dNodesO
                    zTr(1:nPoints)= normrnd(zAvg, zStd);
                else
                    nNodes = floor(nPoints/dNodesO)+1;
                    if mod(nPoints, dNodesO) == 0
                        nodes = zeros(1, nNodes);
                        nodes(1:end-1) = 1:dNodesO:nPoints;
                        nodes(end) = nPoints;
                    elseif mod(nPoints, dNodesO) == 1
                        nodes = 1:dNodesO:nPoints;
                    else
                        nNodes = nNodes+1;
                        nodes = zeros(1, nNodes);
                        nodes(1:end-1) = 1:dNodesO:nPoints;
                        nodes(end) = nPoints;
                    end
                    zNodes = normrnd(zAvg, zStd, 1, nNodes);
                    zTr = spline(nodes, zNodes, 1:nPoints);
                end
                zTr(zTr<0) = 0;
                tracks(nTr).xTr = xTr;
                tracks(nTr).yTr = yTr;
                tracks(nTr).zTr = zTr;
                tracks(nTr).frame = [(n-1)*nCycleFrames+1, (n-1)*nCycleFrames + ...
                    numel(xTr)];
                for j = 2:oligState
                    nTr = nTr + 1;
                    dx = dxLast(1:tau(j));
                    dy = dyLast(1:tau(j));
                    xTr = x0(i+(j-1)*nOlig)+cumsum(dx);
                    yTr = y0(i+(j-1)*nOlig)+cumsum(dy);
                    tracks(nTr).xTr = xTr;
                    tracks(nTr).yTr = yTr;
                    tracks(nTr).zTr = zTr(1:length(xTr));
                    tracks(nTr).frame = [(n-1)*nCycleFrames+1, ...
                        (n-1)*nCycleFrames + numel(xTr)];
                end
            end
            
            for i = 1:nMono
                nTr = nTr + 1;
                tau = round(random('exp', ltDye));
                if (tau > maxTau)
                    maxTau = tau(1);
                end
                dx = randn(tau,1)*mdm;
                dy = randn(tau,1)*mdm;
                xTr = x0(nOlig*oligState+i)+cumsum(dx);
                yTr = y0(nOlig*oligState+i)+cumsum(dy);
                nPoints = length(xTr);
                if nPoints<=dNodesM
                    zTr(1:nPoints)= normrnd(zAvg, zStd);
                else
                    nNodes = floor(nPoints/dNodesM)+1;
                    if mod(nPoints, dNodesM) == 0
                        nodes = zeros(1, nNodes);
                        nodes(1:end-1) = 1:dNodesM:nPoints;
                        nodes(end) = nPoints;
                    elseif mod(nPoints, dNodesM) == 1
                        nodes = 1:dNodesM:nPoints;
                    else
                        nNodes = nNodes+1;
                        nodes = zeros(1, nNodes);
                        nodes(1:end-1) = 1:dNodesM:nPoints;
                        nodes(end) = nPoints;
                    end
                    zNodes = normrnd(zAvg, zStd, 1, nNodes);
                    zTr = spline(nodes, zNodes, 1:nPoints);
                end
                zTr(zTr<0) = 0;
                tracks(nTr).xTr = xTr;
                tracks(nTr).yTr = yTr;
                tracks(nTr).zTr = zTr;
                tracks(nTr).frame = [(n-1)*nCycleFrames+1, (n-1)*nCycleFrames + ...
                    numel(xTr)];
            end
        end
            
        nFrames = nCycles*nCycleFrames;
        locs = struct('xCoord', [], 'yCoord', [], 'zCoord', [], 'amp', []);
        
        for i=1:nFrames
            set(handles.status, 'String', ['Data set ',int2str(iSet),...
            ': Generation of Locs. ', 'Frame ', int2str(i), ...
            ' of ', int2str(nFrames),' done.']);
            drawnow;
            locs(i).xCoord = [];
            for j=1:nFluoroTotal
                minTrFrame = tracks(j).frame(1);
                maxTrFrame = tracks(j).frame(2);
                if i >= minTrFrame && i <= maxTrFrame
                    locs(i).xCoord = vertcat(locs(i).xCoord, tracks(j).xTr(...
                        i-minTrFrame+1));
                    locs(i).yCoord = vertcat(locs(i).yCoord, tracks(j).yTr(...
                        i-minTrFrame+1));
                    locs(i).zCoord = vertcat(locs(i).zCoord, tracks(j).zTr(...
                        i-minTrFrame+1));
                    locs(i).amp = vertcat(locs(i).amp, sigEcorr);
                end
            end
        end
        
        %cell array for exporting ground truth data
        header = {'"frame"','"x [nm]"','"y [nm]"','"z [nm]"','"sigma [nm]"',...
            '"intensity [photon]"','"background [photon]"'};
        GT = cell(size(locs,2), 6);
        id = 1;
        for i=1:nFrames
            im = zeros(mImSR, nImSR);
            xLocs = locs(i).xCoord;
            if any(xLocs)
                nLocs = numel(xLocs);
                for j=1:nLocs
                    %get xy coordinates in units of SR pixels
                    xLoc = round(locs(i).xCoord(j)/pxlSzSR);
                    yLoc = round(locs(i).yCoord(j)/pxlSzSR);
                    %get z coordinate in units of nanometers
                    zLoc = round(locs(i).zCoord(j));
                    %calc amplitude for PSF model, do illumination correction
                    %if requested
                    if useGauss
                        illCorrF = illGauss(yLoc-illCenterY, xLoc-illCenterX)...
                            /illProfileGaussAvg;
                    else
                        illCorrF = 1;
                    end
                    amp = illCorrF*locs(i).amp(j).*exp(-zLoc/pd);
                    psfIm = gaussPsf(psfx,psfy);
                    psfIm = random('poisson',amp)*(psfIm./(sum(psfIm(:))));
                    xInd = xLoc - psfWinHalfSz:xLoc + psfWinHalfSz;
                    yInd = yLoc - psfWinHalfSz:yLoc + psfWinHalfSz;
                    xIndPSF=xInd>0 & xInd<= nImSR;
                    yIndPSF=yInd>0 & yInd<= mImSR;
                    xInd(xInd<=0 | xInd>nImSR)=[];
                    yInd(yInd<=0 | yInd>mImSR)=[];
                    im(yInd, xInd) = im(yInd, xInd) + psfIm(yIndPSF, xIndPSF);
                    %write ground truth data
                    GT{id, 1}=i;
                    GT{id, 2}=(locs(i).xCoord(j))*1000;
                    GT{id, 3}=(locs(i).yCoord(j))*1000;
                    GT{id, 4}=(locs(i).zCoord(j));
                    GT{id, 5}=sigx*1000;
                    GT{id, 6}=amp;
                    GT{id, 7}=sigBG;
                    id = id + 1;
                end
            end
            %digital binning of emitter signals
            im = sum(reshape(im,bin,[]) ,1);
            im = reshape(im,mImSR/bin,[]).';
            im = sum(reshape(im,bin,[]),1);
            im = reshape(im,nImSR/bin,[]).';
            %adding mean fluorescent background signal
            im = im + sigBG;
            %convert to photoelectrons and add detection noise
            im = im*QE+im*QE*exNoise^2-random('poisson',im*QE*exNoise^2);
            %calc image in units of digital counts
            im = im.*ph2cnts;
            %add camera offset and noise
            pause(0.005);
            im = im + normrnd(camOffset, camNoise, size(im,1), size(im,2));
            imwrite(uint16(im), [outputPath,'\', fileName, '_', num2str(iSet, ...
                '%0.2d'),'.tif'], 'WriteMode', 'Append');
            set(handles.status, 'String', ['Data set ',int2str(iSet),...
            ': Generation of image frames. Frame ', ...
                int2str(i), ' of ', int2str(nFrames),' done.']);
            drawnow;
            if mod(i-1,10)==0
                imagesc(handles.im, im);
                axis(handles.im, 'image');
                drawnow;
            end
        end
        set(handles.status, 'String', ['Data set ',int2str(iSet),...
            ': Save ground-truth data...']);
        drawnow;
        GT_table = [header; GT];
        writetable(cell2table(GT_table), [outputPath, '\', fileName, '_GT_', ...
            num2str(iSet, '%0.2d'), '.csv'], 'WriteVariableNames', false);
        save([outputPath,'\', fileName, '_', num2str(iSet, '%0.2d'),'.mat'], ...
            'dt', 'wAreaPxl', 'hAreaPxl', 'pxlSzIm','nCycles', ...
            'nCycleFrames', 'lambda', 'NA', 'QE', 'ph2cnts', ...
            'camOffset', 'camNoise', 'useGauss', 'sigmaGauss', 'pd', ...
            'sigE', 'sigBG', 'pDen', 'ltDye', 'oligState', 'frOlig', ...
            'Dm', 'Do', 'zAvg', 'zStd', 'fileName','tracks', 'locs', ...
            'GT_table', 'exNoise', 'mdm', 'mdo', 'pxlSzSR', 'nPar', ...
            'psfFWHM', 'sigx', 'sigy', 'mImSR', 'nImSR', 'sigEcorr', ...
            'dNodesM', 'dNodesO', 'nFrames', '-mat');
        set(handles.status, 'String', ['Data set ',int2str(iSet),...
            ': Save ground-truth data... Done.']);
        drawnow;
    end
    clear all;
end


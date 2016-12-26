%
% make_roi_image_mask.m
%
% Use this program to create an ROI in the SPM-environment. The program
%   has a functionality that is similar to the froi-function make_roi of
%   N. Knouf (Kanwisher lab, MIT)
%
% Usage:
%
% Use SPM to call up xSPM structure for a specific contrast
%       (e.g., objects-scrambled)
% After that you can type make_roi to show these voxels on coronals
% You have to change anatImage for your subject (always wmr*.img)
%    and maybe the rangeOfSlices too, depending on what you want to see
%       (now it only shows the posterior part of the brain)
%
% All significant voxels are selected (shown in red) by default
% This is what you can do:
% - press 'c' on the keyboard: de-select all voxels
%       (de-selected significant voxels are shown in yellow)
% - click the left mouse button to select individual voxels/regions
% - press 'b' on the keyboard:
%       change the size of the region selected with a left mouse click
% - press 'i' on the keyboard:
%       compute the intersection between the selected region and
%       the original set of significant voxels
% - press 's' on the keyboard: save matlab-file (do NOT type extension)
% - press 'e' on the keyboard to end the program
%
% The saved structure makeROI contains all the selected voxels in the field
%   makeROI.selected, with the anatomical coordinates in makeROI.selected.anatXYZ.
%   You can also find all originally significant voxels in other fields.
%
% Written by Hans Op de Beeck
% All questions and bug reports should be addressed to Hans Op de Beeck
%       (hans.opdebeeck@psy.kuleuven.be)
%
% Edited by Jessica Bulthé
% This script now automatically
%

clear anatIme
clear makeROI

% spm fmri;

subjectIDs = [{'C1' 'C2' 'C3' 'C4' 'C5' 'C6' 'C7' 'C8' 'C9' 'C10' 'C11' 'C12' 'C13' 'C14' 'C15' 'C16' 'C18' 'C19' 'C20' 'C21' 'C22' 'C25' 'C26' 'C27' ...
    'D1' 'D2' 'D3' 'D4' 'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14' 'D15' 'D16' 'D18' 'D19' 'D20' 'D21' 'D22' 'D25' 'D26' 'D27'}];     

treshold = 0.1;

for s = 1:size(subjectIDs,2)
    subject = subjectIDs(s); 

    HomeDir = 'E:\Research\Dyscalculie Studie\fMRI\';
    AnatomieDir = [HomeDir 'Anatomical Scans\'];
    SPMDir = [HomeDir 'Statistiek\' char(subject) '\Loc\'];
    ROIDir = [HomeDir 'ROIs\Anatomisch ' num2str(treshold) '\' char(subject) '\'];
    MaskDir = [HomeDir 'ROIs\Masks\'];
    ScriptDir = [HomeDir 'scripts\make anatomical rois\'];
    
    mkdir(ROIDir);

    % Which masks are present
    Masks = dir(MaskDir);
    Masks(1:2) = [];
    num_Masks = size(Masks, 1);

    % Load in SPM the SPM.mat file with the masks
%     for m = 1:num_Masks
    for m = 18:18

        current_mask = Masks(m).name;

        % put the spm.mat file masked with the current_mask in graphics window
        % of SPM
        run([ScriptDir 'load_spm_graphics.m']);

        % Load in anatomy file
        anatomie_file = dir([AnatomieDir 'wr' char(subject) '_ANATOMIE_*.nii']);
        anatImage = spm_vol([AnatomieDir anatomie_file.name]);

        % The slices
        rangeOfSlices = [5 95];

        % Leave the rest to peace
        volumeSize= [size(anatImage.private.dat,1) size(anatImage.private.dat,2) size(anatImage.private.dat,3)];
        telActive = size(xSPM.XYZ,2);

        % process slice numbers and positions
        numSlices = rangeOfSlices(2) - rangeOfSlices(1) + 1;
        perDim = sqrt(numSlices);
        if round(perDim) == perDim
            XimNum = round(perDim);
            YimNum = round(perDim);
        else
            XimNum = round(perDim);
            YimNum = round(perDim) + 1;
        end;

        telSlice = 0;
        for slY=1:YimNum,
            for slX=1:XimNum,
                telSlice = telSlice + 1;
                if telSlice <= numSlices
                    sliceXYimPos(telSlice,1:2) = [slX slY];
                end;
            end;
        end;

        for sl=1:numSlices,
            startXpos = (sliceXYimPos(sl,1)-1)*volumeSize(1);
            startYpos = (sliceXYimPos(sl,2)-1)*volumeSize(3);
            hulp = flipdim(flipdim(squeeze(anatImage.private.dat(:,rangeOfSlices(1)+sl-1,:))',1),2);
            anatIm(startYpos+1:startYpos+volumeSize(3),startXpos+1:startXpos+volumeSize(1)) = hulp;
        end;
        for y=1:size(anatIm,1),
            for x=1:size(anatIm,2),
                if isnan(anatIm(y,x))
                    anatIm(y,x) = 0;
                end;
            end;
        end;
        colorRange = [0 max(max(anatIm))]
        figure(6)
        %A=imshow(anatIm,[0 1700]);
        A=imshow(anatIm,colorRange);
        hold on;

        makeROI.allPoints = xSPM.XYZ;
        makeROI.nAllPoints = size(makeROI.allPoints,2);
        telVox = 0;
        for p=1:makeROI.nAllPoints,
            if xSPM.XYZ(2,p)>=rangeOfSlices(1) &  xSPM.XYZ(2,p)<=rangeOfSlices(2);
                telVox = telVox + 1;
                makeROI.allInFOV(1:3,telVox) = xSPM.XYZ(:,p);
            end;
        end;
        makeROI.nInFov = telVox

        % if there are less then 20 voxels, continue to next mask
        if makeROI.nInFov < 20
            continue;
        end

        for p=1:makeROI.nInFov,
            sl = makeROI.allInFOV(2,p)-rangeOfSlices(1)+1;
            startXpos = (sliceXYimPos(sl,1)-1)*volumeSize(1);
            startYpos = (sliceXYimPos(sl,2)-1)*volumeSize(3);
            makeROI.figCoord(:,p) = [startYpos+(volumeSize(3)-makeROI.allInFOV(3,p)+1)  startXpos+(volumeSize(1)-makeROI.allInFOV(1,p)+1)];
        end;

        P = plot(makeROI.figCoord(2,:),makeROI.figCoord(1,:),'r.');
        set(P,'MarkerSize',1);
        makeROI.selected.figCoord = makeROI.figCoord;
        makeROI.selected.number = size(makeROI.selected.figCoord,2);

        boxSize=3;

        fprintf('Computing intersection\n');
        hulpSelected = makeROI.selected;
        makeROI.selected.figCoord = [];
        makeROI.selected.number = 0;
        for sel=1:hulpSelected.number,
            for v=1:makeROI.nInFov,
                if hulpSelected.figCoord(1,sel)==makeROI.figCoord(1,v) &  hulpSelected.figCoord(2,sel)==makeROI.figCoord(2,v)
                    makeROI.selected.number = makeROI.selected.number + 1;
                    makeROI.selected.figCoord(1:2,makeROI.selected.number) = hulpSelected.figCoord(:,sel);
                end;
            end;
        end;
        A=imshow(anatIm,colorRange);
        hold on;
        P = plot(makeROI.figCoord(2,:),makeROI.figCoord(1,:),'y.');
        set(P,'MarkerSize',1);
        P = plot(makeROI.selected.figCoord(2,:),makeROI.selected.figCoord(1,:),'r.');
        set(P,'MarkerSize',1);

        % save
        % first clean up doubles
        double(1:makeROI.selected.number)=0;
        for i=1:makeROI.selected.number-1,
            for j=i+1:makeROI.selected.number,
                if makeROI.selected.figCoord(1,i)==makeROI.selected.figCoord(1,j) & makeROI.selected.figCoord(2,i)==makeROI.selected.figCoord(2,j)
                    double(j) = 1;
                end;
            end;
        end;

        hulpSelected = makeROI.selected;
        makeROI.selected.figCoord = [];
        makeROI.selected.number = 0;
        for sel=1:hulpSelected.number,
            if double(sel) == 0
                makeROI.selected.number = makeROI.selected.number + 1;
                makeROI.selected.figCoord(1:2,makeROI.selected.number) = hulpSelected.figCoord(:,sel);
            end;
        end;

        % than add anatomical coordinates to the figCoord
        makeROI.selected.anatXYZ = [];
        for sel=1:makeROI.selected.number,
            for v=1:makeROI.nInFov,
                if makeROI.selected.figCoord(1,sel)==makeROI.figCoord(1,v) &  makeROI.selected.figCoord(2,sel)==makeROI.figCoord(2,v)
                    makeROI.selected.anatXYZ(:,sel) = makeROI.allInFOV(1:3,v);
                end;
            end;
        end;

        fileName = current_mask(1:end-4);
        output=[ROIDir fileName '.mat'];
        save(output, 'makeROI');
        fprintf('ROI saved\n');
    end
end
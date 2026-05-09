%Enviroment setting:
javaaddpath 'C:\Program Files\MATLAB\R2021a\java\mij.jar'. %Path of mij.jar
javaaddpath 'C:\Users\zhang\Downloads\fiji-win64\Fiji.app\jars\ij-1.53f.jar'. %Path of Fiji

MIJ.start('C:\Users\zhang\Downloads\fiji-win64\Fiji.app');
Datafolder = 'D:\260319_possp2026-03-19_16-10-28\S\'; %Path of the raw data folder
addpath(genpath('./AIM'))
N = dir(fullfile(Datafolder,'*.tif'));

% Change the parameters in the following lines:
% Line 36: Whether to remove certain frames (e.g. to avoid saturation) or not
% Line 37: Camera set-up
% Line 38: Gaussian fitting setting in ThunderSTORM
% Line 39: Post-fitting analysis in ThunderSTORM (Merge localisations to avoid replicates)
% Line 40: Post-fitting analysis in ThunderSTORM (Remove uncertainty-derived replicates)
% Line 41: Post-fitting analysis in ThunderSTORM (Remove noise by density filter)
% Line 63: Drift correction: AIM parameter
% Line 78: Final visualisation setting
% Line 98: Radius for localisation clustering

for i = 1:length(N)
    fname= N(i).name;
    candlefile=[Datafolder,fname,'P/mor.csv'];
    if isfile(candlefile)
        continue;
    end
    mkdir(Datafolder,[fname,'P']);
    savefolder=[Datafolder,fname,'P/'];
    fprintf(1, 'Now reading %s\n', fname);
    filepath=[Datafolder,fname];
    loclistpath=[savefolder,'loclist.csv'];
    pathcommand=['path=',filepath];
    loclistpathadd=['filepath=',loclistpath];
    loclistpathcommand=[loclistpathadd,' fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=false offset=false saveprotocol=true x=true y=true bkgstd=false id=true uncertainty=true frame=true detections=true")'];
    MIJ.run('Open...', pathcommand);
    MIJ.run("Slice Remover", "first=1 last=50 increment=1");
    MIJ.run("Camera setup", "offset=250 isemgain=false photons2adu=8 pixelsize=103.5");
    MIJ.run("Run analysis","filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Local maximum] connectivity=8-neighbourhood threshold=1.2*std(Wave.F1)estimator=[PSF: Integrated Gaussian] sigma=1.6 fitradius=2 method=[Weighted Least squares] full_image_fitting=false mfaenabled=false renderer=[Averaged shifted histograms] dxforce=false magnification=10.0 dx=5.0 colorizez=false threed=false dzforce=false repaint=50");
    MIJ.run("Show results table", "action=merge zcoordweight=0.1 offframes=1 dist=33.0 framespermolecule=0");
    MIJ.run("Show results table", "action=duplicates distformula=uncertainty");
    MIJ.run("Show results table", "action=density neighbors=5 radius=50.0 dimensions=2D");
    MIJ.run("Export results", loclistpathcommand);
    MIJ.close("*");
    MIJ.closeAllWindows;
    modcsvfile=[savefolder,'loclist.csv'];
    M=readtable(modcsvfile);
    Mx = M.x_nm_;
    My = M.y_nm_;
    Mf = M.frame;
    Mnew=[Mx,My,Mf];
    csvwrite(modcsvfile,Mnew);
    modcsvfile=[savefolder,'loclist.csv'];
    data = readtable(modcsvfile);
    data1=table2array(data);
    Localizations=horzcat(data1(:,3),data1(:,1),data1(:,2))

%% Asign the experimental frame number and localization coordinates (x and y) to localizations variable
    F = Localizations(:,1) ; %frame_id
    X = Localizations(:,2) ; % x position of the localization coordinates
    Y = Localizations(:,3) ; % y position of the localization coordinates

%% AIM drift correction and parameter settings
    trackInterval = 80; % time interval for drift tracking, Unit: frames 
    t_start = tic;
    [LocAIM, AIM_Drift] = AIM(Localizations, trackInterval);
    AIM_time = toc(t_start)
    
    header={'frame','x [nm]', 'y [nm]'};
    combined_data = [header; num2cell(LocAIM)];
    save([savefolder, fname(1:end-4) '_AIM.mat'],'LocAIM', 'AIM_Drift');
    LocAIMoutput=[savefolder,'LocAIMoutput.csv'];
    writetable(array2table(LocAIM, 'VariableNames', header), LocAIMoutput);

%% Save all data
    corcsvfile=[savefolder,'LocAIMoutput.csv'];
    importcommand=['filepath=',corcsvfile,' fileformat=[CSV (comma separated)] livepreview=false rawimagestack= startingframe=1 append=false'];
    MIJ.run("Import results",importcommand);
    MIJ.run("Visualization", "imleft=3.0 imtop=2.4 imwidth=524.2 imheight=524.6 renderer=[Averaged shifted histograms] magnification=10.0 colorizez=false threed=false shifts=2");
    spsavecommand=['Tiff..., path=',savefolder,'superres.tif'];
    MIJ.run('Save', spsavecommand);
    MIJ.run("Label Boundaries");
    MIJ.run("Connected Components Labeling", "connectivity=4 type=[16 bits]");
    MIJ.run("Analyze Regions", "area perimeter circularity centroid centroid equivalent_ellipse ellipse_elong. max._feret_diameter");
    spsavecommand=['path=',savefolder,'mor.csv'];
    MIJ.run('Results...', spsavecommand);
    MIJ.close("*");
    MIJ.closeAllWindows;
    mor = readtable([savefolder,'mor.csv']);
    locs = readmatrix([savefolder,'LocAIMoutput.csv']);
    mor_pos  = [mor.Centroid_X, mor.Centroid_Y] * 1000;
    loc_pos = locs(:, 2:3);
    obj_func = @(s) mean_knn_dist(mor_pos + s, loc_pos);
    initial_guess = [0, 0];
    options = optimset('Display', 'final', 'TolX', 1e-2);
    best_shift = fminsearch(obj_func, initial_guess, options);
    shifted_mor = mor_pos + best_shift;
    [idx, dists] = knnsearch(shifted_mor, loc_pos);
    max_radius = 350; %cluster setting 
    valid_mask = dists < max_radius;
    counts = histcounts(idx(valid_mask), 1:height(mor)+1);
    mor.LocalizationCount = counts';
    mor_with_counts=[savefolder,'mor_with_counts.csv'];
    writetable(mor, mor_with_counts);

   
end

function d_mean = mean_knn_dist(centroids, points)
    % Finds distance to the single nearest neighbor for every point
    [~, d] = knnsearch(centroids, points);
    d_mean = mean(d);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This code constitutes an important tool for obtaining the results in
%    the manuscript entitled  "Phase analysis simulating the Takeda method 
%    to obtain a 3D profile of SARS-CoV-2 cells", submitted to Transactions
%    on Pattern Analysis and Applications-Springer.

%    We recall that the ROIs are built from 3D tomographic data in the 
%    research work "Live imaging of SARS-CoV-2 infection in mice reveals
%    that neutralizing antibodies require Fc function for optimal efficacy"
%    with DOI: 10.1016/j.immuni.2021.08.015, where we build several frames 
%    (video image captures of the tomographic images)and the phase was
%    obtained for each frame simulating the Takeda method. By integrating 
%    the phase of different frames of the SARS-CoV-2 cell we obtain the 3D
%    profile.
%
%    Correspondings Authors:
%    jesus.arriagahdz@correo.buap.mx                  b.cuevas@irya.unam.mx
%
%    This algorithm is a simple routine innovating the main files 
%    required to obtain the results along with the application of our main 
%    function "Results.m".


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Fringes pattern reading (we include the matrix containing the used 
%    fringes pattern data (data.m). The "load()" instruction works according
%    the the file location in each computer.
load('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Programas\franjas_150.mat')

%    Variable "Z" call in Results(.,.,z)
z = data;

%    Essential parameters construction according to a first ROI (regarding
%    the size of the medical image) containing the remaining ROIs
image = imread('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-1 Frames\3D_Cell_1_11.png');
imagen = rgb2gray(image);
imagen = double(imagen);
CC = 512;
[a,b] = size(imagen);
C = 5;
A = a*C;
B = b*C;

%    All the ROI frames are previously stored, segmenting the SARS-CoV-2
%    cells obtained from the tomographic image video.  Subsequently, we
%    open the ROI image files using the following routine (opening all the 
%    ordered images with a single instruction in the following loop).

%    Medical image , Cell 1 Frames 
lee_archivos = dir('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-1 Frames\*.png'); 

%    Medical image, Cell 2 Frames
%lee_archivos = dir('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-2 Frames\*.png'); 

J = length(lee_archivos);

for k = 1:J
    %    Instruction for the medical image or cell 1 frames
    jpgFilename = strcat('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-1 Frames\3D_Cell_1_', num2str(k), '.png');
    
    %    Instruction for the medical image or cell 2 frames
    %jpgFilename = strcat('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-2 Frames\3D_Cell_2_', num2str(k), '.png');
    
    imageData = imread(jpgFilename);
    Frames(:,:,:,k) = imageData;
end

%    All the masks are previously stored, generated for each ROI frame,
%    segmenting the cell. Subsequently, we use the following routine to open
%    the files of each mask synchronized with each image frame.

%    Binary mask for each cell 1 Frame 
lee_mascaras = dir('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-1 Mask Frames\*.mat'); 

%    Binary mask for each cell 2 Frame
%lee_mascaras = dir('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-2 Mask Frames\*.mat'); 

K = length(lee_archivos);
for m = 1:K
    %    Instruction for the Binary mask corresponding to cell 1
    mascara = strcat('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-1 Mask Frames\Mask_3D_Cell_1_', num2str(m), '.mat');

    %    Instruction for the Binary mask corresponding to cell 2
    %mascara = strcat('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\Cel-2 Mask Frames\Mask_3D_Cell_2_', num2str(m), '.mat');
        
    mascaras = load(mascara);
    Masks(:,:,m) = cell2mat(struct2cell(mascaras));
end

%    In the following routines previous parameters are generalized, which 
%    are regarded as essential ones. Direct segmentation elements are
%    applied over each capture performed in the video of the tomographic
%    image. We segment two times to delimit the ROI of the SARS-CoV-2 cell.
%    All the ROI frames were previously stored, corresponding to the ROI
%    segmenting the cell.

for k = 1:J
    imag = rgb2gray(Frames(:,:,:,k));
    Imag = double(imag);
    Image = (imresize(Imag, [A B]));

    %   ROI segmentation containing the Cell 1 information only
    cropp = imcrop(Image, [3067.5 2978.5 772 803]);
    cropp2 = imcrop(cropp, [46.5 59.5 705 705]);
    Crop = (imresize(cropp2, [CC CC]));

    %   ROI segmentation containing the Cell 2 information only
    %cropp = imcrop(Image, [214.5 1889.5 1132 1139]);
    %cropp2 = imcrop(cropp, [40.5 50 988 988]);
    %Crop = (imresize(cropp2, [CC CC]));
    
    %    Binary masks synchronization for each frame.  The latter
    %    delimitates specifically the information of the cell in each frame.
    %    This constitutes our optical segmentation.
    Mask_Crop = Crop.*Masks(:,:,k);

    %    Results from the phase unwrapped from the information of the
    %    SARS-CoV-2 cell for each frame. We refer to several instructions 
    %    since the function admits as an input variable in X the images 
    %    with masks in double format and GT, whereas in variable Y could be
    %    constant. However, variable Z should be always the proposed pattern.
    %Result(:,:,k) = Results(Crop,Masks(:,:,k),z);
    Result(:,:,k) = Results(Mask_Crop,1,z);
    ResultB(:,:,k) = Result(:,:,k).*Masks(:,:,k);


    %    In the following routine the edges with value equal to zero are 
    %    removed from the binary mask segmenting precisely the cell 
    %    information. The latter constitutes the first version of the 
    %    3D model from accurate cell data.
   
    Aux = ResultB(:,:,k);
    for i=1:CC
        for j=1:CC
            if Aux(i,j) == 0
                Result_0(i,j,k) = NaN;
            else
                Result_0(i,j,k) = Aux(i,j);
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Routine for generating a preliminary version of a matrix with size
%	mxmxz with the frames obtained from Video. We build the elements 
%	from the corresponding level curves or isophotes, integrating the
%	data in a 3D matrix.
%	3D matrix 512x512x21, for cell 1 and 512x512x22 for cell 2

%	Rule fo three for the scale, 59nm for cell 1 and 83nm for cell 2
a = 512;

%    Cell 1
XX = (a/465)*59;

%    Cell 2
%XX = (a/440)*83;

[X, Y]=meshgrid(0:(XX/a):XX);
X=imresize(X,[a a]);    Y=imresize(Y,[a a]);

figure;
hold on
 for ii = 1 : J
     zz = (3.5*J) - (3.5*ii);        %  separation between several z planes
     [~,h] = contour(X,Y,Result_0(:,:,ii),30); % plot contour at the bottom
     h.ContourZLevel = zz;
     
 end
hold off
view(3);
xlabel('X [nm]'); ylabel('Y [nm]'); zlabel('Z [nm]', 'Rotation',0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Routine for reading the interpolation files, performed by our algorithm
%    called "2DFrame_Fit" developed in Python, described in a first version 
%    in our previous work "Two-dimensional Legendre polynomials as a basis 
%    for interpolation of data to optimize the solution of the irradiance 
%    transport equation analyzed as a boundary problem on surfaces testing"
%    with DOI: 10.1364/AO.58.005057.  3000 files are read, generated by the
%    interpolator. However, for the sake of functionality and optimal PC
%    operation, we consider 500 interpolation images files, generating first
%    50 files from the 21 frames, subsequently 100, then 170, reaching 250
%    and finally 500 interpolation images with the least percentage error.
%    These generate a matrix with size mxmxm (cubic matrix) to visualiza a
%    3D model for the SARS-CoV-2 cell, the calculadte virion by our proposal.

%    Images for Cell 1 fit
lee_FitImages = dir('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\LegendrePol_3D images\Cell_1\*.mat'); 

%    Images for Cell 2 fit
%lee_FitImages = dir('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\LegendrePol_3D images\Cell_2\*.mat'); 
M = length(lee_FitImages);


for m = 1:M
    
    %    Instruction for opening the fit files of cell 1
    file = strcat('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\LegendrePol_3D images\Cell_1\MaskTemp', num2str(m), '.mat');
    
    %    Instruction for opening the fit files of cell 2
    %file = strcat('C:\Users\jesus\Documents\Gerchberg-Saxton COVID-19 Micro\Takeda simulado en covid\Fig 3D analy\LegendrePol_3D images\Cell_2\MaskTemp2_', num2str(m), '.mat');        
    
    EstrAux = load(file);
    fitimage = EstrAux.mat;
    %    Identification of the images generated by LegendrePol_2D_Fit
    t=m;
    aux_0 = fitimage;
    aux_1 = imresize(aux_0, [M M]);
    FitImage(:,:,t) = aux_1;
    %    Images with edges and zeros in NaN variables
    
    Aux = FitImage(:,:,t);
    for i=1:M
        for j=1:M
            if Aux(i,j) == 0
                NanImages(i,j,t) = NaN;
            else
                NanImages(i,j,t) = Aux(i,j);
            end
        end
    end

end

%   Rescaling the measurements, 59nm for cell 1 and 83nm for cell 2, 
%   cubic matrix
%   Cell 1
XX = (M/410)*59;

%    Cell 2
%XX = (M/425)*83;

[X, Y]=meshgrid(0:(XX/a):XX);
X=imresize(X,[M M]);    Y=imresize(Y,[M M]);

figure;
hold on
 for ii = 1 : M
     zz = (J) - (ii);
     [~,h] = contour(X,Y,NanImages(:,:,ii),30);% plot contour at the bottom
     h.ContourZLevel = zz;
     
 end
hold off
view(3);
xlabel('X [nm]'); ylabel('Y [nm]'); zlabel('Z [nm]', 'Rotation',0)

%    Instruction to see the 3D model of the SARS-CoV-2 cell
volumeViewer(NanImages);
 




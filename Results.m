function [Res] = Results(x,y,z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This code constitutes an important tool for obtaining the results in
%    the manuscript entitled  "Phase analysis simulating the Takeda method
%    to obtain a 3D profile of SARS-CoV-2 cells", submitted to Transactions
%    on Image Processing-IEEE.
%
%    Correspondings Authors:
%    Dr. Jesus Alonso Arriaga Hernandez
%    jesus.arriagahdz@correo.buap.mx;    dr.j.a.arriaga.hernandez@gmail.com
%    Dra. Bolivia Teresa Cuevas Otahola
%    b.cuevas@irya.unam.mx;                      b.cuevas.otahola@gmail.com
%
%    This function contains three input data of square matrices of the form 
%    MxM (in our results we studied the M=512 case), which should be images
%    transformed into Gray Tones (GT) of double type. In variable "X" we 
%    have the ROI (Region of Interest) in the medical image, in "Y" we have
%    a binary mask delimiting the ROI cell (SARS-CoV-2 cell). In variable 
%    "Z" we have the periodic pattern inducing periodicity, which are as 
%    important as the variable X, since otherwise it is not possible to 
%    simulate the Takeda's Method. The output from the function (Res) is
%    the unwrapped phase of the cell in the ROI of the medica image. We 
%    bear in mind that we aim to obtain the phase in all frames from the 
%    obtained phase of the information from a single frame in a 3D 
%    tomographic image, to finally obtain a 3D model.
%
%    In this work, we analyze the majority of objects digitization or 3D
%    reconstruction techniques, considering from these the Takeda method 
%    (Fringes Projection or Fourier Transform Method) by its simplicity. We
%    induce a periodicity in the images to subsequently simulate the Takeda
%    Method, filtering the periodic noise and finally unwrapping the phase.
%    This work was carried out in collaboration with Dra. Lilia Cedillo and 
%    Dr. Ygnacio Martínez, who greatly supported the interpretation of the 
%    results. The videos of the tomographic image were obtained in the
%    research entitled "Live imaging of SARS-CoV-2 infection in mice
%    reveals that neutralizing antibodies require Fc function for optimal
%    efficacy";  with DOI: 10.1016/j.immuni.2021.08.015


%%
%    ROI variables introduction (x) and binary mask of SARS-CoV-2 cell in 
%    the ROI (y)
    
Mask = y;
image = x;

CC = 512;
a = CC;                      b= CC;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imag_test = Mask.*image;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2D Fourier Tranform ordered towards the center
fft_TestIma = fftshift(fft2(Imag_test));

[X, Y]=meshgrid(-((a)/2):1:(a)/2);
d1 = 13;
X=imresize(X,[a a]);    Y=imresize(Y,[a a]);

%    2D simple filter to remove periodic noise in the image ROI (X)
%    applying the 2D Fourier transform.
circ1=((X).^2+(Y).^2<=(d1).^2);
Filtro = double(circ1);
FT_fil = (fft_TestIma).*Filtro;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY_ifftimag = (ifft2(FT_fil));
XY_IfftImag = abs(XY_ifftimag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Removal and identification of the periodic noise
A = max(max(XY_IfftImag));
A2 = max(max(Imag_test));
if A < A2
    ruido = XY_IfftImag - double(Imag_test);
else
    ruido = double(Imag_test) - XY_IfftImag;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                      The Takeda Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Periodicity induction with 85 fringes pattern, Variable Z

z2 = z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b] = size(z2);
imagtest = (imresize(XY_IfftImag, [a b]));
z = z2.*imagtest;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2D Fourier transform ordered towards the center of the Takeda method,
%    for the image ROI with induced periodicity by the variable Z and only 
%    the periodic image Z
FT_period_image = fftshift(fft2(z));
FT_patron = fftshift(fft2(z2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y]=meshgrid(-((a)/2):1:(a)/2);
d2 = d1;

%    Image filtering by means of the Takeda method, for the ROI image
%    with periodicity and only periodic image Z
X=imresize(X,[a a]);    Y=imresize(Y,[a a]);
circ1=((X - (342-(a/2))).^2+(Y).^2<=(d2).^2);
Filtro = double(circ1);
FT_fil_image = (FT_period_image).*Filtro;
FT_fil_patron = (FT_patron).*Filtro;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY_ift_Filimag = (ifft2(FT_fil_image));
XY_ift_patron = z2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Takeda Shiffted
for i = 1 : a
    for j = 1 : b
        if z2(i,j) == 1;
            Data_shiff(i,j) = 0;
        else
            Data_shiff(i,j) = 1;
        end
    end
end

FT_shif_patron = fftshift(fft2(Data_shiff));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y]=meshgrid(-((a)/2):1:(a)/2);
d3 = d2;
X=imresize(X,[a a]);    Y=imresize(Y,[a a]);
circ1=((X - (342-(a/2))).^2+(Y).^2<=(d3).^2);
Filtro = double(circ1);
FT_fil_image = (FT_period_image).*Filtro;
FT_fil_shif_patron = (FT_shif_patron).*Filtro;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY_ift_shif_patron = Data_shiff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Phase acquisition by means of the Taked Method, application of Eq. 5 
%    of the manuscript

S_A = ((XY_ift_shif_patron).*(XY_ift_Filimag));
S_B = ((XY_ift_patron).*(XY_ift_Filimag));
Super_product = (S_A)-(S_B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Removal with edges with values equal to zero of the built masks.
%    Accurate segmentation procedure of the phase of the SARS-CoV-2 cell
%    delimited by NAN variables
result_0 = Mask.*abs(Super_product);
for i=1:CC
    for j=1:CC
        if result_0(i,j) == 0
            Result_0(i,j) = NaN;
        else
            Result_0(i,j) = result_0(i,j);
        end
    end
end

%    The following line is optional, since they could be removed, given 
%    that they only show the resulting image.  Moreover, the graphic type 
%    can be changed by "mesh();", instead of the "imagesc();"
figure; imagesc(Result_0);
title('First final result with NAN');

%    Final result and output variable.
Res = Result_0;

end
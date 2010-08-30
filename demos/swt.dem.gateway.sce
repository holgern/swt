
demopath = get_absolute_file_path("swt.dem.gateway.sce");

subdemolist=['one dimensional wavelet transform', 'swt1d.sce';..
 '1-D Wavelet decomposition using wavedec()','wavedec.dem.sce';..
 '1-D de-noising of test signal using wden()','wden.dem.sce';..
 '1-D scale 2 freq demo using scal2frq()','scale2freq.dem.sce';..
 '1-D dwtmode demo', 'dwtmode.dem.sce';..
'two dimensional wavelet transform', 'swt2d.sce';..
'three dimensional wavelet transform', 'swt3d.sce';..
'one dimensional dualtree complex wavelet transform','cowt1d.sce';..
'two dimensional dualtree complex wavelet transform','cowt2d.sce';..
'continous wavelet transform', 'cwt.sce';..
//'stationary wavelet transform', 'swtswt.sce';..
'image processing', 'image.sce';];


subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================
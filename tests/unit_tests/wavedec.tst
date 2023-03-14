// -------------------------------------------------------------------------
// SWT - Scilab wavelet toolbox
// Copyright (C) 2010-2014  Holger Nahrstaedt
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//-------------------------------------------------------------------------
//
//  <-- NO CHECK ERROR OUTPUT -->


// dwt1d  Test

loadmatfile("-mat",fullfile(get_swt_path(),"tests/unit_tests/Data.mat"));

// wavedec
// haar
[cA1,cD1]=dwt(s1,'haar');
[cA2,cD2]=dwt(cA1,'haar');
[cA3,cD3]=dwt(cA2,'haar');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,'haar');
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );



// db1
for N=1:36
wname="db"+sprintf("%d",N);
[cA1,cD1]=dwt(s1,wname);
[cA2,cD2]=dwt(cA1,wname);
[cA3,cD3]=dwt(cA2,wname);
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,wname);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );
end;



// coif1
for N=1:17
wname="coif"+sprintf("%d",N);
[cA1,cD1]=dwt(s1,wname);
[cA2,cD2]=dwt(cA1,wname);
[cA3,cD3]=dwt(cA2,wname);
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,wname);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );
end;
for N=2:20
wname="sym"+sprintf("%d",N);

// sym4
[cA1,cD1]=dwt(s1,wname);
[cA2,cD2]=dwt(cA1,wname);
[cA3,cD3]=dwt(cA2,wname);
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,wname);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );
end;



// bior1.1
bior_fam=["bior1.1","bior1.3", "bior1.5", "bior2.2", "bior2.4", "bior2.6",...
"bior2.8", "bior3.1", "bior3.3", "bior3.5", "bior3.7",...
"bior3.9", "bior4.4", "bior5.5", "bior6.8"];
// bior1.1
for i=1:max(size(bior_fam))
[cA1,cD1]=dwt(s1,bior_fam(i));
[cA2,cD2]=dwt(cA1,bior_fam(i));
[cA3,cD3]=dwt(cA2,bior_fam(i));
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,bior_fam(i));
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );
end;

for i=1:max(size(bior_fam))
// wavedec type 2 iput
Lo_D=rand(1,20,'normal');
Hi_D=rand(1,20,'normal');
[cA1,cD1]=dwt(s1,Lo_D,Hi_D);
[cA2,cD2]=dwt(cA1,Lo_D,Hi_D);
[cA3,cD3]=dwt(cA2,Lo_D,Hi_D);
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );
[cA1,cD1]=dwt(s1,bior_fam(i));
[cA2,cD2]=dwt(cA1,bior_fam(i));
[cA3,cD3]=dwt(cA2,bior_fam(i));
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );
end;

// type 3
// asymh
for i=1:max(size(bior_fam))
dwtmode('asymh','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','asymh');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','asymh');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','asymh');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );

// symw
dwtmode('symw','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','symw');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','symw');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','symw');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );


// asymw
dwtmode('asymw','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','asymw');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','asymw');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','asymw');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );


// zpd
dwtmode('zpd','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','zpd');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','zpd');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','zpd');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );


// sp0
dwtmode('sp0','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','sp0');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','sp0');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','sp0');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );


// sp1
dwtmode('sp1','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','sp1');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','sp1');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','sp1');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );


// ppd
dwtmode('ppd','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','ppd');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','ppd');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','ppd');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );

// per
dwtmode('per','nodisp');
[cA1,cD1]=dwt(s1,bior_fam(i),'mode','per');
[cA2,cD2]=dwt(cA1,bior_fam(i),'mode','per');
[cA3,cD3]=dwt(cA2,bior_fam(i),'mode','per');
[Lo_D,Hi_D]=wfilters(bior_fam(i),'d');
c0=[cA3 cD3 cD2 cD1];
l0=[length(cA3) length(cD3) length(cD2) length(cD1) length(s1)];
[c,l]=wavedec(s1,3,Lo_D,Hi_D);
assert_checkalmostequal ( c , c0 , %eps );
assert_checkalmostequal ( l , l0 , %eps );

dwtmode("symh",'nodisp');
end;

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

// idwt
// haar
[cA,cD]=dwt(x1,'haar');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('haar');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
r=idwt(cA,cD,'haar');
assert_checkalmostequal ( r , x0 , %eps, %eps*10);
[cA,cD]=dwt(x2,'haar');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('haar');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
r=idwt(cA,cD,'haar');
assert_checkalmostequal ( r , x0 , %eps, %eps*10);
[cA,cD]=dwt(s1,'haar');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('haar');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
r=idwt(cA,cD,'haar');
assert_checkalmostequal ( r , x0 , %eps, %eps*10);


// db1 - db10
for N=1:10
  wname="db"+sprintf("%d",N);
  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10);
  [cA,cD]=dwt(x2,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10);
  [cA,cD]=dwt(s1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10);
end;






// coif1 - coif5
for N=1:5
  wname="coif"+sprintf("%d",N);
  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10  );
  [cA,cD]=dwt(x2,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10  );
  [cA,cD]=dwt(s1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10  );
end





// sym4 - sym10
for N=4:10
  wname="sym"+sprintf("%d",N);
  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10 );
  [cA,cD]=dwt(x2,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10 );
  [cA,cD]=dwt(s1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps , %eps*10);
end;



// bior1.1
wnames=['bior1.1','bior1.3','bior1.5','bior2.2','bior2.4','bior2.6','bior2.8','bior3.1','bior3.3','bior3.5','bior3.7','bior3.9'];
for N=1:12
  wname=wnames(N);
  [cA,cD]=dwt(x1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10);
  [cA,cD]=dwt(x2,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10);
  [cA,cD]=dwt(s1,wname);
  [Lo_D,Hi_D,Lo_R,Hi_R]=wfilters(wname);
  a0=conv(dyadup(cA),Lo_R);
  d0=conv(dyadup(cD),Hi_R);
  x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
  r=idwt(cA,cD,wname);
  assert_checkalmostequal ( r , x0 , %eps, %eps*10);
end;



// // type 2
[cA,cD]=dwt(x1,'bior3.9');
Lo_R=rand(1,20,'normal');
Hi_R=rand(1,20,'normal');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
r=idwt(cA,cD,Lo_R,Hi_R);
assert_checkalmostequal ( r , x0 , %eps, %eps*10);
[cA,cD]=dwt(x2,'bior3.9');
Lo_R=rand(1,20,'normal');
Hi_R=rand(1,20,'normal');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
r=idwt(cA,cD,Lo_R,Hi_R);
assert_checkalmostequal ( r , x0 , %eps, %eps*10);
[cA,cD]=dwt(s1,'bior3.9');
Lo_R=rand(1,20,'normal');
Hi_R=rand(1,20,'normal');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,2*length(cA)-length(Lo_R)+2);
r=idwt(cA,cD,Lo_R,Hi_R);
assert_checkalmostequal ( r , x0 , %eps, %eps*10);


// type 3
[cA,cD]=dwt(x1,'sym8');
r=idwt(cA,cD,'sym8',50);
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym8');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,50);
assert_checkalmostequal ( r , x0 , %eps, %eps*10);
[cA,cD]=dwt(x2,'sym8');
r=idwt(cA,cD,'sym8',50);
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym8');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,50);
assert_checkalmostequal ( r , x0 , %eps, %eps*10);
[cA,cD]=dwt(s1,'sym8');
r=idwt(cA,cD,'sym8',50);
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym8');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,50);
assert_checkalmostequal ( r , x0 , %eps, %eps*10);
Lo_R=rand(1,50,'normal');
Hi_R=rand(1,50,'normal');
a0=conv(dyadup(cA),Lo_R);
d0=conv(dyadup(cD),Hi_R);
x0=wkeep(a0+d0,50);
r=idwt(cA,cD,Lo_R,Hi_R,50);
assert_checkalmostequal ( r , x0 , %eps, %eps*10);



// type 4
[cA,cD]=dwt(x1,'db7');
a0=idwt(cA,cD,'db7','mode','symh');
a1=idwt(cA,cD,'db7','mode','symw');
a2=idwt(cA,cD,'db7','mode','asymh');
a3=idwt(cA,cD,'db7','mode','asymw');
a4=idwt(cA,cD,'db7','mode','sp0');
a5=idwt(cA,cD,'db7','mode','sp1');
a6=idwt(cA,cD,'db7','mode','zpd');
a7=idwt(cA,cD,'db7','mode','ppd');
a8=idwt(cA,cD,'db7','mode','per');
[Lo_R,Hi_R]=wfilters('db7','r');
aa0=conv(dyadup(cA),Lo_R);
dd0=conv(dyadup(cD),Hi_R);
x0=wkeep(aa0+dd0,2*length(cA)-length(Lo_R)+2);
xx1=wkeep(aa0+dd0,2*length(cA));
assert_checkalmostequal ( a0 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a1 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a3 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a4 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a5 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a6 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a7 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a8 , xx1 , %eps, %eps*10);



// type 5
[cA,cD]=dwt(x1,'db7');
a0=idwt(cA,cD,'db7',50,'mode','symh');
a1=idwt(cA,cD,'db7',50,'mode','symw');
a2=idwt(cA,cD,'db7',50,'mode','asymh');
a3=idwt(cA,cD,'db7',50,'mode','asymw');
a4=idwt(cA,cD,'db7',50,'mode','sp0');
a5=idwt(cA,cD,'db7',50,'mode','sp1');
a6=idwt(cA,cD,'db7',50,'mode','zpd');
a7=idwt(cA,cD,'db7',50,'mode','ppd');
a8=idwt(cA,cD,'db7',50,'mode','per');
[Lo_R,Hi_R]=wfilters('db7','r');
aa0=conv(dyadup(cA),Lo_R);
dd0=conv(dyadup(cD),Hi_R);
x0=wkeep(aa0+dd0,50);
assert_checkalmostequal ( a0 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a1 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a3 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a4 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a5 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a6 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a7 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a8 , x0 , %eps, %eps*10);


// type 6
[cA,cD]=dwt(x1,'db7');
Lo_R=rand(1,14,'normal');
Hi_R=rand(1,14,'normal');
a0=idwt(cA,cD,Lo_R,Hi_R,50,'mode','symh');
a1=idwt(cA,cD,Lo_R,Hi_R,50,'mode','symw');
a2=idwt(cA,cD,Lo_R,Hi_R,50,'mode','asymh');
a3=idwt(cA,cD,Lo_R,Hi_R,50,'mode','asymw');
a4=idwt(cA,cD,Lo_R,Hi_R,50,'mode','sp0');
a5=idwt(cA,cD,Lo_R,Hi_R,50,'mode','sp1');
a6=idwt(cA,cD,Lo_R,Hi_R,50,'mode','zpd');
a7=idwt(cA,cD,Lo_R,Hi_R,50,'mode','ppd');
a8=idwt(cA,cD,Lo_R,Hi_R,50,'mode','per');
aa0=conv(dyadup(cA),Lo_R);
dd0=conv(dyadup(cD),Hi_R);
x0=wkeep(aa0+dd0,50);
assert_checkalmostequal ( a0 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a1 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a3 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a4 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a5 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a6 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a7 , x0 , %eps, %eps*10);
assert_checkalmostequal ( a8 , x0 , %eps, %eps*10);


// column vector
[cA,cD]=dwt(x1','sym9');
[cA1,cD1]=dwt(x1,'sym9');
a0=idwt(cA,cD,'sym9');
a1=idwt(cA',cD,'sym9');
a2=idwt(cA,cD','sym9');
a3=idwt(cA',cD','sym9');

assert_checkalmostequal ( cA , cA1 , %eps, %eps*10);
assert_checkalmostequal ( cD , cD1 , %eps, %eps*10);
assert_checkalmostequal ( a1 , a0 , %eps, %eps*10);
assert_checkalmostequal ( a2 , a0 , %eps, %eps*10);
assert_checkalmostequal ( a3 , a0 , %eps, %eps*10);

// void coef
[cA,cD]=dwt(x1,'sym9');
a0=idwt(cA,[],'sym9');
d0=idwt([],cD,'sym9');
[Lo_R,Hi_R]=wfilters('sym9','r');
aa0=wkeep(conv(dyadup(cA),Lo_R),2*length(cA)-length(Lo_R)+2);
dd0=wkeep(conv(dyadup(cD),Hi_R),2*length(cA)-length(Lo_R)+2);

assert_checkalmostequal ( a0 , aa0 , %eps, %eps*10);
assert_checkalmostequal ( d0 , dd0 , %eps, %eps*10);

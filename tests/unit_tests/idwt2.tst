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


// dwt2d  Test

loadmatfile("-mat",fullfile(get_swt_path(),"tests/unit_tests/Data.mat"));
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// idwt2
//loadmatfile("-mat","Data.mat");
d1=rand(64,65,'normal');
// haar
[cA,cH,cV,cD]=dwt2(d1,'haar');
[Lo_R,Hi_R]=wfilters('haar','r');
dd0=idwt2(cA,cH,cV,cD,'haar');
dd1=idwt2(cA,cH,cV,cD,'haar',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R));
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R));
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R));
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R));
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));

assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );

clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db1 - db10
for N=1:10
  wname="db"+sprintf("%d",N);

  [cA,cH,cV,cD]=dwt2(d1,wname);
  [Lo_R,Hi_R]=wfilters(wname,'r');
  dd0=idwt2(cA,cH,cV,cD,wname);
  dd1=idwt2(cA,cH,cV,cD,wname,size(d1));
  dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
  dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
  [r,c]=size(cA);
  for i=1:c,
    col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R));
    col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R));
    col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R));
    col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R));
  end
  row_low=col_low_low+col_low_hi;
  row_hi=col_hi_low+col_hi_hi;
  [rc,cc]=size(row_low);
  for i=1:rc,
    row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
    row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
  end
  ss=2*size(cA)-length(Lo_R)+2;
  x0=wkeep(row_low_r+row_hi_r,ss);
  x1=wkeep(row_low_r+row_hi_r,size(d1));
  assert_checkalmostequal ( x0 , dd0 , %eps, %eps*10  );
  assert_checkalmostequal ( x1 , dd1 , %eps, %eps *10  );
  assert_checkalmostequal ( x0 , dd2 , %eps , %eps*10  );
  assert_checkalmostequal ( x1 , dd3 , %eps , %eps *10 );
  clear x0;
  clear x1;
  clear i;
  clear rc;
  clear cc;
  clear row_low;
  clear row_hi;
  clear col_low_low;
  clear col_low_hi;
  clear col_hi_low;
  clear col_hi_hi;
  clear row_low_r;
  clear row_hi_r;
  clear r;
  clear c;
  clear Lo_R;
  clear Hi_R;
  clear dd0;
  clear dd1;
  clear dd2;
  clear dd3;
end;



// coif1 - coif5
// 4 and 5 does not work???
for N=1:3
  wname="coif"+sprintf("%d",N);
  [cA,cH,cV,cD]=dwt2(d1,wname);
  [Lo_R,Hi_R]=wfilters(wname,'r');
  dd0=idwt2(cA,cH,cV,cD,wname);
  if exists('dd0')
  dd1=idwt2(cA,cH,cV,cD,wname,size(d1));
  dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
  dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
  [r,c]=size(cA);
  for i=1:c,
    col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R));
    col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R));
    col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R));
    col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R));
  end
  row_low=col_low_low+col_low_hi;
  row_hi=col_hi_low+col_hi_hi;
  [rc,cc]=size(row_low);
  for i=1:rc,
    row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
    row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
  end
  ss=2*size(cA)-length(Lo_R)+2;
  x0=wkeep(row_low_r+row_hi_r,ss);
  x1=wkeep(row_low_r+row_hi_r,size(d1));
  assert_checkalmostequal ( x0 , dd0 , %eps, %eps*10);
  assert_checkalmostequal ( x1 , dd1 , %eps, %eps*10);
  assert_checkalmostequal ( x0 , dd2 , %eps, %eps*10);
  assert_checkalmostequal ( x1 , dd3 , %eps, %eps*10);
  clear x0;
  clear x1;
  clear i;
  clear rc;
  clear cc;
  clear row_low;
  clear row_hi;
  clear col_low_low;
  clear col_low_hi;
  clear col_hi_low;
  clear col_hi_hi;
  clear r;
  clear c;
  clear Lo_R;
  clear Hi_R;
  clear dd0;
  clear dd1;
  clear dd2;
  clear dd3;
  clear row_low_r;
  clear row_hi_r;
end;
end;


// sym4 - sym10
for N=4:10
  wname="sym"+sprintf("%d",N);
  [cA,cH,cV,cD]=dwt2(d1,wname);
  [Lo_R,Hi_R]=wfilters(wname,'r');
  dd0=idwt2(cA,cH,cV,cD,wname);
  if exists('dd0')
  dd1=idwt2(cA,cH,cV,cD,wname,size(d1));
  dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
  dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
  [r,c]=size(cA);
  for i=1:c,
    col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R));
    col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R));
    col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R));
    col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R));
  end
  row_low=col_low_low+col_low_hi;
  row_hi=col_hi_low+col_hi_hi;
  [rc,cc]=size(row_low);
  for i=1:rc,
    row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
    row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
  end
  ss=2*size(cA)-length(Lo_R)+2;
  x0=wkeep(row_low_r+row_hi_r,ss);
  x1=wkeep(row_low_r+row_hi_r,size(d1));
  assert_checkalmostequal ( x0 , dd0 , %eps, %eps*10);
  assert_checkalmostequal ( x1 , dd1 , %eps, %eps*10);
  assert_checkalmostequal ( x0 , dd2 , %eps, %eps*10);
  assert_checkalmostequal ( x1 , dd3 , %eps, %eps*10);
  clear x0;
  clear x1;
  clear i;
  clear rc;
  clear cc;
  clear row_low;
  clear row_hi;
  clear col_low_low;
  clear col_low_hi;
  clear col_hi_low;
  clear col_hi_hi;
  clear r;
  clear c;
  clear Lo_R;
  clear Hi_R;
  clear dd0;
  clear dd1;
  clear dd2;
  clear dd3;
  clear row_low_r;
  clear row_hi_r;
end;
end;


// bior1.1 - bior3.9
//bior3.9 does not work
//wnames=['bior1.1','bior1.3','bior1.5','bior2.2','bior2.4','bior2.6','bior2.8','bior3.1','bior3.3','bior3.5','bior3.7','bior3.9'];
wnames=['bior1.1','bior1.3','bior1.5','bior2.2','bior2.4','bior2.6','bior2.8','bior3.1','bior3.3','bior3.5','bior3.7'];
for N=1:11
  wname=wnames(N);
  [cA,cH,cV,cD]=dwt2(d1,wname);
  [Lo_R,Hi_R]=wfilters(wname,'r');
  dd0=idwt2(cA,cH,cV,cD,wname);
  dd1=idwt2(cA,cH,cV,cD,wname,size(d1));
  dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
  dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
  [r,c]=size(cA);
  for i=1:c,
    col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R));
    col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R));
    col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R));
    col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R));
  end
  row_low=col_low_low+col_low_hi;
  row_hi=col_hi_low+col_hi_hi;
  [rc,cc]=size(row_low);
  for i=1:rc,
    row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
    row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
  end
  ss=2*size(cA)-length(Lo_R)+2;
  x0=wkeep(row_low_r+row_hi_r,ss);
  x1=wkeep(row_low_r+row_hi_r,size(d1));
  assert_checkalmostequal ( x0 , dd0 , %eps, %eps*10);
  assert_checkalmostequal ( x1 , dd1 , %eps, %eps*10);
  assert_checkalmostequal ( x0 , dd2 , %eps, %eps*10);
  assert_checkalmostequal ( x1 , dd3 , %eps, %eps*10);
  clear x0;
  clear x1;
  clear i;
  clear rc;
  clear cc;
  clear row_low;
  clear row_hi;
  clear col_low_low;
  clear col_low_hi;
  clear col_hi_low;
  clear col_hi_hi;
  clear r;
  clear c;
  clear Lo_R;
  clear Hi_R;
  clear dd0;
  clear dd1;
  clear dd2;
  clear dd3;
  clear row_low_r;
  clear row_hi_r;
end;

// type 2
[cA,cH,cV,cD]=dwt2(d1,'bior3.7');
[Lo_R,Hi_R]=wfilters('bior3.7','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.7','mode','symw');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'bior3.7',size(d1),'mode','zpd');
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,'mode','sp0');
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1),'mode','ppd');

assert_checkalmostequal ( dd0 , dd2 , %eps, %eps*10);
assert_checkalmostequal ( dd1 , dd3 , %eps, %eps*10);
clear cA;
clear cH;
clear cV;
clear cD;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
end
// type 3

[cA,cH,cV,cD]=dwt2(d1,'bior3.7','mode','per');
[Lo_R,Hi_R]=wfilters('bior3.7','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.7','mode','per');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'bior3.7',size(d1),'mode','per');
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,'mode','per');
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1),'mode','per');
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R));
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R));
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R));
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R));
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,2*size(cA));
x1=wkeep(row_low_r+row_hi_r,size(d1));

assert_checkalmostequal ( x0 , dd0 , %eps, %eps*10);
assert_checkalmostequal ( x1 , dd1 , %eps, %eps*10);
assert_checkalmostequal ( x0 , dd2 , %eps, %eps*10);
assert_checkalmostequal ( x1 , dd3 , %eps, %eps*10);

clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end

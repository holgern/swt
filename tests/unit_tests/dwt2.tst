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

// dwt 2d

loadmatfile("-mat",fullfile(get_swt_path(),"tests/unit_tests/Data.mat"));
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;

// haar
[cA,cH,cV,cD]=dwt2(d1,'haar');
[Lo_D,Hi_D]=wfilters('haar','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');

assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db1
[cA,cH,cV,cD]=dwt2(d1,'db1');
[Lo_D,Hi_D]=wfilters('db1','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;
clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db2
[cA,cH,cV,cD]=dwt2(d1,'db2');
[Lo_D,Hi_D]=wfilters('db2','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps*1000);
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db3
[cA,cH,cV,cD]=dwt2(d1,'db3');
[Lo_D,Hi_D]=wfilters('db3','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps*10000);
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db4
[cA,cH,cV,cD]=dwt2(d1,'db4');
[Lo_D,Hi_D]=wfilters('db4','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps*10000000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db5
[cA,cH,cV,cD]=dwt2(d1,'db5');
[Lo_D,Hi_D]=wfilters('db5','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db6
[cA,cH,cV,cD]=dwt2(d1,'db6');
[Lo_D,Hi_D]=wfilters('db6','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');

clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db7
[cA,cH,cV,cD]=dwt2(d1,'db7');
[Lo_D,Hi_D]=wfilters('db7','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db8
[cA,cH,cV,cD]=dwt2(d1,'db8');
[Lo_D,Hi_D]=wfilters('db8','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db9
[cA,cH,cV,cD]=dwt2(d1,'db9');
[Lo_D,Hi_D]=wfilters('db9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps *10000);
assert_checkalmostequal ( ccA , cca , %eps *10000);
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps *10000);
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// db10
[cA,cH,cV,cD]=dwt2(d1,'db10');
[Lo_D,Hi_D]=wfilters('db10','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps*1000);
assert_checkalmostequal ( cD , ccd , %eps*100000);
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps*100000);
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// coif1
[cA,cH,cV,cD]=dwt2(d1,'coif1');
[Lo_D,Hi_D]=wfilters('coif1','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;


// coif2
[cA,cH,cV,cD]=dwt2(d1,'coif2');
[Lo_D,Hi_D]=wfilters('coif2','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps*100000);
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps*100000);
assert_checkalmostequal ( cD , ccd , %eps*100000);
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps*100000);
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// coif3
[cA,cH,cV,cD]=dwt2(d1,'coif3');
[Lo_D,Hi_D]=wfilters('coif3','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps*100000);
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// coif4
[cA,cH,cV,cD]=dwt2(d1,'coif4');
[Lo_D,Hi_D]=wfilters('coif4','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// coif5
[cA,cH,cV,cD]=dwt2(d1,'coif5');
[Lo_D,Hi_D]=wfilters('coif5','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sym4
[cA,cH,cV,cD]=dwt2(d1,'sym4');
[Lo_D,Hi_D]=wfilters('sym4','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sym5
[cA,cH,cV,cD]=dwt2(d1,'sym5');
[Lo_D,Hi_D]=wfilters('sym5','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sym6
[cA,cH,cV,cD]=dwt2(d1,'sym6');
[Lo_D,Hi_D]=wfilters('sym6','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sym7
[cA,cH,cV,cD]=dwt2(d1,'sym7');
[Lo_D,Hi_D]=wfilters('sym7','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sym8
[cA,cH,cV,cD]=dwt2(d1,'sym8');
[Lo_D,Hi_D]=wfilters('sym8','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sym9
[cA,cH,cV,cD]=dwt2(d1,'sym9');
[Lo_D,Hi_D]=wfilters('sym9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sym10
[cA,cH,cV,cD]=dwt2(d1,'sym10');
[Lo_D,Hi_D]=wfilters('sym10','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior1.1
[cA,cH,cV,cD]=dwt2(d1,'bior1.1');
[Lo_D,Hi_D]=wfilters('bior1.1','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior1.3
[cA,cH,cV,cD]=dwt2(d1,'bior1.3');
[Lo_D,Hi_D]=wfilters('bior1.3','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior1.5
[cA,cH,cV,cD]=dwt2(d1,'bior1.5');
[Lo_D,Hi_D]=wfilters('bior1.5','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior2.2
[cA,cH,cV,cD]=dwt2(d1,'bior2.2');
[Lo_D,Hi_D]=wfilters('bior2.2','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior2.4
[cA,cH,cV,cD]=dwt2(d1,'bior2.4');
[Lo_D,Hi_D]=wfilters('bior2.4','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior2.6
[cA,cH,cV,cD]=dwt2(d1,'bior2.6');
[Lo_D,Hi_D]=wfilters('bior2.6','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior2.8
[cA,cH,cV,cD]=dwt2(d1,'bior2.8');
[Lo_D,Hi_D]=wfilters('bior2.8','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior3.1
[cA,cH,cV,cD]=dwt2(d1,'bior3.1');
[Lo_D,Hi_D]=wfilters('bior3.1','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior3.3
[cA,cH,cV,cD]=dwt2(d1,'bior3.3');
[Lo_D,Hi_D]=wfilters('bior3.3','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior3.5
[cA,cH,cV,cD]=dwt2(d1,'bior3.5');
[Lo_D,Hi_D]=wfilters('bior3.5','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior3.7
[cA,cH,cV,cD]=dwt2(d1,'bior3.7');
[Lo_D,Hi_D]=wfilters('bior3.7','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// bior3.9
[cA,cH,cV,cD]=dwt2(d1,'bior3.9');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D);
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// symh

[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','symh');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','symh');
m_ex=wextend(2,'symh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// symw
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','symw');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','symw');
m_ex=wextend(2,'symw',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// asymh
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','asymh');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','asymh');
m_ex=wextend(2,'asymh',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// asymw
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','asymw');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','asymw');
m_ex=wextend(2,'asymw',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// zpd
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','zpd');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','zpd');
m_ex=wextend(2,'zpd',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sp0
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','sp0');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','sp0');
m_ex=wextend(2,'sp0',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// sp1
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','sp1');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','sp1');
m_ex=wextend(2,'sp1',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// ppd
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','ppd');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','ppd');
m_ex=wextend(2,'ppd',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+length(Lo_D)-1 c+length(Lo_D)-1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

// per
[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','per');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(d1,Lo_D,Hi_D,'mode','per');
m_ex=wextend(2,'per',d1,length(Lo_D));
[r,c]=size(d1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r c]),'m');
cch=dyaddown(wkeep(col_low_hi,[r c]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r c]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r c]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;

dd1=rand(51,51,'normal');
[cA,cH,cV,cD]=dwt2(dd1,'bior3.9','mode','per');
[Lo_D,Hi_D]=wfilters('bior3.9','d');
[ccA,ccH,ccV,ccD]=dwt2(dd1,Lo_D,Hi_D,'mode','per');
m_ex=wextend(2,'per',dd1,length(Lo_D));
[r,c]=size(dd1);
[rex,cex]=size(m_ex);
for i=1:rex,
  row_low(i,:)=conv(Lo_D,m_ex(i,:));
  row_hi(i,:)=conv(Hi_D,m_ex(i,:));
end
[rrex,ccex]=size(row_low);
for i=1:ccex,
  col_low_low(:,i)=conv(Lo_D,row_low(:,i));
  col_low_hi(:,i)=(conv(Hi_D,row_low(:,i)));
  col_hi_low(:,i)=(conv(Lo_D,row_hi(:,i)));
  col_hi_hi(:,i)=(conv(Hi_D,row_hi(:,i)));
end
cca=dyaddown(wkeep(col_low_low,[r+1 c+1]),'m');
cch=dyaddown(wkeep(col_low_hi,[r+1 c+1]),'m');
ccv=dyaddown(wkeep(col_hi_low,[r+1 c+1]),'m');
ccd=dyaddown(wkeep(col_hi_hi,[r+1 c+1]),'m');
assert_checkalmostequal ( cA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( cH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( cV , ccv , %eps,  %eps*1000 );
assert_checkalmostequal ( cD , ccd , %eps, %eps*1000 );
assert_checkalmostequal ( ccA , cca , %eps, %eps*1000 );
assert_checkalmostequal ( ccH , cch , %eps, %eps*1000 );
assert_checkalmostequal ( ccV , ccv , %eps, %eps*1000 );
assert_checkalmostequal ( ccD , ccd , %eps, %eps*1000 );
clear cA;
clear cH;
clear cV;
clear cD;

clear m_ex;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear cca;
clear cch;
clear ccv;
clear ccd;
clear r;
clear c;
clear rex;
clear cex;
clear rrex;
clear ccex;
clear dd1;
clear x1;
clear x2;
clear d1;
clear d2;

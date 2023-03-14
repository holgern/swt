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
// detcoef
level=10;
[c,l]=wavedec(s1,level,'sym10');
cA=list();
cD=list();
[cA(1),cD(1)]=dwt(s1,'sym10');
for (i=2:level)
	[cA(i),cD(i)]=dwt(cA(i-1),'sym10');
end;

cddetMax=detcoef(c,l);
cdet=list();
for (i=1:level)
cdet(i)=detcoef(c,l,i);
end;


assert_checkalmostequal ( cddetMax , cD(level) , %eps );
for (i=1:level)
	assert_checkalmostequal ( cdet(i) , cD(i) , %eps );
end;

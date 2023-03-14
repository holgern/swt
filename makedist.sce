//  Copyright (C) 2020 Stéphane Mottelet.
//
//  SCIIPOPT is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Scilab libraries and it can be 
//  used within the Scilab environment.
//
//  SCIIPOPT is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.

path_makedist = get_absolute_file_path("makedist.sce");
cd(path_makedist);

deletefile(".gitignore");
[sci,v]=getversion()
http_get("https://atoms.scilab.org/toolboxes/swt/"+mgetl("VERSION")+"/DESCRIPTION","DESCRIPTION")
if getos() <> "Windows" then
    filename = "swt_"+mgetl("VERSION")+".bin."+v(2)+"."+getos()+".tar.gz"
    cd ..
    movefile("swt/makedist.sce","makedist.sce")
    unix_g("tar cvzf "+filename+" swt")
else
    filename = "swt_"+mgetl("VERSION")+".bin."+v(2)+"."+getos()+".zip"
    cd ..
    movefile("swt/makedist.sce","makedist.sce")
    host(SCI+"/tools/zip/zip.exe -r "+filename+" sci-ipopt")
end
movefile("makedist.sce","swt/makedist.sce")
cd(path_makedist);

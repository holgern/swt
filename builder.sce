// Copyright (C) 2008 - INRIA
// Copyright (C) 2009-2011 - DIGITEO

// This file is released under the 3-clause BSD license. See COPYING-BSD.

mode(-1);
function builder_main()

  TOOLBOX_NAME  = "swt";
  TOOLBOX_TITLE = "wavelet toolbox";
  toolbox_dir   = get_absolute_file_path("builder.sce");

  // Check Scilab's version
  // =============================================================================

  try
	  v = getversion("scilab");
  catch
	  error(gettext("Scilab 5 or more is required."));
  end
  if v(1)==5 & v(2) < 3 then
    // new API in scilab 5.3
    error(gettext('Scilab 5.3 or more is required.'));
  end

  // =============================================================================
  // Uncomment this line to make a debug version of the Toolbox
  // setenv("DEBUG_SCILAB_DYNAMIC_LINK","YES")
  //setenv("__USE_DEPRECATED_STACK_FUNCTIONS__","YES")

  // Action
  // =============================================================================

  tbx_builder_macros(toolbox_dir);
  tbx_builder_src(toolbox_dir);
  tbx_builder_gateway(toolbox_dir);
  if (getscilabmode() ~= "NWNI") then
	tbx_builder_help(toolbox_dir);
  end;
  tbx_build_loader(toolbox_dir);
  tbx_build_cleaner(toolbox_dir);

endfunction 
// =============================================================================
builder_main();
clear builder_main;
// =============================================================================


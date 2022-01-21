#!/usr/bin/env python
#
# Script to build the files found in this directory.
#------------------------------------------------------------------------------
import os, re, string, subprocess, sys, importlib
Import('env')
sys.path.append(os.getenv("MUSE_WORK_DIR")+'/site_scons')
from SCons.Script import *

print("mu2eii/SConscript:emoe: PWD:"+os.getenv("PWD"))

x = subprocess.call(os.getenv("MUSE_WORK_DIR")+'/mu2eii/scripts/build_config_muse',shell=True)

stntuple_env = env.Clone()

stntuple_env['CPPPATH' ].append('-I'+os.environ['MUSE_WORK_DIR']+'/build/include');
stntuple_env['CXXFLAGS'].append('-I'+os.environ['MUSE_WORK_DIR']+'/build/include');
stntuple_env.Append(FORTRANFLAGS = ['-I'+os.environ['MUSE_WORK_DIR']+'/build/include']);
#------------------------------------------------------------------------------
# done
#------------------------------------------------------------------------------
exec(open(os.environ['MUSE_WORK_DIR']+"/site_scons/mu2eii_site_init.py").read())

from stntuple_helper    import *

stntuple_env.Append(BUILDERS = {'StntupleCodegen'  : stntuple_codegen })
stntuple_env.Append(BUILDERS = {'StntupleRootCint' : stntuple_rootcint})

Export('stntuple_env')
Export('stntuple_helper')

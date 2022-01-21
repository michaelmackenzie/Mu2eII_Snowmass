import os, re, string, subprocess
#------------------------------------------------------------------------------
def stntuple_codegen_script(source, target, env, for_signature):
    # print("\n>>> mod_codegen called: pwd = "+os.environ['PWD'])
    # print(">>> source = ",len(source),str(source[0]))
    # print(">>> target = ",len(target),str(target[0])+'\n')
    # print(">>> env[\'MUSE_WORK_DIR\']:"+os.environ['MUSE_WORK_DIR'])
    
    tmpdir = os.path.dirname(str(target[0]));

    cmd    = 'if [ ! -d '+tmpdir+' ] ; then mkdir -p '+tmpdir+'; fi ; ';
    cmd    = cmd+'source '+str(source[0]) + ' ' +str(target[0]);

    # print( ">>> cmd = %s"%cmd)
    return cmd

stntuple_codegen = Builder(generator     = stntuple_codegen_script,
                           single_source = 0,
                           suffix        = '.cc',
                           src_suffix    = '.sh')

#------------------------------------------------------------------------------
def stntuple_gen_rootcint(source, target, env, for_signature):
    # print("\n>>> stntuple_gen_rootcint called:")
    # print(">>> source = ",len(source),str(source[0]), str(source[1]))
    # print(">>> target = ",len(target),str(target[0]))

    class_include = str(source[1]);
    linkdef       = str(source[0]);

    # print ("[stntuple_gen_rootcint] class_include = %s"%class_include)
    # print ("[stntuple_gen_rootcint] linkdef       = %s"%linkdef)
    #------------------------------------------------------------------------------
    # building in a satellite release requires both env vars 
    #------------------------------------------------------------------------------
    base_dir = os.environ.get('MUSE_WORK_DIR');
    includes = '-I'+base_dir+'/build/include';

    includes = includes + " -I"+os.environ['MUSE_WORK_DIR'];
    includes = includes + " -I"+os.environ['ART_INC'     ];
    includes = includes + " -I"+os.environ['ART_ROOT_IO_INC'];
    includes = includes + " -I"+os.environ['BTRK_INC'    ];
    includes = includes + " -I"+os.environ['CETLIB_INC'  ];
    includes = includes + " -I"+os.environ['CETLIB_EXCEPT_INC'];
    includes = includes + " -I"+os.environ['CANVAS_INC'  ];
    includes = includes + " -I"+os.environ['FHICLCPP_INC'];
    includes = includes + " -I"+os.environ.get('HEP_CONCURRENCY_INC');
    includes = includes + " -I"+os.environ['CLHEP_INC'   ];
    includes = includes + " -I"+os.environ['BOOST_INC'   ];
    includes = includes + " -I"+os.environ.get('TBB_INC');

    dict     = str(target[0]);
    tmp_lib_dir = os.path.dirname(dict);
    
    pcm_file = dict.replace(".cxx","_rdict.pcm");

    # print( "[stntuple_gen_rootcint] dict:"+dict + "   pcm_file:"+pcm_file)
    #------------------------------------------------------------------------------
    # if building a satellite release, BUILD_BASE points to a remote directory
    # need to use MU2E_SATELLITE_RELEASE
    #------------------------------------------------------------------------------
    dir = os.environ.get('MUSE_BUILD_DIR')+'/'+'mu2eii';

    lib_dir = dir+"/lib";
    cmd     = 'if [ ! -d '+tmp_lib_dir+' ] ; then mkdir -p '+tmp_lib_dir+'; fi;';
#------------------------------------------------------------------------------
# export $HOME is, for an unknown reason, needed to build in a DOCKER container
# the same is needed to use genreflex - that change goes into SConstruct
# without that, ROOT doesn't find the home directory.
# The mechanism of failure is not understood.
#------------------------------------------------------------------------------
    cmd = cmd+" export HOME="+os.environ["HOME"]+";";
    cmd = cmd+" rootcint -f "+dict+" -c  -D_CODEGEN_ -DMU2E "+includes+" "+class_include+" "+linkdef+"; ";
    cmd = cmd+'if [ ! -d '+lib_dir+' ] ; then mkdir '+lib_dir+' ; fi ; ';
    cmd = cmd+"mv "+pcm_file+" "+lib_dir+'/.'; 
    # print( "[stntuple_gen_rootcint] >>> cmd = %s"%cmd)
    return cmd

stntuple_rootcint = Builder(generator     = stntuple_gen_rootcint,
                            single_source = 0,
                            suffix        = '.o',
                            src_suffix    = '.h')

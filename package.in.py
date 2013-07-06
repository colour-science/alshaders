#! /usr/bin/python

import sys
import os
import errno
import shutil
import glob
import tarfile
import platform

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

MAJOR_VERSION = @ALS_MAJOR_VERSION@
MINOR_VERSION = @ALS_MINOR_VERSION@
PATCH_VERSION = @ALS_PATCH_VERSION@

subdirs = [
    'alBlackbody',
    'alColorSpace',
    'alCombine',
    'alCurvature',
	'alHair',
	'alInputVector',
    'alLayer',
    'alNoise',
    'alPattern',
    'alPhotometric',
    'alRemap',
    'alSurface',
    'test',
    'common'
] 

# pattern for files to include in the source distribution
name_src = 'alShaders-src-%d.%d.%d' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION)
ptrn_src_src = ['*.cpp', '*.h', '*.mtd', '*.txt', '*.py', '*.ass', '*.ass.gz', '*.nk']
ptrn_src_bin = []
files_src = ['BUILD_INSTRUCTIONS.txt', 'CMakeLists.txt', 'package.in.py', 'test.in.py', 'README', 'TODO.txt']

# patterns for files to include in the osx binary distribution
name_osx = 'alShaders-osx-%d.%d.%d' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION)
ptrn_osx_src = ['*.mtd', '*.py']
ptrn_osx_bin  = ['*.dylib']
files_osx = ['BUILD_INSTRUCTIONS.txt', 'README']

# patterns for files to include in the windows binary distribution
name_win = 'alShaders-win-%d.%d.%d' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION)
ptrn_win_src = ['*.mtd', '*.py']
ptrn_win_bin  = ['*.dll']
files_win = ['BUILD_INSTRUCTIONS.txt', 'README']

# patterns for files to include in the linux binary distribution
name_linux = 'alShaders-linux-%d.%d.%d' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION)
ptrn_linux_src = ['*.mtd', '*.py']
ptrn_linux_bin  = ['*.so']
files_linux = ['BUILD_INSTRUCTIONS.txt', 'README']

def copyPatternsToDistDir(subDirs, subDirPrefix, filePatterns, distDir):
    for dir in subDirs:
        destdir = os.path.join(distDir, dir)
        mkdir_p(destdir)
        subdir = os.path.join(subDirPrefix, dir)
        for ptrn in filePatterns:
            for fn in glob.iglob(os.path.join(subdir, ptrn)):
                shutil.copy(fn, destdir)
                        
def copyFilesToDistDir(files, distDir):
    for fn in files:
        shutil.copy(fn, distDir)
        
def createArchive(distDir, name):
    f = tarfile.open(os.path.join('..', '%s.tar.gz' % name), 'w:gz')
    f.add(distDir, arcname = name)
    f.close()
    
def createDistribution(name, ptrn_src, ptrn_bin, files):
    distdir = 'build/%s' % name
    shutil.rmtree(distdir, ignore_errors=True)
    os.mkdir(distdir)
    if len(ptrn_src):
        copyPatternsToDistDir(subdirs, '.', ptrn_src, distdir)
    if len(ptrn_bin):
        copyPatternsToDistDir(subdirs, 'build', ptrn_bin, distdir)
    copyFilesToDistDir(files, distdir)  
    createArchive(distdir, name)

# source distribution
createDistribution(name_src, ptrn_src_src, ptrn_src_bin, files_src)
if platform.system() == "Darwin":
    # OS X distribution
    createDistribution(name_osx, ptrn_osx_src, ptrn_osx_bin, files_osx)
elif platform.system() == "Windows":
    # Windows
    createDistribution(name_win, ptrn_win_src, ptrn_win_bin, files_win)
elif platform.system() == "Linux":
    # Linux
    createDistribution(name_linux, ptrn_linux_src, ptrn_linux_bin, files_linux)
else:
    print 'Warning: unknown system "%s", not creating binary package' % platform.system()

    
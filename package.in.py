#! /usr/bin/python

import sys
import os
import errno
import shutil
import glob
import tarfile
import zipfile
import platform

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

MAJOR_VERSION = '@ALS_MAJOR_VERSION@'
MINOR_VERSION = '@ALS_MINOR_VERSION@'
PATCH_VERSION = '@ALS_PATCH_VERSION@'

ARNOLD_VERSION = '@ARNOLD_VERSION@'

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
    'common'
] 

# pattern for files to include in the source distribution
name_src = 'alShaders-src-%s.%s.%s-ai-%s' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION, ARNOLD_VERSION)
ptrn_src_src = ['*.cpp', '*.h', '*.txt', '*.py']
ptrn_src_bin = []
files_src = ['BUILD_INSTRUCTIONS.txt', 'CMakeLists.txt', 'package.in.py', 'test.in.py', 'README', 'TODO.txt', 'local.cmake.win']

# patterns for files to include in the osx binary distribution
name_osx = 'alShaders-osx-%s.%s.%s-ai-%s' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION, ARNOLD_VERSION)
ptrn_osx_src = ['*.py']
ptrn_osx_bin  = ['*.dylib', '*.mtd', '*.py']
files_osx = ['BUILD_INSTRUCTIONS.txt', 'README']

# patterns for files to include in the windows binary distribution
name_win = 'alShaders-win-%s.%s.%s-ai-%s' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION, ARNOLD_VERSION)
ptrn_win_src = ['*.py']
ptrn_win_bin  = ['*.dll', '*.mtd', '*.py', '*.spdl']
files_win = ['BUILD_INSTRUCTIONS.txt', 'README']

# patterns for files to include in the linux binary distribution
name_linux = 'alShaders-linux-%s.%s.%s-ai-%s' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION, ARNOLD_VERSION)
ptrn_linux_src = ['*.py']
ptrn_linux_bin  = ['*.so', '*.mtd', '*.py', '*.spdl']
files_linux = ['BUILD_INSTRUCTIONS.txt', 'README']

def copyPatternsToDistDir(subDirs, subDirPrefix, filePatterns, distDir, isSrc=False):
    for dir in subDirs:
        if isSrc:
            destdir = os.path.join(distDir, dir)
        else:
            destdir = distDir
        mkdir_p(destdir)
        subdir = os.path.join(subDirPrefix, dir)
        for ptrn in filePatterns:
            sdir = subdir
            if ptrn == '*.dll' and platform.system() == "Windows":
                sdir = os.path.join(sdir, 'Release')
            for fn in glob.iglob(os.path.join(sdir, ptrn)):
                shutil.copy(fn, destdir)
                        
def copyFilesToDistDir(files, distDir):
    for fn in files:
        shutil.copy(fn, distDir)
        
def createArchive(distDir, name, isSrc=False):
    if not isSrc and platform.system() == "Windows":
        f = zipfile.ZipFile(os.path.join('..', '%s.zip' % name), 'w')
        for fn in glob.iglob(os.path.join(distDir, '*')):
            f.write(fn, arcname=os.path.join(name, os.path.basename(fn)))
        f.close()
    else:
        f = tarfile.open(os.path.join('..', '%s.tar.gz' % name), 'w:gz')
        f.add(distDir, arcname = name)
        f.close()
    
def createDistribution(name, ptrn_src, ptrn_bin, files, isSrc=False):
    distdir = 'build/%s' % name
    shutil.rmtree(distdir, ignore_errors=True)
    os.mkdir(distdir)
    if len(ptrn_src):
        copyPatternsToDistDir(subdirs, '.', ptrn_src, distdir, isSrc)
    if len(ptrn_bin):
        copyPatternsToDistDir(subdirs, 'build', ptrn_bin, distdir, isSrc)
    copyFilesToDistDir(files, distdir)  
    createArchive(distdir, name, isSrc)

# source distribution
createDistribution(name_src, ptrn_src_src, ptrn_src_bin, files_src, True)
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

    
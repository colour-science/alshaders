#! /usr/bin/python

import sys
import os
import errno
import shutil
import glob
import tarfile

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

MAJOR_VERSION = 0
MINOR_VERSION = 2
PATCH_VERSION = 1

subdirs = [
    'alBlackbody', 
    'alCombine',
    'alCurvature',
    'alLayer',
    'alNoise',
    'alPattern',
    'alPhotometric',
    'alRemap',
    'alSurface',
    'alInputVector',
    'alHair',
    'common'
] 

# pattern for files to include in the source distribution
name_src = 'alShaders-src-%d.%d.%d' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION)
ptrn_src_src = ['*.cpp', '*.h', '*.mtd', '*.mel', '*.txt']
ptrn_src_bin = []
files_src = ['BUILD_INSTRUCTIONS.txt', 'CMakeLists.txt', 'package.py', 'README', 'TODO.txt']

# patterns for files to include in the osx binary distribution
name_osx = 'alShaders-osx-%d.%d.%d' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION)
ptrn_osx_src = ['*.mtd', '*.mel']
ptrn_osx_bin  = ['*.dylib']
files_osx = ['BUILD_INSTRUCTIONS.txt', 'README']

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
# OS X distribution
createDistribution(name_osx, ptrn_osx_src, ptrn_osx_bin, files_osx)

    
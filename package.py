#! /usr/bin/python

import sys
import os
import shutil

MAJOR_VERSION = 0
MINOR_VERSION = 1
PATCH_VERSION = 0

versionedname = 'alShaders-%d.%d.%d' % (MAJOR_VERSION,MINOR_VERSION,PATCH_VERSION)
distdir = 'build/%s' % versionedname

shutil.rmtree(distdir, ignore_errors=True)
os.mkdir(distdir)

subdirs = os.walk('.').next()[1]

# pattern for files to include in the source distribution
filepatterns = ['*.cpp', '*.h', '*.mtd', '*.mel', '*.txt']

for dir in subdirs:
	if dir[0] == '.' or dir == 'build':
		continue
	
	destdir = os.path.join(distdir, dir)
	os.mkdir(destdir)
	for ptrn in filepatterns:
		cpcmd = 'cp %s %s' % (os.path.join(dir,ptrn), destdir)
		#print cpcmd
		os.system(cpcmd)
		
# now copy the top-level files
shutil.copy('BUILD_INSTRUCTIONS.txt', distdir)
shutil.copy('CMakeLists.txt', distdir)
shutil.copy('package.py', distdir)
shutil.copy('README', distdir)
shutil.copy('TODO.txt', distdir)

# now create the archive
os.system('mv %s .' % distdir)
tarcmd = 'tar cf %s.tar %s' % (versionedname, versionedname)
gzcmd = 'gzip -9 %s.tar' % versionedname
os.system(tarcmd)
os.system(gzcmd)
os.system('rm -r %s' % versionedname)
os.system('mv %s.tar.gz ..' % versionedname)
	
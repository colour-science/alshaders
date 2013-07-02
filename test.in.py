# test.py
# run the test suite of renders and report timings etc
import subprocess
import time
import os

VERSION = "@ALS_VERSION@"

class Timer:    
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start

# names of the tests
tests = {
			'als_bplastic_rough':1,
			'als_lyr_glass_goldleaf':0
		}

# first make the output directory where the results will be stored
output_dir = os.path.join(os.getcwd(), "test/output/%s" % VERSION)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)


for test_name,opaque in tests.items():
	header = open("test/test_header.ass", "r").read()
	output = 'include "test/test_%s.ass"\n' % test_name
	output += header
	test_tmp = open('test/test_tmp.ass', 'w')
	test_tmp.write(output)
	test_tmp.close()

	cmd = 'kick -dp -dw -set ringShape.opaque %d -set shellShape.opaque %d -set driver_exr_beauty.filename "test/output/%s/%s.exr" test/test_tmp.ass' % (opaque, opaque, VERSION, test_name)

	with Timer() as t:
		proc = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		rc = proc.wait()
	
	if rc == 0:
		print 'Test %s completed successfully in %.02f seconds' % (test_name, t.interval)
	else:
		print 'Test %s FAILED after %.02f seconds' % (test_name, t.interval)
		for line in proc.stdout:
			print line[:-1]

# clean up temporary files
os.remove('test/test_tmp.ass')
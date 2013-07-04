# test.py
# run the test suite of renders and report timings etc
import subprocess
import time
import os
import sys

VERSION = "@ALS_VERSION@"

class Timer:    
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start

# names of the tests and whether the surfaces should be opaque or not
tests = {
			'als_bplastic_shiny':1,
			'als_bplastic_rough':1,
			'als_goldleaf':1,
			'als_glass':0,
			'als_lyr_glass_goldleaf':0,
			'als_merlot':0,
			'als_ccover_bump':1,
			'als_metallicpaint':1,
			'als_skin':1,
			'als_jade':0
		}

# first make the output directory where the results will be stored
output_dir = os.path.join(os.getcwd(), "test/output/%s" % VERSION)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

# loop over the tests, combine the asses and render them
test_results = {}
for test_name,opaque in tests.items():
	header = open("test/test_header.ass", "r").read()
	output = 'include "test/test_%s.ass"\n' % test_name
	output += header
	test_tmp = open('test/test_tmp.ass', 'w')
	test_tmp.write(output)
	test_tmp.close()

	log = open('%s/log_%s.txt' % (output_dir, test_name), 'w')

	cmd = 'kick -v 6 -t 2 -dp -dw -set ringShape.opaque %d -set shellShape.opaque %d -set driver_exr_beauty.filename "test/output/%s/%s.exr" test/test_tmp.ass' % (opaque, opaque, VERSION, test_name)

	sys.stdout.write('%s...' % test_name)
	sys.stdout.flush()

	with Timer() as t:
		proc = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
		rc = proc.wait()
	
	if rc == 0:
		print 'OK in %.02f seconds' % t.interval
		result = t.interval
	else:
		print 'FAILED after %.02f seconds' % t.interval
		result = 0

	for line in proc.stdout:
		log.write(line)
	log.close()

	test_results[test_name] = result

# write out the results to a single file
result_total = 0.0
results_file = open('test/output/%s/results_%s.dat' % (VERSION,VERSION), 'w')
for test_name, result in test_results.items():
	results_file.write('%s %f\n' % (test_name, result))
	result_total += result
results_file.write('\nTOTAL: %f\n' % result_total)
results_file.close()

# clean up temporary files
os.remove('test/test_tmp.ass')
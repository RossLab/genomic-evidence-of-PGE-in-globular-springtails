#!/usr/bin/env python

import sys
triallic = 0

for line in sys.stdin:
    line_elements = line.split('\t')
    scf = line_elements[0]
    pos = line_elements[1]
    cov_supports = [int(cov) for cov in line_elements[3].split(':')[0:4]]
    cov_supports.sort(reverse = True)
    if cov_supports[2] > 0:
        triallic += 1
    if cov_supports[1] > 0 and cov_supports[2] == 0:
        sys.stdout.write('{}\t{}\t{}\t{}\n'.format(scf, pos, cov_supports[0], cov_supports[1]))

sys.stderr.write('discarded in total ' + str(triallic) + ' sites\n')
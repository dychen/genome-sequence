"""
Requires the following files in the same directory:
    T03-1_S10_L001_R1_001.fastq
    T03-1_S10_L001_R2_001.fastq
Outputs the following files:
    Target_1_Mismatch_1.txt
    Target_1_Mismatch_2.txt
    ...
    Target_4_Mismatch_4.txt
"""

FILE_1 = 'T03-1_S10_L001_R1_001.fastq'
FILE_2 = 'T03-1_S10_L001_R2_001.fastq'

ADAPTOR = 'AGATCGGAAGAGCGTCGTGTAGGGAAA'
TARGET_1 = 'AGCTGACGTTTGTACTCCAGCG'
TARGET_2 = 'AGCTGACGTTTGTACTCCAGCGTCTCATCTTTATGCGTCAGCAGAGATTTCTGCT'
TARGET_3 = 'AGCAGAAATCTCTGCTGACGCATAAAGATGAGA'
TARGET_4 = 'AGCAGAAATCTCTGCTGACGCATAAAGATGAGACGCTGGAGTACAAACGTCAGCT'
TARGET_MAP = {
    1: TARGET_1,
    2: TARGET_2,
    3: TARGET_3,
    4: TARGET_4,
}


def match_sequence(sequence, target, distance=0):
    target_subsequence = sequence[-len(target):]
    # Edge case: Adaptor subsquence is found too close to the beginning of the
    #            sequence. Length of target subsequence is less then the length
    #            of target.
    if len(target_subsequence) == len(target):
        err_ct = 0
        for i in xrange(len(target_subsequence)):
            if target_subsequence[i] != target[i]:
                err_ct += 1
        if err_ct == distance:
            return target_subsequence
        else:
            return None
    return None

def adaptor_match(sequence):
    try:
        return sequence.index(ADAPTOR)
    except ValueError: # No match
        return -1


def analyze_sequence(header, sequence, counts_dict, subsequences_dict,
                     debug=False):
    adaptor_index = adaptor_match(sequence)
    if adaptor_index > -1:
        subsequence = sequence[:adaptor_index]
        for target_idx, target in enumerate([TARGET_1, TARGET_2, TARGET_3,
                                             TARGET_4]):
            for distance in [0, 1, 2, 3, 4]:
                subsequence_match = match_sequence(subsequence, target,
                                                   distance)
                if subsequence_match:
                    if debug:
                        print '====='
                        print header
                        print 'Target', target_idx+1
                        print 'Differences', distance
                        print subsequence_match
                        print target
                        print '====='
                    if counts_dict:
                        counts_dict[target_idx+1][distance] += 1
                    if subsequences_dict:
                        if distance > 0:
                            subsequences_dict[target_idx+1][distance].append(
                                subsequence_match
                            )


def print_counts_dict(counts_dict):
    for target_idx in [1, 2, 3, 4]:
        for distance in [0, 1, 2, 3, 4]:
            print ('Target: %d, Differences: %d, Count: %d'
                   % (target_idx, distance, counts_dict[target_idx][distance]))


def write_subsequences_dict(subsequences_dict):
    for target_idx in [1, 2, 3, 4]:
        for distance in [1, 2, 3, 4]:
            if len(subsequences_dict[target_idx][distance]) > 0:
                filename = ('Target_%d_Mismatch_%d.txt'
                            % (target_idx, distance))
                with open(filename, 'w') as f:
                    f.write('%s\n' % TARGET_MAP[target_idx])
                    f.write('=====\n')
                    for subsequence in subsequences_dict[target_idx][distance]:
                        f.write('%s\n' % subsequence)


def analyze(filename):
    """
    Example format (group lines in groups of 3):
        @M05354:6:000000000-B87B4:1:1102:9664:2850 2:N:0:10
        AGCTGACGTTTGTACTCCAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG...
        +
        @M05354:6:000000000-B87B4:1:1102:21618:2852 2:N:0:10
        AGCTGACGTTTGTACTCCAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG...
        +
    Counts dictionary format:
        {
            [Target index]: {
                [Distance number (# edits)]: [int],
                ...
            },
            ...
        }
    """
    counts_dict = {
        1: { 0: 0, 1: 0, 2: 0, 3: 0, 4: 0 },
        2: { 0: 0, 1: 0, 2: 0, 3: 0, 4: 0 },
        3: { 0: 0, 1: 0, 2: 0, 3: 0, 4: 0 },
        4: { 0: 0, 1: 0, 2: 0, 3: 0, 4: 0 },
    }
    # Store matching subsequences for all targets with edit distances 1-4
    subsequences_dict = {
        1: { 1: [], 2: [], 3: [], 4: [] },
        2: { 1: [], 2: [], 3: [], 4: [] },
        3: { 1: [], 2: [], 3: [], 4: [] },
        4: { 1: [], 2: [], 3: [], 4: [] },
    }
    ct_total = 0
    with open(filename, 'r') as f:
        header, sequence = '', ''
        for line in f:
            # Expected: line of length 2 with characters: '+'
            if len(line) < 10 and '+' in line:
                analyze_sequence(header, sequence, counts_dict,
                                 subsequences_dict)
                ct_total += 1
            # Expected: line of length 301 with characters: 'ATCG...'
            elif len(line) > 200:
            # Expected: @M05354:6:000000000-B87B4:1:1102:21618:2852 2:N:0:10
                sequence = line.strip()
            else:
                header = line.strip()
    print_counts_dict(counts_dict)
    print 'Total: %d' % ct_total
    write_subsequences_dict(subsequences_dict)


if __name__=='__main__':
    #analyze(FILE_1)
    analyze(FILE_2)

from HMM import *

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub)

if __name__ == '__main__':
    model1 = create_restriction_model(site="ACGT", p_site=0.01)
    model2 = create_restriction_model(site="ACGT", p_site=0.001)
    test_sequence = 'gataggattatcattcataagtttcagagcaatgtccttattctggaacttggatttatggctcttttggtttaatttcgcctgattcttgatctcctttagcttctcgacgtgggcctttttcttgccatatggatccgctgcacggtcctgttccctagcatgtacgtgagcgtatttccttttaaaccacgacgctttgtcttcattcaacgtttcccattgtttttttctactattgctttgctgtgggaaaaacttatcgaaagatgacgactttttcttaattctcgttttaagagcttggtgagcgctaggagtcactgccag'
    test_sequence = test_sequence.upper()
    res1 = model1.viterbi(test_sequence)
    res2 = model2.viterbi(test_sequence)
    # check if model with p_site=0.01 found all sites in right place
    true_sites = list(find_all(test_sequence, 'ACGT'))
    assert list(find_all(res1, 'SSSS')) == true_sites
    print('Model with p_site=0.01 found all sites in right place!')
    # check if model did not mark wrong sites
    assert len(true_sites) * 4 == res1.count('S')
    print('Model with p_site=0.01 did not mark wrong nucleotides as site!')
    print()
    # mark all true sites with `*`
    for i in true_sites.copy():
        true_sites.extend(list(range(i + 1, i + 4)))
    sites_marks = ''.join(' ' if i not in true_sites else '*' for i in range(len(test_sequence)))
    print('Model with p_site = 0.001: ', res2)
    print('Model with p_site = 0,01:  ',res1)
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>', test_sequence)
    print('         SITES:            ', sites_marks)
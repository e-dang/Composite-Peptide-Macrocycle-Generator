import itertools

length = 3
stereo = 'D'
aa = ['Gly', 'Ala', 'Leu', 'Met', 'Phe', 'Trp', 'Lys', 'Gln', 'Glu', 'Ser', 'Pro', 'Val', 'Ile', 'Cys', 'Tyr', 'His', 'Arg', 'Asn', 'Asp', 'Thr']

with open('~/Documents/UCLA_Research/aa_seq/aa' + str(length) + stereo + '.txt', 'w') as f:
    ls = []
    for combo in itertools.permutations(aa, length):
        if 'Cys' in combo or 'His' in combo or 'Lys' in combo or 'Ser' in combo or 'Trp' in combo or 'Asn' in combo or 'Tyr' in combo:
            ls.append(combo)

    f.write(str(len(ls)) + '\n')
    for item in ls:
        rep = ''
        for amino in item:
            if amino == 'Gly':
                rep += amino + '-'
            else:
                rep += stereo + '-' + amino + '-'
        rep = rep[0:-1]
        f.write(rep)
        f.write('\n')

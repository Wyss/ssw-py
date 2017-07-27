def printer(res, ref, read):
    print("optimal_score: %d\tsub-optimal_score: %d\t" %
        (res['optimal_score'], res["sub-optimal_score"]))
    target_begin = res['target_begin']
    if target_begin >= 0:
        print("target_begin: %d\t" % (target_begin), end='')
    print("target_end: %d\t" % (res['target_end']) )
    query_begin = res['query_begin']
    if query_begin >= 0:
        print("query_begin: %d\t" % (query_begin), end='')
    print("query_end: %d\n" % (res['query_end']))
    cigar = res['CIGAR']

    if cigar is not None:
        e = 0
        left = 0
        pb = query_begin
        qb = target_begin
        STEP2 = False
        STEP3 = False
        STEPEND = False
        lim = len(cigar)
        while e < lim or left > 0:
            q = qb
            p = pb
            count = 0
            print("Target: %8d    " % (q), end='')
            for c in range(e, lim, 2):
                length = int(cigar[c])
                letter = cigar[c + 1]
                L1 = left if (count == 0 and left > 0) else length
                for j in range(L1):
                    if letter == 'I':
                        print("-", end='')
                    else:
                        print(ref[q], end='')
                        q += 1
                    count += 1
                    if count == 60:
                        STEP2 = True
                        break # goto STEP 2
                if STEP2:
                    STEP2 = False
                    break
            # end for
            # STEP 2
            print("    %d\n                    " % (q - 1), end='')
            q = qb
            count = 0
            for c in range(e, lim, 2):
                length = int(cigar[c])
                letter = cigar[c + 1]
                L1 = left if (count == 0 and left > 0) else length
                for j in range(L1):
                    if letter == 'M':
                        if ref[q] == read[p]:
                            print('|', end='')
                        else:
                            print("*", end='')
                        q += 1
                        p += 1
                    else:
                        print('*', end='')
                        if letter == 'I':
                            p += 1
                        else:
                            q += 1
                    count += 1
                    if count == 60:
                        qb = q
                        STEP3 = True
                        break # goto STEP 3
                if STEP3:
                    STEP3 = False
                    break
            # end for
            # STEP 3
            p = pb
            print("\nQuery:  %8d    " % (p), end='')
            count = 0
            for c in range(e, lim, 2):
                length = int(cigar[c])
                letter = cigar[c + 1]
                L1 = left if (count == 0 and left > 0) else length
                for j in range(L1):
                    if letter == 'D':
                        print('-', end='')
                    else:
                        print(read[p], end='')
                        p += 1
                    count += 1
                    if count == 60:
                        pb = p
                        left = L1 - j - 1
                        e = c + 2 if left == 0 else c
                        STEPEND = True
                        break # goto STEPEND
                if STEPEND:
                    STEPEND = False
                    break
            e = c + 2
            left = 0
            # STEPEND
            print("    %d\n" % (p - 1))
        # end while
    #end if
# end def
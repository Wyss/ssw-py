from ssw import SSW
from ssw.printer import printer
import imp
print(imp.find_module('ssw'))

if __name__ == '__main__':
    a = SSW()
    a.setRead(b"ACGT")

    print("1. test exact")
    ref = b"TTTTACGTCCCCC"
    a.setReference(ref)
    res = a.align()
    a.printResult(res)
    print("2. test deletion")
    ref = b"TTTTACAGTCCCCC"
    a.setReference(ref)
    res = a.align()
    a.printResult(res)
    print("3a. test insertion with gap_open penalty=3")
    ref = b"TTTTACTCCCCC"
    a.setReference(ref)
    res = a.align()
    a.printResult(res)
    print("3b. test insertion with gap_open penalty=0")
    res = a.align(gap_open=0)
    a.printResult(res)
    print(")))))))))))))))))))")
    printer(res, ref.decode('utf8'), "ACGT")

    print("4. start_idx test")
    b = SSW()
    b.setRead("ACTG")
    b.setReference("ACTCACTG")
    res = b.align(start_idx=4)
    b.printResult(res, start_idx=4)

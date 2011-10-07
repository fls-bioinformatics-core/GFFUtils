import sys
if __name__ == "__main__":
    f1 = open(sys.argv[1],'r')
    f2 = open(sys.argv[2],'r')
    line1 = f1.readline()
    line2 = f2.readline()
    if line1 == line2:
        print "First lines match"
    else:
        print "First lines don't match"
        if len(line1) != len(line2):
            print "Lines have different lengths"
            if len(line1) > len(line2):
                print "Line 1 is longer"
            else:
                print "Line 2 is longer"
        print "%s" % line1
        print "%s" % line2
        for i in range(min(len(line1),len(line2))):
            c1 = line1[i]
            c2 = line2[i]
            if c1 == c2:
                print ' ',
            else:
                print '^',
    f1.close()
    f2.close()

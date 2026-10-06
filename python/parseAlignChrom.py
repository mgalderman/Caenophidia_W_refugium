import sys

chromlist = []

for line in open(sys.argv[1], 'r'):
        chrom = line.split()[0]
        chromlist.append(chrom)

#for chrom in chromlist:
#       print "gathering aligments for"+str(chrom)

current = "NA"

# for line in open(sys.argv[2], 'r'):
#       if len(line.strip()) > 0:
#               check = line.split()[0]
#               if check[0].isdigit():
# #                     print line.rstrip()
#                       current = str(line.split()[4])
#                       if current in chromlist:
#                               print line.rstrip()

for line in open(sys.argv[2], 'r'):
        if len(line.strip()) > 0:
                check = line.split()[0]
                if check[0].isdigit():
                        elements = line.split()
                        if len(elements) >= 5:
#                               print line.rstrip()
                                current = str(line.split()[4])
                                if current in chromlist:
                                        print line.rstrip()
                        else:
                                if str(current) in chromlist:
                                        print line.rstrip()
                else:
                        if str(current) in chromlist:
                                print line.rstrip()
        else:
                if str(current) in chromlist:
                        print line.rstrip()
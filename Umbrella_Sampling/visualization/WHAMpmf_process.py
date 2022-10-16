from sys import argv

def process_pmf(fnamein, fnameout, *args):
    fhandle = open(fnamein,'r')
    incontent = fhandle.readlines()
    fhandle.close()
    cleancontent = []
    for i in incontent:
        temp = i.split()
        if len(temp) == 0:
            continue
        tempout = []
        try:

            for j in temp:
                    tempout.append(j)
            cleancontent.append(tempout)
        except:
            continue
    fouthandle = open(fnameout, 'w')
    for i in cleancontent:
        i.append('\n')
        temp = "\t".join(i)
        fouthandle.write(temp)
    fouthandle.close()


process_pmf(*argv[1:])
#process_pmf('out2.pmf', 'out_processed.pmf')



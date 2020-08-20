"""Runs various checks to **TRY** to ensure this file won't crash BUFKIT. I'm not
sure we'll ever be 100% sure BUFKIT won't die, however!
"""
import matplotlib.pyplot as plt
import numpy as np

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))

filename = './soundings/ACARSvapor_MDW.buf'
#filename = './soundings/rap_kmdw.buf'
with open(filename) as f: data = f.readlines()

ret = False
indices = [index for index, value in enumerate(data) if value == \
           'PRES TMPC TMWC DWPC THTE DRCT SKNT OMEG\n']
sum = 0
max_p = -9999
if len(indices) >= 5:
    for i in range(len(indices)-1):
    #for i in range(10):
        start = indices[i]+2

        # Figure out the end of this section
        pres = []
        hght = []
        tmpc = []
        dwpc = []
        thte = []
        tmwc = []
        line_val = data[start].strip().split()
        line = start
        line_length = len(line_val)
        while line_length != 0:
            line_val = data[line].strip().split()
            next_line_val = data[line+1]
            idx = next_line_val.find('STN YYMMDD/HHMM PMSL PRES SKTC STC1 SNFL WTNS\n')
            if idx != -1: break

            line_length = len(line_val)
            if line_length > 2:
                pres.append(float(line_val[0]))
                tmpc.append(float(line_val[1]))
                tmwc.append(float(line_val[2]))
                thte.append(float(line_val[4]))
                dwpc.append(float(line_val[3]))
            elif line_length != 0:
                hght.append(float(line_val[1]))
            line += 1
        end = line - 2
        num_plevs = int(((end - start) / 2) + 1)
        sum = sum + num_plevs

        min_p = np.maximum(max_p, max(pres))
        print("================== %s ==================" % (str(i)))
        print(min(pres), max(pres))
        #print(min(tmpc), max(tmpc), len(tmpc))
        #print(min(dwpc), max(dwpc), len(dwpc))
        #print(min(thte), max(thte), len(thte))
        #print(min(tmwc), max(tmwc), len(tmwc))
        print(strictly_decreasing(pres))
        print(strictly_increasing(hght))
        plt.plot(pres)
    #plt.ylim(0,500)
    plt.show()
    print(min_p)
    idx = data.index('STN YYMMDD/HHMM PMSL PRES SKTC STC1 SNFL WTNS\n')
    num_plevs = int((idx-2-indices[i+1])/2)
    sum = sum + num_plevs
else:
    ret = True
if sum < 201:
    ret = True

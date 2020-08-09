# --------------------------------------------------------------------------------------
#
# Standard header and other BUFKIT-related file data that we don't need floating around
# in the main codes.
#
# --------------------------------------------------------------------------------------
FILE_HEADER = [
    "SNPARM = PRES;TMPC;TMWC;DWPC;THTE;DRCT;SKNT;OMEG;CFRL;HGHT\r\n",
    "STNPRM = SHOW;LIFT;SWET;KINX;LCLP;PWAT;TOTL;CAPE;LCLT;CINS;EQLV;LFCT;BRCH\r\n"
]

SFC_HEADER = [
    "STN YYMMDD/HHMM PMSL PRES SKTC STC1 SNFL WTNS\r\n",
    "P01M C01M STC2 LCLD MCLD HCLD\r\n",
    "SNRA UWND VWND R01M BFGR T2MS\r\n",
    "Q2MS WXTS WXTP WXTZ WXTR USTM\r\n",
    "VSTM HLCY SLLH WSYM CDBP VSBK\r\n",
    "TD2M\r\n"
]

SFC_VALS = [
    "%s 0.00 0.00 0.00 0.00 0.00\r\n" % (0.00),
    "-9999.00 0.00 0.00 -9999.00 -9999.00 0.00\r\n",
    "0.00 0.00 0.00 0.00 0.00 0.00\r\n",
    "0.00 0.00 -0.00 999.00 -9999.00 0.00\r\n",
    "0.00\r\n"
]

# In the cases where we need additional "faked" data to avoid crashing BUFKIT. I can't
# imagine a scenario where a jet provides data above 16,000 m. This data is only used
# if RAP data is not supplied to augment the vertical ACARS profile.
EXTRA_LINES = [
    "107.80 -70.06 -70.14 -102.17 383.53 117.26 7.22 -9999.0\r\n",
    "0.00 16244.36\r\n",
    "94.90 -71.86 -71.93 -102.76 394.21 90.70 15.94 -9999.0\r\n",
    "0.00 16998.72\r\n",
    "83.20 -73.36 -73.43 -103.36 406.25 101.31 23.79 -9999.0\r\n",
    "0.00 17771.12\r\n",
    "72.90 -71.56 -71.66 -103.96 425.67 111.55 24.87 -9999.0\r\n",
    "0.00 18547.52\r\n",
    "63.70 -66.66 -66.87 -104.57 453.13 114.48 23.92 -9999.0\r\n",
    "0.00 19353.28\r\n",
    "55.50 -57.96 -58.63 -105.19 491.16 90.00 13.61 -9999.0\r\n",
    "0.00 20294.94\r\n",
    "48.20 -57.16 -57.99 -105.81 513.23 96.17 14.47 -9999.0\r\n",
    "0.00 21184.93\r\n",
    "41.80 -55.06 -56.25 -106.44 539.73 93.18 14.02 -9999.0\r\n",
    "0.00 22090.05\r\n",
    "36.10 -54.16 -55.64 -107.07 565.11 71.82 13.71 -9999.0\r\n",
    "0.00 23027.91\r\n",
    "31.20 -52.16 -54.19 -107.70 594.51 74.51 20.37 -9999.0\r\n",
    "0.00 23967.30\r\n",
    "26.80 -51.06 -53.60 -108.35 623.96 88.51 22.36 -9999.0\r\n",
    "0.00 24953.14\r\n",
    "22.90 -50.16 -53.25 -109.01 655.25 91.85 24.12 -9999.0\r\n",
    "0.00 25977.62\r\n",
    "19.50 -48.06 -52.17 -109.69 692.47 85.66 28.27 -9999.0\r\n",
    "0.00 27031.67\r\n",
    "16.50 -45.26 -50.89 -110.38 735.32 87.32 33.28 -9999.0\r\n",
    "0.00 28139.24\r\n"
]

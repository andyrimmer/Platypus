from matplotlib import pyplot
import sys

fieldToPlot = sys.argv[1]
passFieldVals = []
failFieldVals = []
nBins = 20

if len(sys.argv) > 2:
    nBins = int(sys.argv[2])

for line  in sys.stdin:
    if line.startswith("#"):
        continue
    else:
        cols = line.strip().split("\t")
        info = cols[7]
        theFilter = cols[6]

        for infoField in info.split(";"):
            key,vals = infoField.split("=")
            vals = vals.split(",")

            if key == fieldToPlot:
                for val in vals:
                    if theFilter == "PASS":
                        passFieldVals.append(float(val))
                    else:
                        failFieldVals.append(float(val))

pyplot.hist(passFieldVals, bins=nBins, label="PASS_" + fieldToPlot, normed=True)
pyplot.hist(failFieldVals, bins=nBins, label="FAIL_" + fieldToPlot, normed=True, alpha=0.5)
pyplot.legend()
pyplot.show()

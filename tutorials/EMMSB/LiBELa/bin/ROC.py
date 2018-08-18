#!/usr/bin/python
import sys, math, os.path, getopt


if (len(sys.argv) < 2):
        print ("Usage: ROC.py  -i input -a number_of_actives -d number_of_decoys -o output_file  [-v]");
        sys.exit(2);
try:
        opts, args = getopt.getopt(sys.argv[1:], "vi:a:d:o:")
except getopt.GetoptError, err:
        print str(err)
	print ("Usage: ROC.py  -i input -a number_of_actives -d number_of_decoys -o output_file [-v]");
        sys.exit(2);

verbose = False
actives = 0
decoys = 0
ligand = None
infile = None
output = "roc.dat";
for o,a in opts:
        if o == "-v":
            verbose = True
        elif o == "-a":
            actives = int(a)
        elif o == "-d":
            decoys = int(a)
	elif o == "-i":
	    infile = a
	elif o == "-o":
	    output = a;
        else:
            assert False, "unhandled option"

print "ROCcing data for %s : %d ligands and %d decoys ..." % (infile, actives, decoys)

bychance = (actives*1.0/(actives+decoys))

#print "Active ratio by chance: %7.4f" % (bychance)
# changes the previous
#os.system("more %s | sort -n -k 9 | awk '{print $1}' > tmp.dat " % (infile));
os.system("more %s | sort -g -k 9 | awk '{print $1}' > tmp.dat " % (infile));
output=open(output,'w');
output.write("# ROC for file %s\n" % (infile));
output.write("# NLIG: %d\n" % (actives))
output.write("# NDEC: %d\n" % (decoys))
data=open('tmp.dat', 'r');
count=0;
dec_found=0;
lig_found=0;
x=[];
y=[];
for line in data:
	count=count+1;
	mol = int(line)
    	if (mol <= actives):
		lig_found=lig_found+1
	else:
		dec_found = dec_found+1
		x.append((dec_found*1.0)/decoys);
		y.append((lig_found*1.0)/actives);
		output.write("%7.4f %7.4f %7.4f\n" % ((dec_found*1.0)/decoys, (lig_found*1.0)/actives, ((lig_found*1.0)/count)/(bychance)));
data.close()
os.system("rm tmp.dat")

# Computing logAUC

bi=0.00;
auc=0.00;
lauc=0.00;
for i in range(0, len(x)-1):
	auc = auc + ((x[i+1]-x[i])*y[i]); 				# Usual AUC
	bi = y[i+1] - x[i+1]*((y[i+1]-y[i])/(x[i+1]-x[i])); 		# Adjusted logAUC
	if (x[i] >= 0.001):
		lauc = lauc + ( ((y[i+1]-y[i])/math.log(10)) + bi*(math.log10(x[i+1]) - math.log10(x[i])));
lauc = lauc / (math.log10(1.0/0.001));

output.write("# AUC = %7.3f%s\n" % ( (auc*100), "%"));
output.write("# logAUC = %7.3f%s\n" % ((lauc*100.0), "%") );
output.write("# Adjusted logAUC = %7.3f%s\n" % ( ((lauc*100)-14.5), "%"));
output.close();

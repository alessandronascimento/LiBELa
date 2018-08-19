#!/usr/bin/env python

from random import randrange

mylist=[];
total=2156
desired=300;

for i in range(1, (total+1)):
    mol='../mol2/fda.{0}.mol2.gz'.format(i)
    mylist.append(mol)


nclean=total-desired;
for i in range(1, nclean):
    rindex=randrange(0, len(mylist))
    mylist.pop(rindex)

output=open('multimol.dat','w');
for i in range(0, len(mylist)):
    output.write("%s\n" % mylist[i]);

output.close();


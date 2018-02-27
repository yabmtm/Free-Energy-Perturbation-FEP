import os, sys, string
import numpy as np

usage = """Usage: mklambdas.sh run.mdp topol.top conf.gro itpfile.itp ....
    
    Makes directories with altered mdp file, and copied topology, gro, itp files etc."""

if len(sys.argv) < 3:
    print usage
    sys.exit(1)

mdpfile = sys.argv[1]
files_to_be_copied = sys.argv[2:]

lambdas = np.arange(0.0, 1.05, 0.050)
#vdw_lambdas = lambdas.tolist()[0:-1] + [1.0 for lam in lambdas]
#coul_lambdas =  [0.0 for lam in lambdas[0:-1]] + lambdas.tolist() 

#print 'vdw_lambdas', vdw_lambdas
#print 'coul_lambdas', coul_lambdas

# read in the mdp file contents
fin = open(mdpfile, 'r')
mdp_text = fin.read()
fin.close()

### NOT NEEDED FOR ALCHEMICAL FEP ###
## add list of vdw and coul lambdas
#mdp_text = mdp_text.replace( '$VDW-LAMBDAS$', string.joinfields(['%2.2f'%i for i in vdw_lambdas], ' ') )
#mdp_text = mdp_text.replace( '$COUL-LAMBDAS$', string.joinfields(['%2.2f'%i for i in coul_lambdas], ' ') )


for i in range(len(lambdas)):

    newdir="lambda_%d"%i
    print  "Making directory %s, and populating it"%newdir
    if not os.path.exists(newdir):
        os.mkdir(newdir)

    # add initial lambda
    custom_mdp_text = mdp_text.replace('$LAMBDA$', str(i))

    # write the files
    fout = open( os.path.join(newdir, 'grompp.mdp'), 'w')
    fout.write(custom_mdp_text)
    fout.close()

    for filename in files_to_be_copied:
        os.system('cp %s %s'%(filename, newdir))


